#!/usr/bin/env Rscript

# Simple ID Mapping: UniProt/ENSP -> Gene Names with Lazy Loading
# Purpose: Convert protein accessions to gene symbols using cache-first, database-fallback approach

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(httr)
  library(jsonlite)
})

# Global variable to store loaded mappings (cache)
.mapping_cache <- NULL

# Load protein to gene mappings from CSV cache
load_mapping_cache <- function(cache_file = "data/cache/protein_to_gene_mappings.csv") {
  if (is.null(.mapping_cache)) {
    if (file.exists(cache_file)) {
      message(sprintf("Loading protein mappings from %s...", cache_file))
      .mapping_cache <<- read_csv(cache_file, show_col_types = FALSE)
      message(sprintf("Loaded %d cached protein → gene mappings", nrow(.mapping_cache)))
    } else {
      # Create empty cache if file doesn't exist
      message("No cache file found, creating new cache...")
      dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
      .mapping_cache <<- data.frame(
        protein_id = character(0),
        gene_symbol = character(0),
        source = character(0),
        description = character(0),
        query_date = character(0),
        stringsAsFactors = FALSE
      )
      write_csv(.mapping_cache, cache_file)
    }
  }
  return(.mapping_cache)
}

# Query UniProt for protein ID to gene mapping
query_uniprot <- function(protein_id) {
  tryCatch({
    # UniProt REST API
    url <- sprintf("https://rest.uniprot.org/uniprotkb/%s.json", protein_id)
    response <- GET(url)
    
    if (status_code(response) == 200) {
      data <- fromJSON(content(response, "text", encoding = "UTF-8"))
      
      # Extract gene name - correct nested structure
      if (!is.null(data$genes) && nrow(data$genes) > 0) {
        gene_symbol <- data$genes$geneName$value[1]  # First gene name from nested structure
        
        # Extract protein description
        description <- ""
        if (!is.null(data$proteinDescription$recommendedName$fullName$value)) {
          description <- data$proteinDescription$recommendedName$fullName$value
        }
        
        return(list(
          gene_symbol = toupper(gene_symbol),
          description = description,
          source = "UniProt_API",
          query_date = Sys.Date()
        ))
      }
    }
    return(NULL)
  }, error = function(e) {
    message(sprintf("UniProt query error for %s: %s", protein_id, e$message))
    return(NULL)
  })
}

# Query Ensembl for ENSP ID to gene mapping
query_ensembl <- function(ensp_id) {
  tryCatch({
    # Step 1: Get protein info to find parent transcript
    protein_url <- sprintf("https://rest.ensembl.org/lookup/id/%s?content-type=application/json", ensp_id)
    protein_response <- GET(protein_url)
    
    if (status_code(protein_response) == 200) {
      protein_data <- fromJSON(content(protein_response, "text", encoding = "UTF-8"))
      
      # Get the parent transcript ID
      if (!is.null(protein_data$Parent)) {
        transcript_id <- protein_data$Parent
        
        # Step 2: Get transcript info to find parent gene
        transcript_url <- sprintf("https://rest.ensembl.org/lookup/id/%s?content-type=application/json", transcript_id)
        transcript_response <- GET(transcript_url)
        
        if (status_code(transcript_response) == 200) {
          transcript_data <- fromJSON(content(transcript_response, "text", encoding = "UTF-8"))
          
          # Get the parent gene ID
          if (!is.null(transcript_data$Parent)) {
            gene_id <- transcript_data$Parent
            
            # Step 3: Get gene info for the clean gene symbol
            gene_url <- sprintf("https://rest.ensembl.org/lookup/id/%s?content-type=application/json", gene_id)
            gene_response <- GET(gene_url)
            
            if (status_code(gene_response) == 200) {
              gene_data <- fromJSON(content(gene_response, "text", encoding = "UTF-8"))
              
              if (!is.null(gene_data$display_name)) {
                gene_symbol <- gene_data$display_name
                
                # Extract description
                description <- ""
                if (!is.null(gene_data$description)) {
                  description <- gene_data$description
                }
                
                return(list(
                  gene_symbol = toupper(gene_symbol),
                  description = description,
                  source = "Ensembl_API",
                  query_date = Sys.Date()
                ))
              }
            }
          }
        }
      }
    }
    return(NULL)
  }, error = function(e) {
    message(sprintf("Ensembl query error for %s: %s", ensp_id, e$message))
    return(NULL)
  })
}

# Add new mapping to cache and save
add_to_cache <- function(protein_id, gene_symbol, source, description, query_date, 
                        cache_file = "data/cache/protein_to_gene_mappings.csv") {
  # Load current cache
  current_cache <- load_mapping_cache(cache_file)
  
  # Add new mapping
  new_row <- data.frame(
    protein_id = protein_id,
    gene_symbol = gene_symbol,
    source = source,
    description = description,
    query_date = as.character(query_date),
    stringsAsFactors = FALSE
  )
  
  # Update cache in memory and file
  .mapping_cache <<- rbind(current_cache, new_row)
  write_csv(.mapping_cache, cache_file)
  
  message(sprintf("Cached: %s → %s (from %s)", protein_id, gene_symbol, source))
}

# Main function: convert any ID to gene symbol with lazy loading
convert_to_gene_symbol <- function(ids, cache_file = "data/cache/protein_to_gene_mappings.csv") {
  if (length(ids) == 0) return(character(0))
  
  message(sprintf("Converting %d IDs to gene symbols...", length(ids)))
  
  # Load mapping cache
  mappings <- load_mapping_cache(cache_file)
  
  results <- rep(NA_character_, length(ids))
  cache_hits <- 0
  db_queries <- 0
  
  for (i in seq_along(ids)) {
    id <- as.character(ids[i])
    
    # Skip if already NA or empty
    if (is.na(id) || id == "") {
      next
    }
    
    # Clean the ID (remove spaces, taxonomy prefixes)
    clean_id <- clean_protein_id(id)
    
    # Already a gene symbol? Keep it
    if (is_gene_symbol(clean_id)) {
      results[i] <- toupper(clean_id)
      next
    }
    
    # Step 1: Check cache first
    mapping_row <- mappings[mappings$protein_id == clean_id, ]
    if (nrow(mapping_row) > 0) {
      results[i] <- mapping_row$gene_symbol[1]
      cache_hits <- cache_hits + 1
      next
    }
    
    # Step 2: Not in cache - query appropriate database
    db_result <- NULL
    
    if (grepl("^[A-Z][0-9][A-Z0-9]{3}[0-9]|^[OPQ][0-9][A-Z0-9]{3}[0-9]", clean_id)) {
      # UniProt accession
      message(sprintf("Querying UniProt for %s...", clean_id))
      db_result <- query_uniprot(clean_id)
      db_queries <- db_queries + 1
      
    } else if (grepl("^ENSP[0-9]+", clean_id)) {
      # Ensembl protein accession
      message(sprintf("Querying Ensembl for %s...", clean_id))
      db_result <- query_ensembl(clean_id)
      db_queries <- db_queries + 1
    }
    
    # Step 3: If found, cache it and use it
    if (!is.null(db_result)) {
      results[i] <- db_result$gene_symbol
      add_to_cache(clean_id, db_result$gene_symbol, db_result$source, 
                  db_result$description, db_result$query_date, cache_file)
    }
    
    # Small delay to be nice to APIs
    if (db_queries %% 10 == 0) {
      Sys.sleep(1)  # 1 second pause every 10 queries
    }
  }
  
  # Count successes
  mapped_count <- sum(!is.na(results))
  message(sprintf("Successfully mapped: %d out of %d IDs (%.1f%%)", 
                 mapped_count, length(ids), 100 * mapped_count / length(ids)))
  message(sprintf("Cache hits: %d, Database queries: %d", cache_hits, db_queries))
  
  return(results)
}

# Clean protein ID (remove prefixes, suffixes)
clean_protein_id <- function(id) {
  # Remove taxonomy prefix (e.g., "9606.ENSP00000295897" -> "ENSP00000295897")
  clean_id <- gsub("^[0-9]+\\.", "", id)
  
  # Remove isoform suffix (e.g., "P02768-1" -> "P02768")
  clean_id <- gsub("-[0-9]+$", "", clean_id)
  
  # Remove version suffix (e.g., "ENSP00000295897.4" -> "ENSP00000295897")
  clean_id <- gsub("\\.[0-9]+$", "", clean_id)
  
  return(trimws(clean_id))
}

# Check if ID is already a gene symbol
is_gene_symbol <- function(id) {
  if (is.na(id) || id == "") return(FALSE)
  
  # Gene symbols: 2-15 characters, start with letter, alphanumeric + underscore/dash
  if (nchar(id) < 2 || nchar(id) > 15) return(FALSE)
  if (!grepl("^[A-Za-z]", id)) return(FALSE)
  if (!grepl("^[A-Za-z][A-Za-z0-9_-]*$", id)) return(FALSE)
  
  # Exclude known protein ID patterns and invalid patterns
  if (grepl("^ENS|^[OPQ][0-9]|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|INVALID", id)) return(FALSE)
  
  return(TRUE)
}

# Function to get mapping statistics
get_mapping_stats <- function(cache_file = "data/cache/protein_to_gene_mappings.csv") {
  mappings <- load_mapping_cache(cache_file)
  
  if (nrow(mappings) == 0) {
    message("No mappings in cache yet")
    return(data.frame())
  }
  
  stats <- mappings %>%
    group_by(source) %>%
    summarise(
      count = n(),
      unique_genes = length(unique(gene_symbol)),
      .groups = "drop"
    ) %>%
    arrange(desc(count))
  
  return(stats)
}

# Test function
test_simple_mapping <- function() {
  message("=== TESTING LAZY-LOADING ID MAPPING ===")
  
  # Test cases - mix of easy and challenging
  test_ids <- c(
    # Gene symbols (should stay as-is)
    "ALB", "A2M",
    # Common UniProt (might be cached)
    "P02768",  # ALB
    # Less common UniProt (probably not cached)
    "P62166",  # NCS1
    # Common ENSP (might be cached)
    "ENSP00000295897",  # ALB
    # Less common ENSP (probably not cached)
    "ENSP00000323929",  # A2M
    # PAXDB format
    "9606.ENSP00000245907",  # C3
    # Invalid
    "INVALID123"
  )
  
  results <- convert_to_gene_symbol(test_ids)
  
  cat("\nTest Results:\n")
  cat("=============\n")
  for (i in seq_along(test_ids)) {
    cat(sprintf("%s -> %s\n", test_ids[i], 
                ifelse(is.na(results[i]), "NOT_MAPPED", results[i])))
  }
  
  # Show cache statistics
  cat("\nCache Statistics:\n")
  stats <- get_mapping_stats()
  if (nrow(stats) > 0) {
    print(stats)
  }
  
  return(results)
}

# Clear cache function (for testing)
clear_cache <- function(cache_file = "data/cache/protein_to_gene_mappings.csv") {
  if (file.exists(cache_file)) {
    file.remove(cache_file)
  }
  .mapping_cache <<- NULL
  message("Cache cleared")
}

# Print usage message
message("Lazy-loading ID mapping loaded!")
message("Usage: gene_symbols <- convert_to_gene_symbol(your_protein_ids)")
message("Test with: test_simple_mapping()")
message("Clear cache: clear_cache()")
message("Note: First run will be slow (database queries), subsequent runs will be fast (cached)") 