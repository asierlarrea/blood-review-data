#!/usr/bin/env Rscript

# Load utilities and set up paths
source(file.path(dirname(dirname(getwd())), "scripts", "utilities", "load_packages.R"))
ensure_output_dirs()


# Biological ID Mapping Functions
# Converts various protein/gene accessions to HGNC gene symbols
# Author: Data Analysis Pipeline
# Date: 2024

# Load required libraries for biological ID mapping
suppressMessages({
  # Try to load biomaRt for Ensembl queries
  if (!require("biomaRt", quietly = TRUE)) {
    message("biomaRt not available - using alternative mapping methods")
  }
  
  # Try to load org.Hs.eg.db for human gene annotations
  if (!require("org.Hs.eg.db", quietly = TRUE)) {
    message("org.Hs.eg.db not available - using alternative mapping methods")
  }
})

# Function to identify accession type
identify_accession_type <- function(accession) {
  acc <- as.character(accession)
  
  if (grepl("^ENSP[0-9]+", acc)) {
    return("ENSEMBL_PROTEIN")
  } else if (grepl("^[A-Z][0-9][A-Z0-9]{3}[0-9]$", acc)) {
    return("UNIPROT_SWISSPROT")
  } else if (grepl("^[A-Z][0-9][A-Z0-9]{3}[0-9]-[0-9]+$", acc)) {
    return("UNIPROT_ISOFORM")
  } else if (grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", acc)) {
    return("UNIPROT_TREMBL")
  } else if (grepl("^[A-Z]{1,2}_[0-9]+", acc)) {
    return("REFSEQ_PROTEIN")
  } else if (grepl("^[A-Z]{2,}[0-9]*$", acc) && nchar(acc) >= 2 && nchar(acc) <= 15) {
    return("GENE_SYMBOL")
  } else {
    return("UNKNOWN")
  }
}

# Function to map ENSP IDs to gene symbols using biomaRt
map_ensp_to_gene <- function(ensp_ids) {
  if (!require("biomaRt", quietly = TRUE)) {
    message("biomaRt not available for ENSP mapping")
    return(rep(NA, length(ensp_ids)))
  }
  
  tryCatch({
    # Connect to Ensembl
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Query for gene symbols
    results <- getBM(
      attributes = c("ensembl_peptide_id", "hgnc_symbol"),
      filters = "ensembl_peptide_id",
      values = ensp_ids,
      mart = ensembl
    )
    
    # Map back to original order
    gene_symbols <- character(length(ensp_ids))
    for (i in seq_along(ensp_ids)) {
      match_idx <- which(results$ensembl_peptide_id == ensp_ids[i])
      if (length(match_idx) > 0 && results$hgnc_symbol[match_idx[1]] != "") {
        gene_symbols[i] <- results$hgnc_symbol[match_idx[1]]
      } else {
        gene_symbols[i] <- NA
      }
    }
    
    return(gene_symbols)
  }, error = function(e) {
    message(paste("Error in ENSP mapping:", e$message))
    return(rep(NA, length(ensp_ids)))
  })
}

# Function to map UniProt IDs to gene symbols using biomaRt
map_uniprot_to_gene <- function(uniprot_ids) {
  if (!require("biomaRt", quietly = TRUE)) {
    message("biomaRt not available for UniProt mapping")
    return(rep(NA, length(uniprot_ids)))
  }
  
  tryCatch({
    # Connect to Ensembl
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Clean UniProt IDs (remove isoform suffixes)
    clean_ids <- gsub("-[0-9]+$", "", uniprot_ids)
    
    # Query for gene symbols
    results <- getBM(
      attributes = c("uniprotswissprot", "hgnc_symbol"),
      filters = "uniprotswissprot",
      values = clean_ids,
      mart = ensembl
    )
    
    # Map back to original order
    gene_symbols <- character(length(uniprot_ids))
    for (i in seq_along(uniprot_ids)) {
      clean_id <- clean_ids[i]
      match_idx <- which(results$uniprotswissprot == clean_id)
      if (length(match_idx) > 0 && results$hgnc_symbol[match_idx[1]] != "") {
        gene_symbols[i] <- results$hgnc_symbol[match_idx[1]]
      } else {
        gene_symbols[i] <- NA
      }
    }
    
    return(gene_symbols)
  }, error = function(e) {
    message(paste("Error in UniProt mapping:", e$message))
    return(rep(NA, length(uniprot_ids)))
  })
}

# Fallback function using org.Hs.eg.db
map_using_orgdb <- function(accessions, from_type) {
  if (!require("org.Hs.eg.db", quietly = TRUE)) {
    return(rep(NA, length(accessions)))
  }
  
  tryCatch({
    if (from_type == "UNIPROT_SWISSPROT") {
      # Map UniProt to gene symbols
      results <- mapIds(org.Hs.eg.db, 
                       keys = gsub("-[0-9]+$", "", accessions),
                       column = "SYMBOL",
                       keytype = "UNIPROT",
                       multiVals = "first")
      return(as.character(results))
    } else {
      return(rep(NA, length(accessions)))
    }
  }, error = function(e) {
    message(paste("Error in org.Hs.eg.db mapping:", e$message))
    return(rep(NA, length(accessions)))
  })
}

# Main mapping function
map_to_gene_symbol <- function(accessions, descriptions = NULL) {
  if (length(accessions) == 0) {
    return(character(0))
  }
  
  gene_symbols <- character(length(accessions))
  
  for (i in seq_along(accessions)) {
    acc <- as.character(accessions[i])
    desc <- if (!is.null(descriptions)) as.character(descriptions[i]) else NULL
    
    # Skip empty accessions
    if (is.na(acc) || acc == "") {
      gene_symbols[i] <- NA
      next
    }
    
    # Identify accession type
    acc_type <- identify_accession_type(acc)
    
    # Apply appropriate mapping
    if (acc_type == "GENE_SYMBOL") {
      # Already a gene symbol
      gene_symbols[i] <- toupper(acc)
    } else if (acc_type == "ENSEMBL_PROTEIN") {
      # Map ENSP to gene symbol
      mapped <- map_ensp_to_gene(acc)
      gene_symbols[i] <- if (!is.na(mapped[1])) toupper(mapped[1]) else NA
    } else if (acc_type %in% c("UNIPROT_SWISSPROT", "UNIPROT_TREMBL", "UNIPROT_ISOFORM")) {
      # Map UniProt to gene symbol
      mapped <- map_uniprot_to_gene(acc)
      if (is.na(mapped[1])) {
        # Try org.Hs.eg.db as fallback
        mapped <- map_using_orgdb(acc, acc_type)
      }
      gene_symbols[i] <- if (!is.na(mapped[1])) toupper(mapped[1]) else NA
    } else {
      # Try to extract from description if available
      if (!is.null(desc) && !is.na(desc) && desc != "") {
        extracted <- extract_gene_from_description(desc)
        gene_symbols[i] <- if (!is.na(extracted)) toupper(extracted) else NA
      } else {
        gene_symbols[i] <- NA
      }
    }
  }
  
  return(gene_symbols)
}

# Helper function to extract gene symbols from descriptions (as fallback)
extract_gene_from_description <- function(description) {
  if (is.na(description) || description == "") {
    return(NA)
  }
  
  desc <- as.character(description)
  
  # Pattern 1: GN=GENE_NAME (UniProt format)
  if (grepl("GN=([A-Za-z0-9_-]+)", desc)) {
    gene <- gsub(".*GN=([A-Za-z0-9_-]+).*", "\\1", desc)
    if (nchar(gene) >= 2 && nchar(gene) <= 15) {
      return(toupper(gene))
    }
  }
  
  # Pattern 2: \GName=GENE_NAME (PeptideAtlas format)
  if (grepl("\\\\GName=([A-Za-z0-9_-]+)", desc)) {
    gene <- gsub(".*\\\\GName=([A-Za-z0-9_-]+).*", "\\1", desc)
    if (nchar(gene) >= 2 && nchar(gene) <= 15) {
      return(toupper(gene))
    }
  }
  
  # Pattern 3: Gene symbols in brackets or after keywords
  patterns <- c(
    "\\[([A-Z][A-Z0-9_-]{1,14})\\]",  # [GENE]
    "Gene: ([A-Z][A-Z0-9_-]{1,14})",   # Gene: GENE
    "Symbol: ([A-Z][A-Z0-9_-]{1,14})", # Symbol: GENE
    "^([A-Z][A-Z0-9_-]{1,14})\\s+"     # GENE at start
  )
  
  for (pattern in patterns) {
    if (grepl(pattern, desc, ignore.case = TRUE)) {
      gene <- gsub(pattern, "\\1", desc, ignore.case = TRUE)
      if (nchar(gene) >= 2 && nchar(gene) <= 15 && 
          !grepl("PROTEIN|FRAGMENT|DOMAIN|FAMILY", gene, ignore.case = TRUE)) {
        return(toupper(gene))
      }
    }
  }
  
  return(NA)
}

# Function to test the mapping system
test_id_mapping <- function() {
  message("Testing ID mapping system...")
  
  # Test cases
  test_ids <- c(
    "ENSP00000295897",  # Ensembl protein
    "P02768",           # UniProt (ALB)
    "Q13438",           # UniProt (OS9)
    "ALB",              # Gene symbol
    "INVALID_ID"        # Invalid
  )
  
  results <- map_to_gene_symbol(test_ids)
  
  for (i in seq_along(test_ids)) {
    message(paste(test_ids[i], "->", 
                  if (is.na(results[i])) "NOT_MAPPED" else results[i]))
  }
  
  return(results)
}

# Function to create a local cache for frequent mappings
create_mapping_cache <- function() {
  # This could be expanded to include common protein->gene mappings
  # to reduce API calls to biomaRt
  cache <- list()
  return(cache)
}

message("ID mapping functions loaded. Key functions:")
message("- map_to_gene_symbol(accessions, descriptions = NULL)")
message("- identify_accession_type(accession)")
message("- test_id_mapping()") 
