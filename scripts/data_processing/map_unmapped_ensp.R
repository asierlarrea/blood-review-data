#!/usr/bin/env Rscript

# Load utilities and set up paths
source(file.path(dirname(dirname(getwd())), "scripts", "utilities", "load_packages.R"))
ensure_output_dirs()


# Map Unmapped ENSP IDs using biomaRt
# Author: Data Analysis Pipeline
# Date: 2024

library(biomaRt)

# Load the fast mapping functions to check which IDs are already mapped
source("fast_id_mapping.R")

# Function to get unmapped ENSP IDs from PaxDB
get_unmapped_ensp_ids <- function() {
  message("Reading PaxDB data...")
  
  # Read PaxDB data
  paxdb_data <- read.csv(get_data_path("paxdb_plasma.csv"), stringsAsFactors = FALSE)
  
  # Extract valid ENSP IDs
  valid_rows <- !is.na(paxdb_data$string_external_id) & 
               !is.na(paxdb_data$abundance) & 
               paxdb_data$abundance > 0 &
               grepl("9606\\.ENSP[0-9]+", paxdb_data$string_external_id)
  
  # Extract clean ENSP IDs (remove organism prefix)
  ensp_ids <- gsub("9606\\.", "", paxdb_data$string_external_id[valid_rows])
  abundances <- paxdb_data$abundance[valid_rows]
  
  message(paste("Found", length(ensp_ids), "valid ENSP IDs"))
  
  # Test which ones can be mapped locally
  local_mapping <- fast_map_to_gene_symbol(ensp_ids)
  
  # Find unmapped IDs
  unmapped_mask <- is.na(local_mapping) | local_mapping == ""
  unmapped_ensp <- ensp_ids[unmapped_mask]
  unmapped_abundances <- abundances[unmapped_mask]
  
  message(paste("Local mapping:", sum(!unmapped_mask), "mapped,", sum(unmapped_mask), "unmapped"))
  
  # Create data frame with unmapped IDs and their abundances
  unmapped_data <- data.frame(
    ensp_id = unmapped_ensp,
    abundance = unmapped_abundances,
    stringsAsFactors = FALSE
  )
  
  # Sort by abundance (highest first)
  unmapped_data <- unmapped_data[order(-unmapped_data$abundance), ]
  
  return(unmapped_data)
}

# Function to map ENSP IDs using biomaRt in batches
map_ensp_with_biomart <- function(ensp_ids, batch_size = 500) {
  message(paste("Mapping", length(ensp_ids), "ENSP IDs using biomaRt..."))
  
  # Initialize results
  results <- data.frame(
    ensp_id = ensp_ids,
    gene_symbol = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Try to connect to biomaRt
  tryCatch({
    message("Connecting to Ensembl biomaRt...")
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
    message("Connected successfully!")
    
    # Process in batches to avoid timeouts
    for(i in seq(1, length(ensp_ids), batch_size)) {
      end_idx <- min(i + batch_size - 1, length(ensp_ids))
      batch_ids <- ensp_ids[i:end_idx]
      
      message(paste("Processing batch", ceiling(i/batch_size), "- IDs", i, "to", end_idx))
      
      tryCatch({
        # Query biomaRt
        batch_results <- getBM(
          attributes = c("ensembl_peptide_id", "hgnc_symbol"),
          filters = "ensembl_peptide_id", 
          values = batch_ids,
          mart = ensembl
        )
        
        # Map results back to original positions
        for(j in 1:length(batch_ids)) {
          match_idx <- which(batch_results$ensembl_peptide_id == batch_ids[j])
          if(length(match_idx) > 0 && batch_results$hgnc_symbol[match_idx[1]] != "") {
            original_idx <- i + j - 1
            results$gene_symbol[original_idx] <- toupper(batch_results$hgnc_symbol[match_idx[1]])
          }
        }
        
        # Progress update
        mapped_so_far <- sum(!is.na(results$gene_symbol))
        message(paste("  Batch completed. Total mapped so far:", mapped_so_far))
        
        # Brief pause to be nice to the server
        Sys.sleep(1)
        
      }, error = function(e) {
        message(paste("  Error in batch", ceiling(i/batch_size), ":", e$message))
      })
    }
    
  }, error = function(e) {
    message(paste("biomaRt connection error:", e$message))
  })
  
  # Summary
  total_mapped <- sum(!is.na(results$gene_symbol))
  mapping_rate <- round(total_mapped / length(ensp_ids) * 100, 1)
  message(paste("biomaRt mapping completed:", total_mapped, "out of", length(ensp_ids), "mapped (", mapping_rate, "%)"))
  
  return(results)
}

# Function to save biomaRt results for future use
save_biomart_mappings <- function(biomart_results, unmapped_data) {
  # Combine with abundance information
  combined_results <- merge(biomart_results, unmapped_data, by = "ensp_id", all.x = TRUE)
  
  # Filter successful mappings
  successful_mappings <- combined_results[!is.na(combined_results$gene_symbol), ]
  
  if(nrow(successful_mappings) > 0) {
    # Sort by abundance
    successful_mappings <- successful_mappings[order(-successful_mappings$abundance), ]
    
    # Save to CSV
    write.csv(successful_mappings, get_output_path("biomart_ensp_mappings.csv", "tables"), row.names = FALSE)
    message(paste("Saved", nrow(successful_mappings), "successful mappings to biomart_ensp_mappings.csv"))
    
    # Create R code for adding to fast_id_mapping.R
    r_code_lines <- character()
    r_code_lines <- c(r_code_lines, "    # Additional ENSP mappings from biomaRt")
    
    for(i in 1:min(nrow(successful_mappings), 100)) {  # Limit to top 100
      ensp <- successful_mappings$ensp_id[i]
      gene <- successful_mappings$gene_symbol[i]
      abundance <- successful_mappings$abundance[i]
      r_code_lines <- c(r_code_lines, paste0('    "', ensp, '" = "', gene, '",  # abundance: ', abundance))
    }
    
    # Save R code to file
    writeLines(r_code_lines, "biomart_ensp_mappings.R")
    message("Generated R code in biomart_ensp_mappings.R for top 100 mappings")
    
    return(successful_mappings)
  } else {
    message("No successful mappings found")
    return(data.frame())
  }
}

# Main execution function
main <- function() {
  message("=== ENSP Mapping with biomaRt ===")
  
  # Step 1: Get unmapped ENSP IDs
  message("\nStep 1: Identifying unmapped ENSP IDs...")
  unmapped_data <- get_unmapped_ensp_ids()
  
  if(nrow(unmapped_data) == 0) {
    message("No unmapped ENSP IDs found!")
    return()
  }
  
  message(paste("Top 10 unmapped ENSP IDs by abundance:"))
  print(head(unmapped_data, 10))
  
  # Step 2: Map using biomaRt
  message("\nStep 2: Mapping with biomaRt...")
  biomart_results <- map_ensp_with_biomart(unmapped_data$ensp_id, batch_size = 200)
  
  # Step 3: Save results
  message("\nStep 3: Saving results...")
  successful_mappings <- save_biomart_mappings(biomart_results, unmapped_data)
  
  # Step 4: Summary
  message("\n=== SUMMARY ===")
  message(paste("Total ENSP IDs in PaxDB:", nrow(unmapped_data) + 20))  # +20 for locally mapped
  message(paste("Already mapped locally:", 20))
  message(paste("Attempted biomaRt mapping:", nrow(unmapped_data)))
  message(paste("Successfully mapped by biomaRt:", nrow(successful_mappings)))
  message(paste("Total mapped:", 20 + nrow(successful_mappings)))
  
  if(nrow(successful_mappings) > 0) {
    total_coverage <- round((20 + nrow(successful_mappings)) / (nrow(unmapped_data) + 20) * 100, 1)
    message(paste("Overall mapping coverage:", total_coverage, "%"))
    
    message("\nTop 10 newly mapped proteins:")
    print(head(successful_mappings[, c("gene_symbol", "abundance")], 10))
  }
}

# Run if called directly
if(!interactive()) {
  main()
} 
