#!/usr/bin/env Rscript

# Efficient Biological ID Mapping with Batching and Fallbacks
# Author: Data Analysis Pipeline
# Date: 2024

# Function for batch ID mapping with fallbacks
efficient_map_to_gene_symbol <- function(accessions, descriptions = NULL, batch_size = 500) {
  if (length(accessions) == 0) {
    return(character(0))
  }
  
  message(paste("Mapping", length(accessions), "IDs to gene symbols..."))
  gene_symbols <- character(length(accessions))
  
  # Group by accession type for efficient batch processing
  acc_types <- sapply(accessions, identify_accession_type_simple)
  
  # Process gene symbols first (no mapping needed)
  gene_mask <- acc_types == "GENE_SYMBOL"
  if(sum(gene_mask) > 0) {
    gene_symbols[gene_mask] <- toupper(accessions[gene_mask])
    message(paste("Direct gene symbols:", sum(gene_mask)))
  }
  
  # Process UniProt IDs in batches
  uniprot_mask <- acc_types %in% c("UNIPROT_SWISSPROT", "UNIPROT_TREMBL")
  if(sum(uniprot_mask) > 0) {
    uniprot_ids <- accessions[uniprot_mask]
    uniprot_results <- batch_map_uniprot(uniprot_ids, batch_size)
    gene_symbols[uniprot_mask] <- uniprot_results
    message(paste("UniProt mapped:", sum(!is.na(uniprot_results))))
  }
  
  # Process ENSP IDs in batches
  ensp_mask <- acc_types == "ENSEMBL_PROTEIN"
  if(sum(ensp_mask) > 0) {
    ensp_ids <- accessions[ensp_mask]
    ensp_results <- batch_map_ensp(ensp_ids, batch_size)
    gene_symbols[ensp_mask] <- ensp_results
    message(paste("ENSP mapped:", sum(!is.na(ensp_results))))
  }
  
  # For remaining unmapped IDs, try description parsing
  unmapped_mask <- is.na(gene_symbols) | gene_symbols == ""
  if(sum(unmapped_mask) > 0 && !is.null(descriptions)) {
    desc_results <- sapply(descriptions[unmapped_mask], extract_gene_from_description_simple)
    gene_symbols[unmapped_mask] <- desc_results
    message(paste("Description parsed:", sum(!is.na(desc_results))))
  }
  
  return(gene_symbols)
}

# Simple accession type identification
identify_accession_type_simple <- function(accession) {
  acc <- as.character(accession)
  
  if (grepl("^ENSP[0-9]+", acc)) {
    return("ENSEMBL_PROTEIN")
  } else if (grepl("^[A-Z][0-9][A-Z0-9]{3}[0-9]", acc)) {
    return("UNIPROT_SWISSPROT")
  } else if (grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]", acc)) {
    return("UNIPROT_TREMBL")
  } else if (grepl("^[A-Z]{2,}[0-9]*$", acc) && nchar(acc) >= 2 && nchar(acc) <= 15) {
    return("GENE_SYMBOL")
  } else {
    return("UNKNOWN")
  }
}

# Batch UniProt mapping with local fallback
batch_map_uniprot <- function(uniprot_ids, batch_size = 500) {
  results <- character(length(uniprot_ids))
  
  # Try biomaRt first
  if (require("biomaRt", quietly = TRUE)) {
    tryCatch({
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
      
      # Process in batches
      for(i in seq(1, length(uniprot_ids), batch_size)) {
        end_idx <- min(i + batch_size - 1, length(uniprot_ids))
        batch_ids <- gsub("-[0-9]+$", "", uniprot_ids[i:end_idx])  # Remove isoforms
        
        batch_results <- getBM(
          attributes = c("uniprotswissprot", "hgnc_symbol"),
          filters = "uniprotswissprot",
          values = batch_ids,
          mart = ensembl
        )
        
        # Map back to original positions
        for(j in 1:length(batch_ids)) {
          match_idx <- which(batch_results$uniprotswissprot == batch_ids[j])
          if(length(match_idx) > 0 && batch_results$hgnc_symbol[match_idx[1]] != "") {
            results[i + j - 1] <- toupper(batch_results$hgnc_symbol[match_idx[1]])
          }
        }
        
        if(i %% (batch_size * 5) == 1) {
          message(paste("Processed", min(end_idx, length(uniprot_ids)), "of", length(uniprot_ids), "UniProt IDs"))
        }
      }
    }, error = function(e) {
      message(paste("biomaRt error for UniProt:", e$message))
    })
  }
  
  # Use local mapping for common proteins as fallback
  unmapped <- is.na(results) | results == ""
  if(sum(unmapped) > 0) {
    results[unmapped] <- sapply(uniprot_ids[unmapped], map_common_uniprot)
  }
  
  return(results)
}

# Batch ENSP mapping
batch_map_ensp <- function(ensp_ids, batch_size = 500) {
  results <- character(length(ensp_ids))
  
  # Try biomaRt first
  if (require("biomaRt", quietly = TRUE)) {
    tryCatch({
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
      
      # Process in batches
      for(i in seq(1, length(ensp_ids), batch_size)) {
        end_idx <- min(i + batch_size - 1, length(ensp_ids))
        batch_ids <- ensp_ids[i:end_idx]
        
        batch_results <- getBM(
          attributes = c("ensembl_peptide_id", "hgnc_symbol"),
          filters = "ensembl_peptide_id",
          values = batch_ids,
          mart = ensembl
        )
        
        # Map back to original positions
        for(j in 1:length(batch_ids)) {
          match_idx <- which(batch_results$ensembl_peptide_id == batch_ids[j])
          if(length(match_idx) > 0 && batch_results$hgnc_symbol[match_idx[1]] != "") {
            results[i + j - 1] <- toupper(batch_results$hgnc_symbol[match_idx[1]])
          }
        }
        
        if(i %% (batch_size * 5) == 1) {
          message(paste("Processed", min(end_idx, length(ensp_ids)), "of", length(ensp_ids), "ENSP IDs"))
        }
      }
    }, error = function(e) {
      message(paste("biomaRt error for ENSP:", e$message))
    })
  }
  
  return(results)
}

# Common protein mappings (fallback for when biomaRt fails)
map_common_uniprot <- function(uniprot_id) {
  common_mappings <- list(
    "P02768" = "ALB",      # Albumin
    "P02787" = "TF",       # Transferrin
    "P00738" = "HP",       # Haptoglobin
    "P02765" = "AHSG",     # Alpha-2-HS-glycoprotein
    "P02647" = "APOA1",    # Apolipoprotein A-I
    "P04217" = "A1BG",     # Alpha-1B-glycoprotein
    "P01023" = "A2M",      # Alpha-2-macroglobulin
    "P02671" = "FGA",      # Fibrinogen alpha chain
    "P02675" = "FGB",      # Fibrinogen beta chain
    "P02679" = "FGG",      # Fibrinogen gamma chain
    "P00450" = "CP",       # Ceruloplasmin
    "P02766" = "TTR",      # Transthyretin
    "P02774" = "GC",       # Vitamin D-binding protein
    "P00747" = "PLG",      # Plasminogen
    "P01024" = "C3",       # Complement C3
    "P0C0L4" = "C4A",      # Complement C4-A
    "P0C0L5" = "C4B",      # Complement C4-B
    "P01031" = "C5",       # Complement C5
    "P13671" = "C6",       # Complement C6
    "P10643" = "C7",       # Complement C7
    "P07357" = "C8A",      # Complement C8 alpha
    "P07358" = "C8B",      # Complement C8 beta
    "P07360" = "C8G",      # Complement C8 gamma
    "P02748" = "C9"        # Complement C9
  )
  
  clean_id <- gsub("-[0-9]+$", "", uniprot_id)
  return(common_mappings[[clean_id]])
}

# Simple description parsing
extract_gene_from_description_simple <- function(description) {
  if (is.na(description) || description == "") {
    return(NA)
  }
  
  desc <- as.character(description)
  
  # Pattern 1: GN=GENE_NAME
  if (grepl("GN=([A-Za-z0-9_-]+)", desc)) {
    gene <- gsub(".*GN=([A-Za-z0-9_-]+).*", "\\1", desc)
    if (nchar(gene) >= 2 && nchar(gene) <= 15) {
      return(toupper(gene))
    }
  }
  
  # Pattern 2: \GName=GENE_NAME
  if (grepl("\\\\GName=([A-Za-z0-9_-]+)", desc)) {
    gene <- gsub(".*\\\\GName=([A-Za-z0-9_-]+).*", "\\1", desc)
    if (nchar(gene) >= 2 && nchar(gene) <= 15) {
      return(toupper(gene))
    }
  }
  
  # Pattern 3: Gene symbols in brackets
  if (grepl("\\[([A-Z][A-Z0-9_-]{1,14})\\]", desc)) {
    gene <- gsub(".*\\[([A-Z][A-Z0-9_-]{1,14})\\].*", "\\1", desc)
    if (!grepl("PROTEIN|FRAGMENT|DOMAIN", gene, ignore.case = TRUE)) {
      return(toupper(gene))
    }
  }
  
  return(NA)
}

# Test function
test_efficient_mapping <- function() {
  message("Testing efficient ID mapping...")
  
  test_ids <- c(
    "P02768",           # UniProt ALB
    "ENSP00000295897",  # Ensembl ALB
    "ALB",              # Gene symbol
    "Q13438",           # UniProt OS9
    "INVALID"           # Invalid
  )
  
  results <- efficient_map_to_gene_symbol(test_ids)
  
  for (i in seq_along(test_ids)) {
    message(paste(test_ids[i], "->", 
                  if (is.na(results[i])) "NOT_MAPPED" else results[i]))
  }
  
  return(results)
}

message("Efficient ID mapping functions loaded:")
message("- efficient_map_to_gene_symbol(accessions, descriptions, batch_size)")
message("- test_efficient_mapping()") 