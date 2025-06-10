#!/usr/bin/env Rscript

# Enhanced Fast Local ID Mapping with biomaRt ENSP mappings
# Author: Data Analysis Pipeline
# Date: 2024

# Load the basic fast mapping functions
source("fast_id_mapping.R")

# Function to load all biomaRt ENSP mappings
load_biomart_ensp_mappings <- function() {
  if (file.exists("biomart_ensp_mappings.csv")) {
    message("Loading biomaRt ENSP mappings...")
    biomart_data <- read.csv("biomart_ensp_mappings.csv", stringsAsFactors = FALSE)
    
    # Create named list for fast lookup
    ensp_mappings <- setNames(toupper(biomart_data$gene_symbol), biomart_data$ensp_id)
    
    message(paste("Loaded", length(ensp_mappings), "ENSP mappings from biomaRt"))
    return(ensp_mappings)
  } else {
    message("biomart_ensp_mappings.csv not found - using local mappings only")
    return(list())
  }
}

# Enhanced mapping function that uses all biomaRt results
enhanced_fast_map_to_gene_symbol <- function(accessions, descriptions = NULL) {
  if (length(accessions) == 0) {
    return(character(0))
  }
  
  message(paste("Enhanced mapping", length(accessions), "IDs to gene symbols..."))
  
  # Load biomaRt ENSP mappings once
  biomart_ensp_mappings <- load_biomart_ensp_mappings()
  
  gene_symbols <- character(length(accessions))
  
  # Process each accession
  for(i in seq_along(accessions)) {
    acc <- as.character(accessions[i])
    desc <- if(!is.null(descriptions)) as.character(descriptions[i]) else ""
    
    # Try enhanced mapping
    gene_symbols[i] <- enhanced_map_single_id(acc, desc, biomart_ensp_mappings)
  }
  
  # Count results
  mapped_count <- sum(!is.na(gene_symbols) & gene_symbols != "")
  message(paste("Enhanced mapped:", mapped_count, "out of", length(accessions), "IDs"))
  
  return(gene_symbols)
}

# Enhanced single ID mapping
enhanced_map_single_id <- function(accession, description = "", biomart_ensp_mappings) {
  acc <- as.character(accession)
  desc <- as.character(description)
  
  # 1. If it's already a gene symbol
  if (grepl("^[A-Z][A-Z0-9_-]{1,14}$", acc) && !grepl("^ENS|^[OPQ][0-9]|^[A-Z][0-9]", acc)) {
    return(toupper(acc))
  }
  
  # 2. Check biomaRt ENSP mappings first (comprehensive)
  if (grepl("^ENSP[0-9]+", acc) && acc %in% names(biomart_ensp_mappings)) {
    return(biomart_ensp_mappings[[acc]])
  }
  
  # 3. Fall back to local common protein mappings
  common_gene <- get_common_protein_mapping(acc)
  if (!is.na(common_gene)) {
    return(common_gene)
  }
  
  # 4. Parse from description
  if (!is.na(desc) && desc != "") {
    desc_gene <- extract_gene_from_description_fast(desc)
    if (!is.na(desc_gene)) {
      return(desc_gene)
    }
  }
  
  # 5. For ENSP IDs without biomaRt mapping, try description parsing
  if (grepl("^ENSP[0-9]+", acc) && !is.na(desc) && desc != "") {
    ensp_gene <- extract_gene_from_ensembl_desc(desc)
    if (!is.na(ensp_gene)) {
      return(ensp_gene)
    }
  }
  
  return(NA)
}

# Test enhanced mapping
test_enhanced_mapping <- function() {
  message("Testing enhanced ID mapping...")
  
  test_ids <- c(
    "P02768",           # UniProt ALB
    "ENSP00000295897",  # Ensembl ALB (should be in biomaRt)
    "ENSP00000245907",  # ENSP from biomaRt (C3)
    "ALB",              # Gene symbol
    "Q13438",           # UniProt OS9
    "INVALID"           # Invalid
  )
  
  test_descriptions <- c(
    "Serum albumin GN=ALB PE=1 SV=2",
    "", 
    "",
    "",
    "Protein OS-9 GN=OS9 PE=1 SV=1",
    ""
  )
  
  results <- enhanced_fast_map_to_gene_symbol(test_ids, test_descriptions)
  
  for (i in seq_along(test_ids)) {
    message(paste(test_ids[i], "->", 
                  if (is.na(results[i])) "NOT_MAPPED" else results[i]))
  }
  
  return(results)
}

message("Enhanced fast ID mapping functions loaded:")
message("- enhanced_fast_map_to_gene_symbol(accessions, descriptions)")
message("- test_enhanced_mapping()")
message("- Uses all 7,013 biomaRt ENSP mappings for maximum coverage") 