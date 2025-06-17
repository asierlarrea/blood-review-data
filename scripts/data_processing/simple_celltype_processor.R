#!/usr/bin/env Rscript

# Simple Cell Type Data Processor
# Purpose: Process ProteomeXchange files with simple cell type columns
# Handles files like pxd025174.csv with Copy_Number_CD4, Copy_Number_CD8 columns

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Function to check if protein IDs map to the same gene
check_same_gene_mapping <- function(protein_ids) {
  # Split protein IDs by semicolon
  ids <- str_split(protein_ids, ";")[[1]]
  ids <- str_trim(ids)
  
  # Skip if no semicolons (single ID)
  if (length(ids) == 1) {
    return(TRUE)
  }
  
  # Map each ID to gene symbol
  gene_symbols <- convert_to_gene_symbol(ids)
  
  # Check if all mappings are successful and map to the same gene
  if (any(is.na(gene_symbols))) {
    return(FALSE)  # Some IDs couldn't be mapped
  }
  
  # Check if all map to the same gene
  return(length(unique(gene_symbols)) == 1)
}

# Function to process simple cell type files
process_simple_celltype_file <- function(file_path, 
                                        celltype_column_mapping = NULL,
                                        output_dir = NULL,
                                        force_mapping = FALSE) {
  
  message(sprintf("Processing simple cell type file: %s", basename(file_path)))
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  
  # Read the data
  data <- read_csv(file_path, show_col_types = FALSE)
  message(sprintf("  Loaded %d rows", nrow(data)))
  
  # Find protein ID column
  protein_id_cols <- c("Majority_Protein_IDs", "Majority protein IDs", "Protein_IDs", "Protein.IDs")
  protein_id_col <- protein_id_cols[protein_id_cols %in% colnames(data)][1]
  
  if (is.na(protein_id_col)) {
    warning("Could not find protein ID column, skipping unique protein filtering")
  } else {
    # Filter proteins based on gene mapping
    message("  Filtering proteins based on gene mapping...")
    data <- data %>%
      filter(check_same_gene_mapping(.data[[protein_id_col]]))
    message(sprintf("  Proteins after filtering: %d", nrow(data)))
  }
  
  # Find gene names column
  gene_cols <- c("Gene_names", "Gene names", "Gene.names", "QGene.names", "Genes")
  gene_col <- gene_cols[gene_cols %in% colnames(data)][1]
  
  if (is.na(gene_col)) {
    stop("Could not find gene names column")
  }
  
  # Auto-detect cell type columns if not provided
  if (is.null(celltype_column_mapping)) {
    # Look for common patterns
    all_cols <- colnames(data)
    
    # Pattern 1: Copy_Number_CELLTYPE
    copy_number_cols <- all_cols[str_detect(all_cols, "^Copy_Number_")]
    
    # Pattern 2: Intensity_CELLTYPE (but not the complex ones)
    simple_intensity_cols <- all_cols[str_detect(all_cols, "^Intensity_[A-Za-z0-9]+$")]
    
    # Pattern 3: CELLTYPE_expression or similar
    expression_cols <- all_cols[str_detect(all_cols, "_expression$|_abundance$|_count$")]
    
    # Pattern 4: LFQ.intensity columns (for CD8 T cell data)
    lfq_intensity_cols <- all_cols[str_detect(all_cols, "^LFQ\\.intensity\\.")]
    
    detected_cols <- c(copy_number_cols, simple_intensity_cols, expression_cols, lfq_intensity_cols)
    
    if (length(detected_cols) == 0) {
      stop("Could not auto-detect cell type columns. Please provide celltype_column_mapping.")
    }
    
    # For LFQ columns, we need to group them by cell type
    if (length(lfq_intensity_cols) > 0) {
      # Determine cell type based on filename and column patterns
      filename <- basename(file_path)
      if (str_detect(filename, "cd8")) {
        celltype_column_mapping <- list("CD8_T_cells" = lfq_intensity_cols)
      } else if (str_detect(filename, "macrophage")) {
        celltype_column_mapping <- list("Macrophages" = lfq_intensity_cols)
      } else {
        # Default: assume all are from same cell type based on filename
        celltype_name <- str_extract(filename, "[a-z_]+") %>% str_to_title()
        celltype_column_mapping <- list()
        celltype_column_mapping[[celltype_name]] <- lfq_intensity_cols
      }
    }
    
    # Create automatic mapping for other column types
    celltype_column_mapping <- list()
    for (col in detected_cols) {
      # Extract cell type name
      if (str_detect(col, "^Copy_Number_")) {
        celltype <- str_extract(col, "(?<=Copy_Number_).*")
      } else if (str_detect(col, "^Intensity_")) {
        celltype <- str_extract(col, "(?<=Intensity_).*")
      } else {
        celltype <- str_extract(col, ".*(?=_)")
      }
      
      # Standardize cell type name
      celltype_std <- case_when(
        str_detect(celltype, "CD4|T4") ~ "CD4_T_cells",
        str_detect(celltype, "CD8|T8") ~ "CD8_T_cells", 
        str_detect(celltype, "B|b_cell") ~ "B_cells",
        str_detect(celltype, "NK|nk") ~ "NK_cells",
        str_detect(celltype, "monocyte|Monocyte") ~ "Monocytes",
        str_detect(celltype, "macrophage|Macrophage") ~ "Macrophages",
        TRUE ~ celltype
      )
      
      celltype_column_mapping[[celltype_std]] <- col
    }
    
    message("  Auto-detected cell type columns:")
    for (ct in names(celltype_column_mapping)) {
      message(sprintf("    %s: %s", ct, celltype_column_mapping[[ct]]))
    }
  }
  
  # Load deduplication utility if not already loaded
  if (!exists("deduplicate_genes")) {
    source("scripts/utilities/gene_deduplication.R")
  }
  
  # Process each cell type
  results <- list()
  
  for (celltype in names(celltype_column_mapping)) {
    col_names <- celltype_column_mapping[[celltype]]
    
    # Handle both single column and multiple columns
    if (length(col_names) == 1) {
      col_name <- col_names
      if (!col_name %in% colnames(data)) {
        warning(sprintf("Column %s not found for cell type %s", col_name, celltype))
        next
      }
      
      message(sprintf("  Processing %s (column: %s)...", celltype, col_name))
      
      # Extract data for this cell type
      celltype_data <- data %>%
        select(all_of(gene_col), all_of(col_name)) %>%
        rename(gene = all_of(gene_col), intensity = all_of(col_name)) %>%
        filter(!is.na(gene) & gene != "" & !is.na(intensity) & intensity > 0)
      
    } else {
      # Multiple columns - calculate median across replicates
      missing_cols <- col_names[!col_names %in% colnames(data)]
      if (length(missing_cols) > 0) {
        warning(sprintf("Columns %s not found for cell type %s", 
                       paste(missing_cols, collapse = ", "), celltype))
        col_names <- col_names[col_names %in% colnames(data)]
      }
      
      if (length(col_names) == 0) {
        warning(sprintf("No valid columns found for cell type %s", celltype))
        next
      }
      
      message(sprintf("  Processing %s (%d columns)...", celltype, length(col_names)))
      
      # Extract data and calculate median intensity across columns
      celltype_data <- data %>%
        select(all_of(gene_col), all_of(col_names)) %>%
        rename(gene = all_of(gene_col))
      
      # Calculate median intensity (excluding zeros and NAs)
      celltype_data$intensity <- apply(celltype_data[, col_names, drop = FALSE], 1, function(x) {
        valid_values <- x[!is.na(x) & x > 0]
        if (length(valid_values) > 0) {
          median(valid_values)
        } else {
          NA
        }
      })
      
      celltype_data <- celltype_data %>%
        select(gene, intensity) %>%
        filter(!is.na(gene) & gene != "" & !is.na(intensity) & intensity > 0)
    }
    
    # For proteins with multiple accessions, use the first successfully mapped gene name
    celltype_data <- celltype_data %>%
      mutate(gene = str_split(gene, ";")) %>%
      unnest(gene) %>%
      mutate(gene = str_trim(gene)) %>%
      filter(gene != "") %>%
      group_by(intensity) %>%
      slice(1) %>%  # Keep only the first gene name for each intensity value
      ungroup()
    
    # Deduplicate genes using median intensity
    celltype_dedup <- deduplicate_genes(celltype_data, "gene", "intensity", 
                                      aggregation_method = "median")
    
    # Add metadata
    celltype_dedup$celltype <- celltype
    celltype_dedup$source <- paste0("ProteomeXchange_", str_remove(basename(file_path), "\\.csv"))
    celltype_dedup$filename <- basename(file_path)
    celltype_dedup$column_name <- ifelse(length(col_names) == 1, col_names, 
                                       paste(length(col_names), "columns"))
    
    results[[celltype]] <- celltype_dedup
    
    message(sprintf("    Final genes for %s: %d", celltype, nrow(celltype_dedup)))
  }
  
  # Combine results
  if (length(results) > 0) {
    combined_results <- bind_rows(results)
    
    # Save results if output directory is provided
    if (!is.null(output_dir)) {
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      
      output_file <- file.path(output_dir, 
                              paste0(str_remove(basename(file_path), "\\.csv"), 
                                    "_processed.csv"))
      write_csv(combined_results, output_file)
      message(sprintf("  Results saved to: %s", output_file))
    }
    
    return(combined_results)
  } else {
    message("  No valid cell type data found")
    return(NULL)
  }
}

# Function to process files with cell type in filename
process_filename_celltype_file <- function(file_path, 
                                         celltype = NULL,
                                         intensity_column = "abundance",
                                         id_column = "string_external_id",
                                         output_dir = NULL,
                                         force_mapping = FALSE) {
  
  message(sprintf("Processing filename-based cell type file: %s", basename(file_path)))
  
  # Auto-detect cell type from filename if not provided
  if (is.null(celltype)) {
    filename <- basename(file_path)
    
    # Extract from paxdb_CELLTYPE.csv pattern
    if (str_detect(filename, "^paxdb_")) {
      celltype_raw <- str_extract(filename, "(?<=paxdb_)[^.]+")
    } else {
      # Try to extract cell type from anywhere in filename
      celltype_raw <- str_extract(filename, "[a-z_]+")
    }
    
    # Standardize cell type name
    celltype <- case_when(
      str_detect(celltype_raw, "b_cell") ~ "B_cells",
      str_detect(celltype_raw, "cd4") ~ "CD4_T_cells",
      str_detect(celltype_raw, "cd8") ~ "CD8_T_cells",
      str_detect(celltype_raw, "nk") ~ "NK_cells",
      str_detect(celltype_raw, "monocyte") ~ "Monocytes",
      str_detect(celltype_raw, "plasma") ~ "Plasma",
      str_detect(celltype_raw, "serum") ~ "Serum",
      TRUE ~ celltype_raw
    )
    
    message(sprintf("  Auto-detected cell type: %s", celltype))
  }
  
  # Read data
  data <- read_csv(file_path, show_col_types = FALSE)
  message(sprintf("  Loaded %d rows", nrow(data)))
  
  # Load gene mapping utility if not already loaded
  if (!exists("convert_to_gene_symbol")) {
    source("scripts/data_processing/simple_id_mapping.R")
  }
  
  # Process ENSP IDs to gene symbols (for PAXDB files)
  if (str_detect(id_column, "string_external_id") && any(str_detect(data[[id_column]], "^9606\\."))) {
    message("  Converting ENSP IDs to gene symbols...")
    data$ensp <- str_replace(data[[id_column]], "^9606\\.", "")
    data$gene <- convert_to_gene_symbol(data$ensp, force_mapping = force_mapping)
    
    # Filter valid genes
    data_clean <- data %>%
      filter(!is.na(gene) & gene != "" & !is.na(.data[[intensity_column]])) %>%
      select(gene, intensity = all_of(intensity_column), ensp)
    
    # Load deduplication utility if not already loaded
    if (!exists("deduplicate_genes")) {
      source("scripts/utilities/gene_deduplication.R")
    }
    
    # Deduplicate genes using median abundance
    data_dedup <- deduplicate_genes(data_clean, "gene", "intensity", 
                                  additional_cols = c("ensp"), 
                                  aggregation_method = "median")
  } else {
    # Assume gene names are already present
    gene_col <- id_column
    data_clean <- data %>%
      filter(!is.na(.data[[gene_col]]) & .data[[gene_col]] != "" & 
             !is.na(.data[[intensity_column]])) %>%
      select(gene = all_of(gene_col), intensity = all_of(intensity_column))
    
    # Load deduplication utility if not already loaded
    if (!exists("deduplicate_genes")) {
      source("scripts/utilities/gene_deduplication.R")
    }
    
    # Deduplicate genes
    data_dedup <- deduplicate_genes(data_clean, "gene", "intensity", 
                                  aggregation_method = "median")
  }
  
  # Add metadata
  data_dedup$celltype <- celltype
  
  # Determine source based on file path
  if (str_detect(file_path, "paxdb")) {
    data_dedup$source <- "PAXDB"
  } else if (str_detect(file_path, "gpmdb")) {
    data_dedup$source <- "GPMDB"
  } else {
    data_dedup$source <- "Unknown"
  }
  
  data_dedup$filename <- basename(file_path)
  
  message(sprintf("  Final genes: %d", nrow(data_dedup)))
  
  # Save results if output directory is provided
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    output_file <- file.path(output_dir, 
                            paste0(str_remove(basename(file_path), "\\.csv"), 
                                  "_processed.csv"))
    write_csv(data_dedup, output_file)
    message(sprintf("  Results saved to: %s", output_file))
  }
  
  return(data_dedup)
}

# Export the functions if this script is sourced
if (!exists("SIMPLE_CELLTYPE_PROCESSOR_LOADED")) {
  SIMPLE_CELLTYPE_PROCESSOR_LOADED <- TRUE
  message("Simple cell type processor functions loaded")
} 