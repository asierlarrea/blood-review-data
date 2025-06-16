#!/usr/bin/env Rscript

# ProteomeXchange Data Processor
# Purpose: Efficiently process large ProteomeXchange files with intensity columns
# Handles files like pxd004352.csv that have multiple intensity columns per cell type

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
})

# Function to process ProteomeXchange files with intensity columns
process_proteomexchange_intensity_file <- function(file_path, 
                                                 output_dir = NULL,
                                                 force_mapping = FALSE,
                                                 chunk_size = 10000) {
  
  message(sprintf("Processing ProteomeXchange file: %s", basename(file_path)))
  
  # Check if the ID mapping utility is available
  if (!exists("convert_to_gene_symbol")) {
    source("scripts/data_processing/simple_id_mapping.R")
  }
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  
  # Read header to understand column structure
  message("  Analyzing column structure...")
  header <- read_csv(file_path, n_max = 0, show_col_types = FALSE)
  all_cols <- colnames(header)
  
  # Find intensity columns
  intensity_cols <- all_cols[str_detect(all_cols, "^Intensity_")]
  message(sprintf("  Found %d intensity columns", length(intensity_cols)))
  
  # Extract cell types from column names
  extract_celltype <- function(col_name) {
    if (str_detect(col_name, "^Intensity_")) {
      celltype <- str_extract(col_name, "(?<=Intensity_)[^_]+")
      return(celltype)
    }
    return(NA)
  }
  
  # Standardize cell type names
  standardize_celltype <- function(celltype) {
    celltype_mapping <- c(
      "B.memory" = "B_cells",
      "B.naive" = "B_cells", 
      "T4.naive" = "CD4_T_cells",
      "T4.CM" = "CD4_T_cells",
      "T4.EM" = "CD4_T_cells",
      "mTregs" = "CD4_T_cells",
      "nTregs" = "CD4_T_cells",
      "T8.naive" = "CD8_T_cells",
      "T8.CM" = "CD8_T_cells", 
      "T8.EM" = "CD8_T_cells"
    )
    
    if (celltype %in% names(celltype_mapping)) {
      return(celltype_mapping[[celltype]])
    }
    return(celltype)
  }
  
  # Create mapping of columns to cell types
  celltype_mapping <- tibble(
    column = intensity_cols,
    celltype_raw = map_chr(intensity_cols, extract_celltype),
    celltype = map_chr(celltype_raw, standardize_celltype)
  ) %>%
    filter(!is.na(celltype_raw))
  
  # Group columns by standardized cell type
  celltype_groups <- celltype_mapping %>%
    group_by(celltype) %>%
    summarise(columns = list(column), .groups = "drop")
  
  message("  Cell types found:")
  for (i in 1:nrow(celltype_groups)) {
    ct <- celltype_groups$celltype[i]
    col_count <- length(celltype_groups$columns[[i]])
    message(sprintf("    %s: %d columns", ct, col_count))
  }
  
  # Function to check if protein ID is unique (no semicolons)
  is_unique_protein <- function(protein_ids) {
    return(!str_detect(protein_ids, ";"))
  }
  
  # Read the full dataset
  message("  Reading full dataset...")
  
  # For very large files, we might want to process in chunks
  # For now, let's try to read the whole file
  tryCatch({
    data <- read_csv(file_path, show_col_types = FALSE)
    message(sprintf("  Loaded %d rows", nrow(data)))
  }, error = function(e) {
    stop(sprintf("Failed to read file: %s", e$message))
  })
  
  # Filter out proteins with multiple IDs
  message("  Filtering unique proteins...")
  protein_id_col <- "Majority protein IDs"
  if (!protein_id_col %in% colnames(data)) {
    # Try alternative column names
    possible_cols <- c("Majority_Protein_IDs", "Protein_IDs", "Protein.IDs")
    found_col <- possible_cols[possible_cols %in% colnames(data)]
    if (length(found_col) > 0) {
      protein_id_col <- found_col[1]
    } else {
      stop("Could not find protein ID column")
    }
  }
  
  data_clean <- data %>%
    filter(is_unique_protein(.data[[protein_id_col]]))
  
  message(sprintf("  Proteins after filtering non-unique IDs: %d", nrow(data_clean)))
  
  # Process each cell type
  results <- list()
  
  for (i in 1:nrow(celltype_groups)) {
    celltype <- celltype_groups$celltype[i]
    columns <- celltype_groups$columns[[i]]
    
    message(sprintf("  Processing %s (%d columns)...", celltype, length(columns)))
    
    # Calculate median intensity across columns for this cell type
    if (length(columns) == 1) {
      data_clean[[paste0(celltype, "_intensity")]] <- data_clean[[columns[1]]]
    } else {
      # Calculate row-wise median across multiple columns
      intensity_matrix <- data_clean[columns]
      data_clean[[paste0(celltype, "_intensity")]] <- apply(intensity_matrix, 1, function(x) {
        valid_values <- x[!is.na(x) & x > 0]
        if (length(valid_values) > 0) {
          return(median(valid_values))
        } else {
          return(NA)
        }
      })
    }
    
    # Extract data for this cell type
    gene_col <- "Gene names"
    if (!gene_col %in% colnames(data_clean)) {
      # Try alternative column names
      possible_cols <- c("Gene_names", "Gene.names", "Genes")
      found_col <- possible_cols[possible_cols %in% colnames(data_clean)]
      if (length(found_col) > 0) {
        gene_col <- found_col[1]
      } else {
        warning("Could not find gene names column, skipping this cell type")
        next
      }
    }
    
    celltype_data <- data_clean %>%
      select(all_of(gene_col), !!paste0(celltype, "_intensity")) %>%
      rename(gene = all_of(gene_col), intensity = !!paste0(celltype, "_intensity")) %>%
      filter(!is.na(gene) & gene != "" & !is.na(intensity) & intensity > 0)
    
    # Handle multiple gene names (separated by semicolons)
    celltype_data <- celltype_data %>%
      mutate(gene = str_split(gene, ";")) %>%
      unnest(gene) %>%
      mutate(gene = str_trim(gene)) %>%
      filter(gene != "")
    
    # Load deduplication utility if not already loaded
    if (!exists("deduplicate_genes")) {
      source("scripts/utilities/gene_deduplication.R")
    }
    
    # Deduplicate genes using median intensity
    celltype_dedup <- deduplicate_genes(celltype_data, "gene", "intensity", 
                                      aggregation_method = "median")
    
    # Add metadata
    celltype_dedup$celltype <- celltype
    celltype_dedup$source <- paste0("ProteomeXchange_", str_remove(basename(file_path), "\\.csv"))
    celltype_dedup$filename <- basename(file_path)
    celltype_dedup$column_count <- length(columns)
    
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

# Function to get a quick preview of ProteomeXchange file structure
preview_proteomexchange_file <- function(file_path, n_rows = 5) {
  
  message(sprintf("Previewing file: %s", basename(file_path)))
  
  # Read header
  header <- read_csv(file_path, n_max = 0, show_col_types = FALSE)
  all_cols <- colnames(header)
  
  message(sprintf("Total columns: %d", length(all_cols)))
  
  # Find intensity columns
  intensity_cols <- all_cols[str_detect(all_cols, "^Intensity_")]
  message(sprintf("Intensity columns: %d", length(intensity_cols)))
  
  if (length(intensity_cols) > 0) {
    message("Sample intensity columns:")
    for (i in 1:min(10, length(intensity_cols))) {
      message(sprintf("  %s", intensity_cols[i]))
    }
    if (length(intensity_cols) > 10) {
      message(sprintf("  ... and %d more", length(intensity_cols) - 10))
    }
  }
  
  # Show other important columns
  important_cols <- all_cols[str_detect(tolower(all_cols), "protein|gene|id")]
  if (length(important_cols) > 0) {
    message("Important columns (protein/gene/id):")
    for (col in important_cols) {
      message(sprintf("  %s", col))
    }
  }
  
  # Read a few sample rows
  if (n_rows > 0) {
    message(sprintf("\nFirst %d rows:", n_rows))
    sample_data <- read_csv(file_path, n_max = n_rows, show_col_types = FALSE)
    print(sample_data[1:min(5, ncol(sample_data))])
  }
}

# Export the functions if this script is sourced
if (!exists("PROTEOMEXCHANGE_PROCESSOR_LOADED")) {
  PROTEOMEXCHANGE_PROCESSOR_LOADED <- TRUE
  message("ProteomeXchange processor functions loaded")
} 