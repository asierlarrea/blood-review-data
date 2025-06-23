# ============================================================================
# UNIFIED DATA LOADER
# ============================================================================
# Description: Centralized data loading for all proteomics data sources
# Eliminates code duplication and ensures consistent data processing
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

# Source dependencies
source("scripts/config/analysis_config.R")
source("scripts/data_processing/simple_id_mapping.R")
source("scripts/utilities/gene_deduplication.R")

#' Load and process a single data source
#' 
#' @param source_name Name of the data source (must be in DATA_SOURCES config)
#' @param force_mapping Force re-mapping of gene IDs (default: FALSE)
#' @return Processed data frame with standardized columns
#' 
load_data_source <- function(source_name, force_mapping = FALSE) {
  
  # Get source configuration
  config <- get_data_source_config(source_name)
  
  # Construct file path
  file_path <- file.path(PROJECT_CONFIG$directories$data_raw, config$file_path)
  
  if (!file.exists(file_path)) {
    stop(paste("Data file not found:", file_path))
  }
  
  message(sprintf("Loading %s data from %s...", config$name, basename(file_path)))
  
  # Load raw data
  raw_data <- read_csv(file_path, 
                       show_col_types = FALSE, 
                       skip = config$skip_rows)
  
  message(sprintf("  - Input rows for %s: %d", config$name, nrow(raw_data)))
  
  # Preprocess ID column if needed
  if (!is.null(config$id_preprocessing)) {
    raw_data[[config$id_column]] <- config$id_preprocessing(raw_data[[config$id_column]])
  }
  
  # Handle gene extraction for sources like GPMDB
  if (!is.null(config$gene_extraction_pattern)) {
    # Use gene_column if specified, otherwise fall back to id_column
    extraction_column <- if (!is.null(config$gene_column)) config$gene_column else config$id_column
    raw_data$gene <- str_extract(raw_data[[extraction_column]], config$gene_extraction_pattern)
  } else if (config$requires_id_mapping) {
    # Perform ID mapping for sources that need it
    raw_data$gene <- convert_to_gene_symbol(raw_data[[config$id_column]], 
                                          force_mapping = force_mapping)
  } else {
    # Use existing gene column or rename
    if (config$id_column != "gene") {
      raw_data$gene <- raw_data[[config$id_column]]
    }
  }
  
  # Rename quantification column for consistency
  if (config$quant_column != "abundance") {
    raw_data$abundance <- raw_data[[config$quant_column]]
  }
  
  # Filter out invalid data
  processed_data <- raw_data %>%
    filter(!is.na(gene), 
           gene != "",
           !is.na(abundance),
           abundance > 0)
  
  # Add metadata columns
  processed_data$source <- config$name
  processed_data$technology <- config$technology
  processed_data$abundance_type <- config$abundance_type
  
  # Deduplicate genes
  if (nrow(processed_data) > 0) {
    processed_data <- deduplicate_genes(
      processed_data, 
      "gene", 
      "abundance",
      additional_cols = c("source", "technology", "abundance_type"),
      aggregation_method = PROJECT_CONFIG$analysis$aggregation_method
    )
  }
  
  # Validate processed data
  validate_data_source(processed_data, config)
  
  message(sprintf("  - Mapped and deduplicated genes for %s: %d", 
                  config$name, nrow(processed_data)))
  
  return(processed_data)
}

#' Load multiple data sources
#' 
#' @param source_names Vector of source names to load
#' @param force_mapping Force re-mapping of gene IDs
#' @return List of processed data frames
#' 
load_multiple_sources <- function(source_names = NULL, force_mapping = FALSE) {
  
  if (is.null(source_names)) {
    source_names <- get_all_data_sources()
  }
  
  # Validate source names
  invalid_sources <- setdiff(source_names, get_all_data_sources())
  if (length(invalid_sources) > 0) {
    stop(paste("Invalid source names:", paste(invalid_sources, collapse = ", ")))
  }
  
  message(sprintf("Loading %d data sources...", length(source_names)))
  
  # Load all sources
  data_list <- lapply(source_names, function(source) {
    tryCatch({
      load_data_source(source, force_mapping)
    }, error = function(e) {
      warning(sprintf("Failed to load %s: %s", source, e$message))
      return(NULL)
    })
  })
  
  # Remove failed loads
  names(data_list) <- source_names
  data_list <- data_list[!sapply(data_list, is.null)]
  
  message(sprintf("Successfully loaded %d out of %d sources", 
                  length(data_list), length(source_names)))
  
  return(data_list)
}

#' Combine multiple data sources into a single data frame
#' 
#' @param data_list List of data frames from load_multiple_sources()
#' @return Combined data frame
#' 
combine_data_sources <- function(data_list) {
  
  if (length(data_list) == 0) {
    stop("No data sources to combine")
  }
  
  # Ensure all data frames have the same essential columns
  required_cols <- c("gene", "abundance", "source", "technology", "abundance_type")
  
  for (i in seq_along(data_list)) {
    missing_cols <- setdiff(required_cols, colnames(data_list[[i]]))
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing columns in %s: %s", 
                   names(data_list)[i], paste(missing_cols, collapse = ", ")))
    }
  }
  
  # Combine data
  combined_data <- do.call(rbind, data_list)
  
  message(sprintf("Combined data: %d total entries, %d unique genes across %d sources",
                  nrow(combined_data), 
                  length(unique(combined_data$gene)),
                  length(unique(combined_data$source))))
  
  return(combined_data)
}

#' Load and process biomarker gene list
#' 
#' @param biomarker_file Path to biomarker file (optional, uses config default)
#' @return Vector of biomarker gene symbols
#' 
load_biomarker_genes <- function(biomarker_file = NULL) {
  
  if (is.null(biomarker_file)) {
    biomarker_file <- file.path(PROJECT_CONFIG$directories$data_metadata, 
                                BIOMARKER_CONFIG$biomarker_file)
  }
  
  if (!file.exists(biomarker_file)) {
    stop(paste("Biomarker file not found:", biomarker_file))
  }
  
  biomarkers <- read_tsv(biomarker_file, show_col_types = FALSE)
  biomarker_genes <- unique(toupper(biomarkers[[BIOMARKER_CONFIG$biomarker_column]]))
  biomarker_genes <- biomarker_genes[!is.na(biomarker_genes) & biomarker_genes != ""]
  
  message(sprintf("Loaded %d biomarker genes", length(biomarker_genes)))
  
  return(biomarker_genes)
}

#' Generate data loading summary statistics
#' 
#' @param data_list List of data frames or combined data frame
#' @return Summary data frame
#' 
generate_loading_summary <- function(data_list) {
  
  if (is.data.frame(data_list)) {
    # Single combined data frame
    summary_stats <- data_list %>%
      group_by(source, technology) %>%
      summarise(
        unique_genes = n_distinct(gene),
        total_entries = n(),
        median_abundance = median(abundance, na.rm = TRUE),
        abundance_type = first(abundance_type),
        .groups = "drop"
      )
  } else {
    # List of data frames
    summary_list <- lapply(names(data_list), function(source) {
      data <- data_list[[source]]
      data.frame(
        source = source,
        technology = data$technology[1],
        unique_genes = length(unique(data$gene)),
        total_entries = nrow(data),
        median_abundance = median(data$abundance, na.rm = TRUE),
        abundance_type = data$abundance_type[1],
        stringsAsFactors = FALSE
      )
    })
    summary_stats <- do.call(rbind, summary_list)
  }
  
  return(summary_stats)
}

#' Save processed data to cache
#' 
#' @param data_list List of processed data frames
#' @param cache_file Output file path
#' 
save_processed_data <- function(data_list, cache_file = NULL) {
  
  if (is.null(cache_file)) {
    cache_file <- file.path(PROJECT_CONFIG$directories$data_cache, 
                           "processed_plasma_data.rds")
  }
  
  # Create cache directory if it doesn't exist
  dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
  
  # Save data
  saveRDS(data_list, cache_file)
  message(sprintf("Saved processed data to %s", cache_file))
}

#' Load processed data from cache
#' 
#' @param cache_file Cache file path
#' @return List of processed data frames
#' 
load_processed_data <- function(cache_file = NULL) {
  
  if (is.null(cache_file)) {
    cache_file <- file.path(PROJECT_CONFIG$directories$data_cache, 
                           "processed_plasma_data.rds")
  }
  
  if (!file.exists(cache_file)) {
    stop(paste("Cache file not found:", cache_file))
  }
  
  data_list <- readRDS(cache_file)
  message(sprintf("Loaded processed data from %s", cache_file))
  
  return(data_list)
}

# Function to process GPMDB data with consistent gene mapping
process_gpmdb_data <- function(data, force_mapping = FALSE) {
  # Load required utilities if not already loaded
  if (!exists("convert_to_gene_symbol")) {
    source("scripts/data_processing/simple_id_mapping.R")
  }
  if (!exists("deduplicate_genes")) {
    source("scripts/utilities/gene_deduplication.R")
  }
  
  message("Processing GPMDB data...")
  message("  Input rows: ", nrow(data))
  
  # Step 1: Try to map transcript accessions to genes first
  message("  Mapping transcript accessions to genes...")
  data$gene_from_accession <- convert_to_gene_symbol(data$accession, force_mapping = force_mapping)
  
  # Step 2: Extract gene names from description as fallback
  message("  Extracting gene names from descriptions...")
  data$gene_from_desc <- stringr::str_extract(data$description, "[A-Z0-9]+(?=,| |$)")
  
  # Step 3: Combine gene mappings - use accession mapping if available, fallback to description
  data$gene <- ifelse(!is.na(data$gene_from_accession), 
                     data$gene_from_accession,
                     data$gene_from_desc)
  
  # Step 4: Clean up and filter invalid entries
  data_clean <- data %>%
    filter(!is.na(gene), gene != "", !is.na(total), total > 0) %>%
    select(gene, total, accession, description)
  
  message("  Valid genes after mapping: ", nrow(data_clean))
  
  # Step 5: Deduplicate genes using median
  data_dedup <- deduplicate_genes(data_clean, "gene", "total", 
                                additional_cols = c("accession", "description"),
                                aggregation_method = "median")
  
  message("  Final unique genes after deduplication: ", nrow(data_dedup))
  
  # Return processed data
  return(data_dedup)
}

#' Load and process QuantMS data from a specific directory
#' 
#' @param sample_type The sample type, e.g., "plasma" or "serum". This corresponds to the subdirectory in `quantms/`.
#' @param force_mapping Force re-mapping of protein IDs to gene symbols
#' @return A processed data frame for QuantMS data
#' 
load_quantms_data <- function(sample_type, force_mapping = FALSE) {
  
  quantms_dir <- file.path(PROJECT_CONFIG$directories$data_raw, "quantms", sample_type)
  
  if (!dir.exists(quantms_dir)) {
    warning(paste("QuantMS directory not found for sample type:", sample_type))
    return(NULL)
  }
  
  message(sprintf("Processing QuantMS %s data...", sample_type))
  
  # Get all CSV files
  quantms_files <- list.files(quantms_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(quantms_files) == 0) {
    warning(paste("No QuantMS CSV files found for sample type:", sample_type))
    return(NULL)
  }
  
  # Read and combine all files
  all_data <- map_dfr(quantms_files, ~{
    read_csv(.x, show_col_types = FALSE) %>%
      select(protein_accession = ProteinName, Ibaq = IbaqNorm, SampleID)
  })
  
  # Map protein accessions to gene symbols on the entire dataset
  all_data$gene <- convert_to_gene_symbol(all_data$protein_accession, force_mapping = force_mapping)
  
  # Filter out proteins that couldn't be mapped
  mapped_data <- all_data %>%
    filter(!is.na(gene) & gene != "" & !is.na(Ibaq) & Ibaq > 0) %>%
    select(gene, Ibaq, SampleID)
  
  # Filter genes that are present in 3 or more samples
  gene_sample_counts <- mapped_data %>%
    group_by(gene) %>%
    summarise(sample_count = n_distinct(SampleID), .groups = "drop")
  
  genes_with_sufficient_samples <- gene_sample_counts %>%
    filter(sample_count >= 3) %>%
    pull(gene)
  
  message(sprintf("  - Genes present in >=3 samples: %d out of %d total genes", 
                  length(genes_with_sufficient_samples), nrow(gene_sample_counts)))
  
  # Filter mapped data to include only genes present in >=3 samples
  filtered_data <- mapped_data %>%
    filter(gene %in% genes_with_sufficient_samples)
  
  # Now, aggregate by gene
  processed_data <- filtered_data %>%
    group_by(gene) %>%
    summarise(
      abundance = median(Ibaq, na.rm = TRUE),
      sample_count = n_distinct(SampleID),
      protein_count = n(),
      .groups = "drop"
    )
    
  message(sprintf("  - Mapped to %d unique genes.", n_distinct(processed_data$gene)))
  
  processed_data <- processed_data %>%
    filter(!is.na(abundance) & abundance > 0)
  
  # Add metadata
  processed_data$source <- "quantms"
  processed_data$technology <- "MS"
  processed_data$abundance_type <- "iBAQ Normalized"
  
  # Keep the standard columns for compatibility (including protein_count to match other sources)
  final_data <- processed_data %>%
    select(gene, abundance, protein_count, source, technology, abundance_type)
  
  return(final_data)
} 