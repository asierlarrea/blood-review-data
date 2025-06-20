# ============================================================================
# CENTRALIZED ANALYSIS CONFIGURATION
# ============================================================================
# Description: Central configuration for all blood proteomics analyses
# Contains: Data sources, file paths, plot themes, analysis parameters
# ============================================================================

# Project Configuration
PROJECT_CONFIG <- list(
  # Project metadata
  project_name = "Blood Review Data Analysis",
  version = "2.0.0",
  
  # Directory structure
  directories = list(
    data_raw = "data/raw",
    data_cache = "data/cache", 
    data_metadata = "data/metadata",
    outputs = "outputs",
    plots = "outputs/plots",
    tables = "outputs/tables",
    reports = "outputs/reports",
    logs = "outputs/logs"
  ),
  
  # Analysis parameters
  analysis = list(
    force_mapping = FALSE,
    aggregation_method = "median", # for gene deduplication
    log_transform_offset = 1,
    z_score_normalize = TRUE,
    remove_outliers = FALSE,
    outlier_threshold = 3 # standard deviations
  )
)

# Data Source Configuration
DATA_SOURCES <- list(
  peptideatlas = list(
    name = "PeptideAtlas",
    technology = "MS",
    file_path = "peptideatlas/peptideatlas.csv",
    id_column = "biosequence_accession",
    gene_column = "biosequence_accession",
    quant_column = "norm_PSMs_per_100K",
    skip_rows = 0,
    requires_id_mapping = TRUE,
    abundance_type = "Normalized PSMs per 100K"
  ),
  
  hpa_ms = list(
    name = "HPA MS",
    technology = "MS", 
    file_path = "hpa/hpa_ms.csv",
    id_column = "Gene",
    gene_column = "Gene",
    quant_column = "Concentration",
    skip_rows = 1,
    requires_id_mapping = FALSE,
    abundance_type = "mg/L"
  ),
  
  hpa_pea = list(
    name = "HPA PEA",
    technology = "PEA",
    file_path = "hpa/hpa_pea.csv", 
    id_column = "Gene",
    gene_column = "Gene",
    quant_column = "Variation between individuals",
    skip_rows = 0,
    requires_id_mapping = FALSE,
    abundance_type = "Coefficient of Variation"
  ),
  
  hpa_immunoassay = list(
    name = "HPA Immunoassay",
    technology = "Immunoassay",
    file_path = "hpa/hpa_immunoassay_plasma.csv",
    id_column = "Gene", 
    gene_column = "Gene",
    quant_column = "Concentration",
    skip_rows = 0,
    requires_id_mapping = FALSE,
    abundance_type = "mg/L"
  ),
  
  gpmdb = list(
    name = "GPMDB",
    technology = "MS",
    file_path = "gpmdb/gpmdb_plasma.csv",
    id_column = "accession",
    gene_column = "description", 
    quant_column = "total",
    skip_rows = 0,
    requires_id_mapping = FALSE,
    abundance_type = "Spectral Count",
    gene_extraction_pattern = "[A-Z0-9]+(?=,| |$)"
  ),
  
  paxdb = list(
    name = "PAXDB",
    technology = "MS",
    file_path = "paxdb/paxdb_plasma.csv",
    id_column = "string_external_id",
    gene_column = "string_external_id",
    quant_column = "abundance", 
    skip_rows = 0,
    requires_id_mapping = TRUE,
    abundance_type = "ppm",
    id_preprocessing = function(x) stringr::str_replace(x, "^9606\\.", "")
  )
)

# Plot Configuration
PLOT_CONFIG <- list(
  # Color palettes
  colors = list(
    technology = c(
      "MS" = "#2E86AB",
      "PEA" = "#A23B72", 
      "Immunoassay" = "#F18F01"
    ),
    biomarker = c(
      "Biomarker" = "#E69F00",
      "Other" = "#69b3a2"
    ),
    databases = c(
      "PeptideAtlas" = "#2E86AB",
      "HPA MS" = "#4A90E2",
      "HPA PEA" = "#A23B72",
      "HPA Immunoassay" = "#F18F01", 
      "GPMDB" = "#50C878",
      "PAXDB" = "#9B59B6",
      "quantms" = "#D55E00",
      "quantms_max" = "#CC79A7"
    )
  ),
  
  # Default theme settings
  theme = list(
    base_size = 12,
    title_size = 16,
    subtitle_size = 12,
    axis_title_size = 12,
    axis_text_size = 10,
    legend_text_size = 10
  ),
  
  # Plot dimensions
  dimensions = list(
    default_width = 10,
    default_height = 6,
    large_width = 12,
    large_height = 8,
    combined_width = 16,
    combined_height = 20,
    dpi = 600
  ),
  
  # Plot output formats (can be overridden by command line)
  formats = c("tiff", "png")  # Default to TIFF and PNG for high-resolution output
)

# Analysis-specific configurations
BIOMARKER_CONFIG <- list(
  biomarker_file = "biomarkers_list.tsv",
  biomarker_column = "To",
  highlight_color = "#E69F00",
  violin_plot_quantiles = c(0.25, 0.5, 0.75)
)

# Utility functions for configuration
get_data_source_config <- function(source_name) {
  if (source_name %in% names(DATA_SOURCES)) {
    return(DATA_SOURCES[[source_name]])
  } else {
    stop(paste("Unknown data source:", source_name))
  }
}

get_plot_colors <- function(type = "technology") {
  return(PLOT_CONFIG$colors[[type]])
}

get_all_data_sources <- function() {
  return(names(DATA_SOURCES))
}

get_ms_sources <- function() {
  ms_sources <- sapply(DATA_SOURCES, function(x) x$technology == "MS")
  return(names(DATA_SOURCES)[ms_sources])
}

# Configuration functions
set_plot_formats <- function(formats_string) {
  if (!is.null(formats_string) && nchar(formats_string) > 0) {
    # Parse comma-separated formats
    formats <- trimws(unlist(strsplit(formats_string, ",")))
    # Validate formats
    valid_formats <- c("png", "svg", "pdf", "jpeg", "tiff")
    invalid_formats <- setdiff(formats, valid_formats)
    if (length(invalid_formats) > 0) {
      warning(paste("Invalid plot formats ignored:", paste(invalid_formats, collapse = ", ")))
      formats <- intersect(formats, valid_formats)
    }
    if (length(formats) > 0) {
      PLOT_CONFIG$formats <<- formats
      message(sprintf("Plot formats set to: %s", paste(formats, collapse = ", ")))
    } else {
      warning("No valid plot formats provided, using default")
    }
  }
}

get_plot_formats <- function() {
  return(PLOT_CONFIG$formats)
}

# Validation functions
validate_data_source <- function(data, source_config) {
  # Check for processed column names that the data loader creates
  required_cols <- c("gene", "abundance", "source", "technology", "abundance_type")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in processed", source_config$name, "data:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # Also validate data content
  if (nrow(data) == 0) {
    stop(paste("No valid data rows in", source_config$name))
  }
  
  return(TRUE)
} 