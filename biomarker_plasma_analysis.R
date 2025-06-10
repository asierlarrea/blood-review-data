#!/usr/bin/env Rscript

# Biomarker Gene Expression Analysis in Plasma Sources
# Analyzes biomarker genes across different plasma databases and shows their position
# in the overall expression distribution for each source

# Set CRAN mirror
if(is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

# Function to install packages if not available
install_if_missing <- function(packages) {
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      install.packages(pkg, dependencies = TRUE)
      if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
        stop(paste("Failed to install package:", pkg))
      }
    }
  }
}

# Install required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "ggridges", "viridis", 
                      "scales", "gridExtra", "stringr", "ggbeeswarm", "readr")
install_if_missing(required_packages)

# Create output directory
output_dir <- "plots/04_Biomarker_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load enhanced ID mapping system
message("Loading enhanced ID mapping system...")
if(file.exists("enhanced_fast_mapping.R")) {
  source("enhanced_fast_mapping.R")
} else {
  stop("Enhanced mapping system not found. Please run the enhanced_fast_mapping.R script first.")
}

# Function to load biomarker genes
load_biomarker_genes <- function() {
  message("Loading biomarker gene list...")
  
  if(!file.exists("biomarkers_list.tsv")) {
    stop("biomarkers_list.tsv not found!")
  }
  
  biomarkers <- read.delim("biomarkers_list.tsv", stringsAsFactors = FALSE)
  
  # Extract unique gene symbols
  biomarker_genes <- unique(toupper(biomarkers$To))
  biomarker_genes <- biomarker_genes[!is.na(biomarker_genes) & biomarker_genes != ""]
  
  message(paste("Loaded", length(biomarker_genes), "unique biomarker genes"))
  return(biomarker_genes)
}

# Function to load plasma data from all sources
load_plasma_data <- function() {
  message("Loading plasma data from all sources...")
  
  plasma_data <- list()
  
  # 1. PeptideAtlas
  message("  - Loading PeptideAtlas plasma data...")
  if(file.exists("PeptideAtlas.csv")) {
    peptideatlas <- read.csv("PeptideAtlas.csv", stringsAsFactors = FALSE)
    
    gene_symbols <- enhanced_fast_map_to_gene_symbol(peptideatlas$biosequence_accession, 
                                                   peptideatlas$biosequence_desc)
    
    # Filter valid mappings
    valid_idx <- !is.na(gene_symbols) & gene_symbols != ""
    
    if(sum(valid_idx) > 0) {
      plasma_data$PeptideAtlas <- data.frame(
        gene_symbol = gene_symbols[valid_idx],
        abundance = as.numeric(peptideatlas$norm_PSMs_per_100K[valid_idx]),
        log_abundance = log10(as.numeric(peptideatlas$norm_PSMs_per_100K[valid_idx]) + 1),
        source = "PeptideAtlas",
        abundance_type = "Normalized PSMs per 100K",
        stringsAsFactors = FALSE
      )
      
      # Remove rows with invalid abundance
      plasma_data$PeptideAtlas <- plasma_data$PeptideAtlas[
        !is.na(plasma_data$PeptideAtlas$abundance) & plasma_data$PeptideAtlas$abundance > 0, ]
    }
  }
  
  # 2. PaxDB Plasma
  message("  - Loading PaxDB plasma data...")
  if(file.exists("PaxDb_Plasma.csv")) {
    lines <- readLines("PaxDb_Plasma.csv")
    data_start <- which(grepl("^string_external_id,abundance", lines))
    if(length(data_start) > 0) {
      paxdb_data <- read.csv("PaxDb_Plasma.csv", skip = data_start - 1, 
                           stringsAsFactors = FALSE, header = TRUE)
      
      # Clean protein IDs and map to gene symbols
      protein_ids <- gsub("^9606\\.", "", paxdb_data$string_external_id)
      gene_symbols <- enhanced_fast_map_to_gene_symbol(protein_ids)
      
      # Filter valid mappings
      valid_idx <- !is.na(gene_symbols) & gene_symbols != "" & 
                  !is.na(paxdb_data$abundance) & paxdb_data$abundance > 0
      
      if(sum(valid_idx) > 0) {
        plasma_data$PaxDB <- data.frame(
          gene_symbol = gene_symbols[valid_idx],
          abundance = paxdb_data$abundance[valid_idx],
          log_abundance = log10(paxdb_data$abundance[valid_idx] + 1),
          source = "PaxDB",
          abundance_type = "ppm",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # 3. HPA Plasma data
  message("  - Loading HPA plasma data...")
  hpa_files <- c("HPA_MS.csv", "HPA_Immunoassay.csv")
  hpa_data_list <- list()
  
  for(file in hpa_files) {
    if(file.exists(file)) {
      technique <- gsub("HPA_|\\.csv", "", file)
      
      hpa_data <- read.csv(file, stringsAsFactors = FALSE)
      
      # Find concentration column
      concentration_col <- NULL
      if("Mean.concentration" %in% colnames(hpa_data)) {
        concentration_col <- "Mean.concentration"
      } else if("Concentration" %in% colnames(hpa_data)) {
        concentration_col <- "Concentration"
      }
      
      if(!is.null(concentration_col)) {
        # Filter rows with valid gene names and concentrations
        valid_rows <- !is.na(hpa_data$Gene) & hpa_data$Gene != "" & 
                     !is.na(hpa_data[[concentration_col]]) & hpa_data[[concentration_col]] != ""
        
        if(sum(valid_rows) > 0) {
          valid_data <- hpa_data[valid_rows, ]
          
          # Convert concentration to numeric
          conc_strings <- as.character(valid_data[[concentration_col]])
          conc_strings <- iconv(conc_strings, to = "UTF-8", sub = "")
          concentration_values <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", conc_strings)))
          
          # Remove rows where concentration conversion failed
          numeric_valid <- !is.na(concentration_values) & concentration_values > 0
          
          if(sum(numeric_valid) > 0) {
            gene_symbols <- toupper(valid_data$Gene[numeric_valid])
            
            hpa_data_list[[technique]] <- data.frame(
              gene_symbol = gene_symbols,
              abundance = concentration_values[numeric_valid],
              log_abundance = log10(concentration_values[numeric_valid] + 0.001),
              source = paste0("HPA_", technique),
              abundance_type = "mg/L",
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  
  # Combine HPA data
  if(length(hpa_data_list) > 0) {
    plasma_data$HPA <- do.call(rbind, hpa_data_list)
    rownames(plasma_data$HPA) <- NULL
  }
  
  # 4. GPMDB Plasma
  message("  - Loading GPMDB plasma data...")
  if(file.exists("GPMDB_Plasma.csv")) {
    lines <- readLines("GPMDB_Plasma.csv")
    header_line <- which(grepl("^#,accession,total", lines))
    if(length(header_line) > 0) {
      gpmdb_data <- read.csv("GPMDB_Plasma.csv", skip = header_line, 
                           stringsAsFactors = FALSE, header = FALSE)
      
      colnames(gpmdb_data) <- c("rank", "accession", "total", "log_e", "EC", "description")
      
      # Filter valid rows
      valid_rows <- !is.na(gpmdb_data$accession) & gpmdb_data$accession != "" & 
                   !is.na(gpmdb_data$total) & gpmdb_data$total > 0
      
      if(sum(valid_rows) > 0) {
        gene_symbols <- enhanced_fast_map_to_gene_symbol(gpmdb_data$accession[valid_rows], 
                                                       gpmdb_data$description[valid_rows])
        
        # Filter valid mappings
        mapped_idx <- !is.na(gene_symbols) & gene_symbols != ""
        
        if(sum(mapped_idx) > 0) {
          plasma_data$GPMDB <- data.frame(
            gene_symbol = gene_symbols[mapped_idx],
            abundance = gpmdb_data$total[valid_rows][mapped_idx],
            log_abundance = log10(gpmdb_data$total[valid_rows][mapped_idx] + 1),
            source = "GPMDB",
            abundance_type = "Spectral Counts",
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  # Combine all plasma data
  all_plasma <- do.call(rbind, plasma_data)
  rownames(all_plasma) <- NULL
  
  message(paste("Total plasma measurements:", nrow(all_plasma)))
  message(paste("Unique genes in plasma:", length(unique(all_plasma$gene_symbol))))
  message(paste("Sources:", paste(unique(all_plasma$source), collapse = ", ")))
  
  return(all_plasma)
}

# Function to analyze biomarkers in plasma data
analyze_biomarkers_in_plasma <- function(plasma_data, biomarker_genes) {
  message("Analyzing biomarkers in plasma data...")
  
  # Mark biomarker genes
  plasma_data$is_biomarker <- plasma_data$gene_symbol %in% biomarker_genes
  
  # Get biomarker subset
  biomarker_data <- plasma_data[plasma_data$is_biomarker, ]
  
  # Summary statistics
  biomarker_summary <- biomarker_data %>%
    group_by(source) %>%
    summarise(
      n_biomarkers = length(unique(gene_symbol)),
      total_biomarkers = length(biomarker_genes),
      coverage_pct = round(100 * n_biomarkers / length(biomarker_genes), 1),
      median_abundance = median(abundance, na.rm = TRUE),
      mean_log_abundance = mean(log_abundance, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Overall summary
  overall_summary <- plasma_data %>%
    group_by(source) %>%
    summarise(
      total_genes = length(unique(gene_symbol)),
      median_abundance = median(abundance, na.rm = TRUE),
      mean_log_abundance = mean(log_abundance, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Combine summaries
  summary_table <- merge(overall_summary, biomarker_summary, 
                        by = "source", suffixes = c("_all", "_biomarkers"))
  
  print(summary_table)
  
  return(list(
    plasma_data = plasma_data,
    biomarker_data = biomarker_data,
    summary = summary_table
  ))
}

# Function to create biomarker expression plots
create_biomarker_plots <- function(analysis_results) {
  plasma_data <- analysis_results$plasma_data
  biomarker_data <- analysis_results$biomarker_data
  
  # Plot 1: Distribution with biomarker overlay
  message("Creating biomarker distribution overlay plot...")
  
  tiff(file.path(output_dir, "biomarker_distribution_overlay.tiff"), 
       width = 16, height = 10, units = "in", res = 600, compression = "lzw")
  
  p1 <- ggplot(plasma_data, aes(x = log_abundance, y = source, fill = source)) +
    geom_density_ridges(alpha = 0.7, scale = 2) +
    geom_point(data = biomarker_data, 
               aes(x = log_abundance, y = source), 
               color = "red", size = 2, alpha = 0.8,
               position = position_jitter(height = 0.1)) +
    scale_fill_viridis_d(name = "Source") +
    theme_minimal() +
    labs(title = "Biomarker Gene Expression in Plasma Sources",
         subtitle = "Red points show biomarker genes overlaid on expression distributions",
         x = "Log10(Abundance)", y = "Data Source") +
    theme(legend.position = "none",
          plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12))
  
  print(p1)
  dev.off()
  
  # Plot 2: Biomarker vs background comparison
  message("Creating biomarker vs background comparison...")
  
  tiff(file.path(output_dir, "biomarker_vs_background.tiff"), 
       width = 14, height = 10, units = "in", res = 600, compression = "lzw")
  
  comparison_data <- plasma_data %>%
    mutate(gene_type = ifelse(is_biomarker, "Biomarker", "Background"))
  
  p2 <- ggplot(comparison_data, aes(x = source, y = log_abundance, fill = gene_type)) +
    geom_boxplot(alpha = 0.8, outlier.alpha = 0.3) +
    scale_fill_manual(values = c("Background" = "lightblue", "Biomarker" = "red"),
                     name = "Gene Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Biomarker vs Background Gene Expression",
         subtitle = "Comparing biomarker genes to overall expression background",
         x = "Data Source", y = "Log10(Abundance)")
  
  print(p2)
  dev.off()
  
  # Plot 3: Biomarker coverage heatmap
  message("Creating biomarker coverage heatmap...")
  
  # Create biomarker presence matrix
  biomarker_genes <- unique(biomarker_data$gene_symbol)
  sources <- unique(plasma_data$source)
  
  coverage_matrix <- expand.grid(gene_symbol = biomarker_genes, 
                               source = sources, 
                               stringsAsFactors = FALSE)
  
  coverage_matrix$detected <- apply(coverage_matrix, 1, function(row) {
    any(biomarker_data$gene_symbol == row[1] & biomarker_data$source == row[2])
  })
  
  # Add abundance values
  coverage_matrix$abundance <- NA
  for(i in 1:nrow(coverage_matrix)) {
    matches <- biomarker_data[biomarker_data$gene_symbol == coverage_matrix$gene_symbol[i] & 
                            biomarker_data$source == coverage_matrix$source[i], ]
    if(nrow(matches) > 0) {
      coverage_matrix$abundance[i] <- matches$log_abundance[1]
    }
  }
  
  tiff(file.path(output_dir, "biomarker_coverage_heatmap.tiff"), 
       width = 12, height = 14, units = "in", res = 600, compression = "lzw")
  
  p3 <- ggplot(coverage_matrix, aes(x = source, y = gene_symbol, fill = abundance)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(name = "Log10\n(Abundance)", na.value = "grey90") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8)) +
    labs(title = "Biomarker Gene Coverage Across Plasma Sources",
         subtitle = "Heatmap showing detection and abundance levels",
         x = "Data Source", y = "Biomarker Gene")
  
  print(p3)
  dev.off()
  
  # Plot 4: Top biomarkers by detection frequency
  message("Creating top biomarkers plot...")
  
  biomarker_freq <- biomarker_data %>%
    group_by(gene_symbol) %>%
    summarise(
      n_sources = n_distinct(source),
      mean_abundance = mean(abundance, na.rm = TRUE),
      mean_log_abundance = mean(log_abundance, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_sources), desc(mean_log_abundance)) %>%
    head(20)
  
  tiff(file.path(output_dir, "top_biomarkers_detection.tiff"), 
       width = 12, height = 8, units = "in", res = 600, compression = "lzw")
  
  p4 <- ggplot(biomarker_freq, aes(x = reorder(gene_symbol, n_sources), 
                                  y = n_sources, fill = mean_log_abundance)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_viridis_c(name = "Mean Log10\n(Abundance)") +
    theme_minimal() +
    labs(title = "Top 20 Biomarkers by Detection Frequency",
         subtitle = "Number of plasma sources detecting each biomarker",
         x = "Biomarker Gene", y = "Number of Sources Detected")
  
  print(p4)
  dev.off()
  
  # Plot 5: Source-specific biomarker rankings
  message("Creating source-specific rankings...")
  
  # Calculate percentiles for each source
  percentile_data <- plasma_data %>%
    group_by(source) %>%
    mutate(abundance_percentile = percent_rank(log_abundance)) %>%
    filter(is_biomarker) %>%
    ungroup()
  
  tiff(file.path(output_dir, "biomarker_percentile_rankings.tiff"), 
       width = 14, height = 10, units = "in", res = 600, compression = "lzw")
  
  p5 <- ggplot(percentile_data, aes(x = source, y = abundance_percentile, 
                                   color = source)) +
    geom_beeswarm(alpha = 0.7, size = 2) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    scale_color_viridis_d(name = "Source") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Biomarker Gene Abundance Percentiles by Source",
         subtitle = "Where biomarkers rank within each source's expression distribution",
         x = "Data Source", y = "Abundance Percentile") +
    scale_y_continuous(labels = scales::percent)
  
  print(p5)
  dev.off()
}

# Main analysis function
main <- function() {
  message("Starting biomarker plasma analysis...")
  
  # Load data
  biomarker_genes <- load_biomarker_genes()
  plasma_data <- load_plasma_data()
  
  # Analyze biomarkers
  analysis_results <- analyze_biomarkers_in_plasma(plasma_data, biomarker_genes)
  
  # Create plots
  create_biomarker_plots(analysis_results)
  
  # Save results
  write.csv(analysis_results$plasma_data, 
           file.path(output_dir, "plasma_expression_data.csv"), 
           row.names = FALSE)
  
  write.csv(analysis_results$biomarker_data, 
           file.path(output_dir, "biomarker_expression_data.csv"), 
           row.names = FALSE)
  
  write.csv(analysis_results$summary, 
           file.path(output_dir, "biomarker_summary_statistics.csv"), 
           row.names = FALSE)
  
  message("Biomarker analysis completed!")
  message(paste("Results saved to:", output_dir))
}

# Run analysis
main() 