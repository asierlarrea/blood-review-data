#!/usr/bin/env Rscript
# ============================================================================
# Biomarker Plasma Analysis (Refined)
# Description: Plots expression distributions for each plasma source, overlays lines for each biomarker gene
# Output: Plots in outputs/plots/biomarker_plasma_analysis
# ============================================================================

# Load utilities and set up output directories
source("scripts/utilities/load_packages.R")
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "ggthemes", "patchwork")
load_packages(required_packages)

# Parse command line arguments (simple approach)
args <- commandArgs(trailingOnly = TRUE)
force_mapping <- "--force-mapping" %in% args

# Load gene mapping utility
source("scripts/data_processing/simple_id_mapping.R")

# Set output directory
output_dir <- get_output_path("biomarker_plasma_analysis", subdir = "plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read biomarker gene list
biomarker_file <- "data/metadata/biomarkers_list.tsv"
biomarkers <- read_tsv(biomarker_file, show_col_types = FALSE)
biomarker_genes <- unique(biomarkers$To)

# Helper: plot distribution with biomarker lines
distribution_plot <- function(df, gene_col, value_col, source_name, log10=TRUE, biomarker_genes, output_dir) {
  # Data preparation
  df <- df %>% filter(!is.na(.data[[gene_col]]), !is.na(.data[[value_col]]))
  if (log10) {
    df <- df %>% mutate(expr = log10(.data[[value_col]] + 1))
    xlab <- paste0("log10(", value_col, "+1)")
  } else {
    df <- df %>% mutate(expr = .data[[value_col]])
    xlab <- value_col
  }
  
  # Biomarker values present in this source
  biomarker_vals <- df %>% filter(.data[[gene_col]] %in% biomarker_genes)
  
  # Calculate statistics
  total_genes <- nrow(df)
  biomarker_count <- nrow(biomarker_vals)
  biomarker_percentage <- round(biomarker_count / total_genes * 100, 1)
  
  # Calculate quantiles for biomarker genes
  biomarker_stats <- biomarker_vals %>%
    summarise(
      min = min(expr, na.rm = TRUE),
      q25 = quantile(expr, 0.25, na.rm = TRUE),
      median = median(expr, na.rm = TRUE),
      q75 = quantile(expr, 0.75, na.rm = TRUE),
      max = max(expr, na.rm = TRUE)
    )
  
  # Create main plot
  p <- ggplot(df, aes(x = expr)) +
    # Add density plot with gradient fill
    geom_density(aes(y = ..density..), fill = "#69b3a2", alpha = 0.3, color = NA) +
    # Add histogram with transparency
    geom_histogram(aes(y = ..density..), bins = 60, fill = "#69b3a2", alpha = 0.4) +
    # Add density line
    geom_density(color = "#1a1a1a", size = 0.8) +
    # Add biomarker gene lines with annotations
    geom_vline(data = biomarker_vals, aes(xintercept = expr), 
               color = "red", linetype = "dashed", alpha = 0.5) +
    # Add biomarker gene labels
    geom_text(data = biomarker_vals, 
              aes(x = expr, y = Inf, label = .data[[gene_col]]),
              angle = 90, hjust = -0.2, vjust = 0.5, size = 3, color = "red") +
    # Add theme and styling
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    # Add labels
    labs(
      title = paste0(source_name, ": Expression Distribution"),
      subtitle = sprintf("Biomarker genes: %d (%.1f%% of total genes)", biomarker_count, biomarker_percentage),
      x = xlab,
      y = "Density"
    )
  
  # Add statistics annotation
  p <- p + annotate(
    "text",
    x = Inf, y = Inf,
    label = sprintf(
      "Biomarker Statistics:\nMin: %.2f\nQ1: %.2f\nMedian: %.2f\nQ3: %.2f\nMax: %.2f",
      biomarker_stats$min, biomarker_stats$q25, biomarker_stats$median,
      biomarker_stats$q75, biomarker_stats$max
    ),
    hjust = 1.1, vjust = 1.1,
    size = 3.5,
    color = "darkred"
  )
  
  # Save plot
  plot_file <- file.path(output_dir, paste0("expression_distribution_", str_replace_all(tolower(source_name), " ", "_"), ".png"))
  ggsave(plot_file, p, width = 12, height = 8, dpi = 300, bg = "white")
  message(paste("[Biomarker Analysis] Plot saved to:", plot_file))
  
  # Return the plot object for potential further modifications
  return(p)
}

# 1. PeptideAtlas
message("[Biomarker Analysis] Processing PeptideAtlas...")
peptideatlas <- read_csv("data/raw/peptideatlas/peptideatlas.csv", show_col_types = FALSE)
# Map accessions to gene symbols
peptideatlas$gene <- convert_to_gene_symbol(peptideatlas$biosequence_accession, force_mapping = force_mapping)
# Use norm_PSMs_per_100K as expression value
peptideatlas <- peptideatlas %>% filter(!is.na(norm_PSMs_per_100K))
p1 <- distribution_plot(peptideatlas, "gene", "norm_PSMs_per_100K", "PeptideAtlas", TRUE, biomarker_genes, output_dir)

# 2. HPA MS
message("[Biomarker Analysis] Processing HPA MS...")
hpa_ms <- read_csv("data/raw/hpa/hpa_ms.csv", show_col_types = FALSE, skip = 1)
hpa_ms <- hpa_ms %>% rename(gene = Gene, expr = Concentration)
p2 <- distribution_plot(hpa_ms, "gene", "expr", "HPA MS", TRUE, biomarker_genes, output_dir)

# 3. HPA PEA
message("[Biomarker Analysis] Processing HPA PEA...")
hpa_pea <- read_csv("data/raw/hpa/hpa_pea.csv", show_col_types = FALSE)
hpa_pea <- hpa_pea %>% rename(gene = Gene, expr = `Variation between individuals`)
p3 <- distribution_plot(hpa_pea, "gene", "expr", "HPA PEA", TRUE, biomarker_genes, output_dir)

# 4. HPA Immunoassay
message("[Biomarker Analysis] Processing HPA Immunoassay...")
hpa_imm <- read_csv("data/raw/hpa/hpa_immunoassay_plasma.csv", show_col_types = FALSE)
hpa_imm <- hpa_imm %>% rename(gene = Gene, expr = Concentration)
p4 <- distribution_plot(hpa_imm, "gene", "expr", "HPA Immunoassay", TRUE, biomarker_genes, output_dir)

# 5. GPMDB
message("[Biomarker Analysis] Processing GPMDB...")
gpmdb <- read_csv("data/raw/gpmdb/gpmdb_plasma.csv", show_col_types = FALSE)
gpmdb$gene <- stringr::str_extract(gpmdb$description, "[A-Z0-9]+(?=,| |$)")
p5 <- distribution_plot(gpmdb, "gene", "total", "GPMDB", TRUE, biomarker_genes, output_dir)

# 6. PAXDB
message("[Biomarker Analysis] Processing PAXDB...")
paxdb <- read_csv("data/raw/paxdb/paxdb_plasma.csv", show_col_types = FALSE)
paxdb$ensp <- stringr::str_replace(paxdb$string_external_id, "^9606\\.", "")
paxdb$gene <- convert_to_gene_symbol(paxdb$ensp, force_mapping = force_mapping)
p6 <- distribution_plot(paxdb, "gene", "abundance", "PAXDB", TRUE, biomarker_genes, output_dir)

# Create a combined plot
message("[Biomarker Analysis] Creating combined plot...")
combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_annotation(
    title = "Biomarker Gene Expression Distributions Across Data Sources",
    subtitle = "Red dashed lines indicate biomarker genes",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  )

# Save combined plot
ggsave(
  file.path(output_dir, "combined_biomarker_distributions.png"),
  combined_plot,
  width = 20,
  height = 24,
  dpi = 300,
  bg = "white"
)

message("[Biomarker Analysis] Completed successfully.") 