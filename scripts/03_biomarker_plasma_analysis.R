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
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "ggthemes", "patchwork", "UpSetR")
load_packages(required_packages)

# Parse command line arguments (simple approach)
args <- commandArgs(trailingOnly = TRUE)
force_mapping <- "--force-mapping" %in% args

# Load gene mapping utility
source("scripts/data_processing/simple_id_mapping.R")

# Load gene deduplication utility
source("scripts/utilities/gene_deduplication.R")

# Set output directory
output_dir <- get_output_path("03_biomarker_plasma_analysis", subdir = "plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read biomarker gene list
biomarker_file <- "data/metadata/biomarkers_list.tsv"
biomarkers <- read_tsv(biomarker_file, show_col_types = FALSE)
biomarker_genes <- unique(biomarkers$To)

# Function to calculate z-score normalization within each database
calculate_zscore_normalization_biomarker <- function(data, intensity_col) {
  data %>%
    mutate(
      log_intensity = log10(.data[[intensity_col]] + 1),
      z_score = (log_intensity - mean(log_intensity, na.rm = TRUE)) / sd(log_intensity, na.rm = TRUE)
    )
}

# Helper function to create distribution and violin plots for each database (with z-score normalization)
create_database_plots <- function(data, gene_col, intensity_col, database_name, biomarker_genes) {
  # Filter out NA and infinite values and apply z-score normalization
  plot_data <- data %>%
    filter(!is.na(.data[[gene_col]]), !is.na(.data[[intensity_col]])) %>%
    calculate_zscore_normalization_biomarker(intensity_col) %>%
    mutate(
      is_biomarker = .data[[gene_col]] %in% biomarker_genes
    )
  
  # Calculate statistics for both log and z-score values
  stats_log <- plot_data %>%
    group_by(is_biomarker) %>%
    summarise(
      n = n(),
      median = median(log_intensity, na.rm = TRUE),
      q1 = quantile(log_intensity, 0.25, na.rm = TRUE),
      q3 = quantile(log_intensity, 0.75, na.rm = TRUE),
      .groups = 'drop'
    )
  
  stats_zscore <- plot_data %>%
    group_by(is_biomarker) %>%
    summarise(
      n = n(),
      median = median(z_score, na.rm = TRUE),
      q1 = quantile(z_score, 0.25, na.rm = TRUE),
      q3 = quantile(z_score, 0.75, na.rm = TRUE),
      .groups = 'drop'
    )
  
  biomarker_count <- sum(plot_data$is_biomarker)
  total_count <- nrow(plot_data)
  percentage <- round(biomarker_count / length(biomarker_genes) * 100, 1)
  
  # Distribution plot (log-transformed)
  dist_plot <- ggplot(plot_data) +
    geom_density(aes(x = log_intensity, fill = is_biomarker), alpha = 0.5) +
    scale_fill_manual(values = c("FALSE" = "#69b3a2", "TRUE" = "#E69F00"),
                     labels = c("All proteins", "Biomarkers")) +
    labs(
      title = paste(database_name, "Protein Distribution (Log10)"),
      subtitle = sprintf("Detected biomarkers: %d (%.1f%%)", biomarker_count, percentage),
      x = "log10(Intensity + 1)",
      y = "Density",
      fill = "Protein Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  # Z-score distribution plot
  dist_plot_zscore <- ggplot(plot_data) +
    geom_density(aes(x = z_score, fill = is_biomarker), alpha = 0.5) +
    scale_fill_manual(values = c("FALSE" = "#69b3a2", "TRUE" = "#E69F00"),
                     labels = c("All proteins", "Biomarkers")) +
    labs(
      title = paste(database_name, "Protein Distribution (Z-Score)"),
      subtitle = sprintf("Z-score normalized for better comparison | Biomarkers: %d (%.1f%%)", biomarker_count, percentage),
      x = "Z-Score (standardized intensity)",
      y = "Density",
      fill = "Protein Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  # Violin plot (log-transformed)
  violin_plot <- ggplot(plot_data, aes(x = is_biomarker, y = log_intensity, fill = is_biomarker)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    scale_fill_manual(values = c("FALSE" = "#69b3a2", "TRUE" = "#E69F00"),
                     labels = c("All proteins", "Biomarkers")) +
    scale_x_discrete(labels = c("FALSE" = "All proteins", "TRUE" = "Biomarkers")) +
    labs(
      title = paste(database_name, "Abundance Distribution (Log10)"),
      subtitle = sprintf("Median (All): %.2f, Median (Bio): %.2f", 
                        stats_log$median[1], stats_log$median[2]),
      x = "",
      y = "log10(Intensity + 1)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "none"
    )
  
  # Z-score violin plot
  violin_plot_zscore <- ggplot(plot_data, aes(x = is_biomarker, y = z_score, fill = is_biomarker)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    scale_fill_manual(values = c("FALSE" = "#69b3a2", "TRUE" = "#E69F00"),
                     labels = c("All proteins", "Biomarkers")) +
    scale_x_discrete(labels = c("FALSE" = "All proteins", "TRUE" = "Biomarkers")) +
    labs(
      title = paste(database_name, "Abundance Distribution (Z-Score)"),
      subtitle = sprintf("Z-score normalized | Median (All): %.2f, Median (Bio): %.2f", 
                        stats_zscore$median[1], stats_zscore$median[2]),
      x = "",
      y = "Z-Score (standardized intensity)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "none"
    )
  
  return(list(
    distribution = dist_plot, 
    violin = violin_plot,
    distribution_zscore = dist_plot_zscore,
    violin_zscore = violin_plot_zscore
  ))
}

# Process each database
message("[Biomarker Analysis] Processing databases...")

# 1. PeptideAtlas
message("Processing PeptideAtlas...")
peptideatlas_raw <- read_csv("data/raw/peptideatlas/peptideatlas.csv", show_col_types = FALSE)
peptideatlas_raw$gene <- convert_to_gene_symbol(peptideatlas_raw$biosequence_accession, force_mapping = force_mapping)

# Deduplicate genes before analysis
peptideatlas <- deduplicate_genes(peptideatlas_raw, "gene", "norm_PSMs_per_100K", 
                                additional_cols = c("biosequence_accession"), 
                                aggregation_method = "median")

peptideatlas_plots <- create_database_plots(peptideatlas, "gene", "norm_PSMs_per_100K", "PeptideAtlas", biomarker_genes)
peptideatlas_biomarkers <- unique(peptideatlas$gene[peptideatlas$gene %in% biomarker_genes])

# 2. HPA MS
message("Processing HPA MS...")
hpa_ms_raw <- read_csv("data/raw/hpa/hpa_ms.csv", show_col_types = FALSE, skip = 1)

# Deduplicate genes before analysis
hpa_ms <- deduplicate_genes(hpa_ms_raw, "Gene", "Concentration", aggregation_method = "median")

hpa_ms_plots <- create_database_plots(hpa_ms, "Gene", "Concentration", "HPA MS", biomarker_genes)
hpa_ms_biomarkers <- unique(hpa_ms$Gene[hpa_ms$Gene %in% biomarker_genes])

# 3. HPA PEA
message("Processing HPA PEA...")
hpa_pea_raw <- read_csv("data/raw/hpa/hpa_pea.csv", show_col_types = FALSE)

# Deduplicate genes before analysis
hpa_pea <- deduplicate_genes(hpa_pea_raw, "Gene", "median_npx", aggregation_method = "median")

hpa_pea_plots <- create_database_plots(hpa_pea, "Gene", "median_npx", "HPA PEA", biomarker_genes)
hpa_pea_biomarkers <- unique(hpa_pea$Gene[hpa_pea$Gene %in% biomarker_genes])

# 4. HPA Immunoassay
message("Processing HPA Immunoassay...")
hpa_imm_raw <- read_csv("data/raw/hpa/hpa_immunoassay_plasma.csv", show_col_types = FALSE)

# Deduplicate genes before analysis
hpa_imm <- deduplicate_genes(hpa_imm_raw, "Gene", "Concentration", aggregation_method = "median")

hpa_imm_plots <- create_database_plots(hpa_imm, "Gene", "Concentration", "HPA Immunoassay", biomarker_genes)
hpa_imm_biomarkers <- unique(hpa_imm$Gene[hpa_imm$Gene %in% biomarker_genes])

# 5. GPMDB
message("Processing GPMDB...")
gpmdb_raw <- read_csv("data/raw/gpmdb/gpmdb_plasma.csv", show_col_types = FALSE)
gpmdb_raw$gene <- stringr::str_extract(gpmdb_raw$description, "[A-Z0-9]+(?=,| |$)")

# Deduplicate genes before analysis
if ("total" %in% colnames(gpmdb_raw)) {
  gpmdb <- deduplicate_genes(gpmdb_raw, "gene", "total", aggregation_method = "median")
} else {
  # Just remove duplicate genes if no quantification column
  gpmdb <- gpmdb_raw %>% 
    filter(!is.na(gene) & gene != "") %>%
    distinct(gene, .keep_all = TRUE)
}

gpmdb_plots <- create_database_plots(gpmdb, "gene", "total", "GPMDB", biomarker_genes)
gpmdb_biomarkers <- unique(gpmdb$gene[gpmdb$gene %in% biomarker_genes])

# 6. PAXDB
message("Processing PAXDB...")
paxdb_raw <- read_csv("data/raw/paxdb/paxdb_plasma.csv", show_col_types = FALSE)
paxdb_raw$ensp <- stringr::str_replace(paxdb_raw$string_external_id, "^9606\\.", "")
paxdb_raw$gene <- convert_to_gene_symbol(paxdb_raw$ensp, force_mapping = force_mapping)

# Deduplicate genes before analysis
if ("abundance" %in% colnames(paxdb_raw)) {
  paxdb <- deduplicate_genes(paxdb_raw, "gene", "abundance", 
                           additional_cols = c("ensp"), 
                           aggregation_method = "median")
} else {
  # Just remove duplicate genes if no quantification column
  paxdb <- paxdb_raw %>% 
    filter(!is.na(gene) & gene != "") %>%
    distinct(gene, .keep_all = TRUE)
}

paxdb_plots <- create_database_plots(paxdb, "gene", "abundance", "PAXDB", biomarker_genes)
paxdb_biomarkers <- unique(paxdb$gene[paxdb$gene %in% biomarker_genes])

# Create comprehensive plots
message("[Biomarker Analysis] Creating comprehensive plots...")

# Log-transformed distribution plots combined
distribution_combined <- (peptideatlas_plots$distribution + hpa_ms_plots$distribution) /
  (hpa_pea_plots$distribution + hpa_imm_plots$distribution) /
  (gpmdb_plots$distribution + paxdb_plots$distribution) +
  plot_annotation(
    title = "Protein Abundance Distribution Across Databases (Log10-transformed)",
    subtitle = "Comparing all proteins vs. biomarkers - Log10 scale",
    tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)')),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave(
  file.path(output_dir, "distributions_by_database.png"),
  distribution_combined,
  width = 15,
  height = 20,
  dpi = 300
)

# Z-score normalized distribution plots combined
distribution_combined_zscore <- (peptideatlas_plots$distribution_zscore + hpa_ms_plots$distribution_zscore) /
  (hpa_pea_plots$distribution_zscore + hpa_imm_plots$distribution_zscore) /
  (gpmdb_plots$distribution_zscore + paxdb_plots$distribution_zscore) +
  plot_annotation(
    title = "Protein Abundance Distribution Across Databases (Z-Score Normalized)",
    subtitle = "Comparing all proteins vs. biomarkers - Z-score normalized for better cross-source comparison",
    tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)')),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave(
  file.path(output_dir, "distributions_by_database_zscore.png"),
  distribution_combined_zscore,
  width = 15,
  height = 20,
  dpi = 300
)

# Log-transformed violin plots combined
violin_combined <- (peptideatlas_plots$violin + hpa_ms_plots$violin) /
  (hpa_pea_plots$violin + hpa_imm_plots$violin) /
  (gpmdb_plots$violin + paxdb_plots$violin) +
  plot_annotation(
    title = "Protein Abundance Distribution Across Databases (Log10-transformed)",
    subtitle = "Comparing all proteins vs. biomarkers - Log10 scale",
    tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)')),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave(
  file.path(output_dir, "abundance_distributions_by_database.png"),
  violin_combined,
  width = 15,
  height = 20,
  dpi = 300
)

# Z-score normalized violin plots combined
violin_combined_zscore <- (peptideatlas_plots$violin_zscore + hpa_ms_plots$violin_zscore) /
  (hpa_pea_plots$violin_zscore + hpa_imm_plots$violin_zscore) /
  (gpmdb_plots$violin_zscore + paxdb_plots$violin_zscore) +
  plot_annotation(
    title = "Protein Abundance Distribution Across Databases (Z-Score Normalized)",
    subtitle = "Comparing all proteins vs. biomarkers - Z-score normalized for better cross-source comparison",
    tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)')),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave(
  file.path(output_dir, "abundance_distributions_by_database_zscore.png"),
  violin_combined_zscore,
  width = 15,
  height = 20,
  dpi = 300
)

# Create presence/absence matrix for biomarkers
biomarker_matrix <- data.frame(
  Gene = biomarker_genes,
  PeptideAtlas = as.integer(biomarker_genes %in% peptideatlas_biomarkers),
  HPA_MS = as.integer(biomarker_genes %in% hpa_ms_biomarkers),
  HPA_PEA = as.integer(biomarker_genes %in% hpa_pea_biomarkers),
  HPA_Immunoassay = as.integer(biomarker_genes %in% hpa_imm_biomarkers),
  GPMDB = as.integer(biomarker_genes %in% gpmdb_biomarkers),
  PAXDB = as.integer(biomarker_genes %in% paxdb_biomarkers)
)

# Create summary of biomarker counts
source_summary <- data.frame(
  Source = c("PeptideAtlas", "HPA MS", "HPA PEA", "HPA Immunoassay", "GPMDB", "PAXDB"),
  Biomarkers = c(
    length(peptideatlas_biomarkers),
    length(hpa_ms_biomarkers),
    length(hpa_pea_biomarkers),
    length(hpa_imm_biomarkers),
    length(gpmdb_biomarkers),
    length(paxdb_biomarkers)
  ),
  Total_Biomarkers = length(biomarker_genes)
)
source_summary$Percentage <- round(source_summary$Biomarkers / source_summary$Total_Biomarkers * 100, 1)

# Create and save UpSet plot
message("[Biomarker Analysis] Creating technology overlap plot...")
png(file.path(output_dir, "biomarker_technology_overlap.png"), 
    width = 12, height = 8, units = "in", res = 300)
upset(
  biomarker_matrix,
  sets = c("PeptideAtlas", "HPA_MS", "HPA_PEA", "HPA_Immunoassay", "GPMDB", "PAXDB"),
  order.by = "freq",
  nsets = 6,
  nintersects = 30,
  point.size = 3,
  line.size = 1,
  mainbar.y.label = "Biomarker Intersections",
  sets.x.label = "Biomarkers per Technology",
  text.scale = 1.2,
  mb.ratio = c(0.6, 0.4),
  keep.order = TRUE,
  main.bar.color = "#69b3a2",
  sets.bar.color = "#E69F00"
)
dev.off()

# Create summary table
write_csv(source_summary, file.path(output_dir, "biomarker_detection_summary.csv"))

# Print summary
message("\nBiomarker Detection Summary:")
message("==========================")
for(i in 1:nrow(source_summary)) {
  message(sprintf("%s: %d biomarkers (%.1f%%)",
                 source_summary$Source[i],
                 source_summary$Biomarkers[i],
                 source_summary$Percentage[i]))
}

message("[Biomarker Analysis] Completed successfully.")

message("\nAnalysis complete. Files saved to: ", output_dir)

# Create distribution plot for PeptideAtlas (as representative dataset)
# Note: After deduplication, we use the deduplicated data
distribution_plot <- ggplot() +
  geom_density(data = data.frame(value = log10(peptideatlas$norm_PSMs_per_100K[!is.na(peptideatlas$norm_PSMs_per_100K)])), 
               aes(x = value), fill = "#69b3a2", alpha = 0.3) +
  geom_density(data = data.frame(value = log10(peptideatlas$norm_PSMs_per_100K[peptideatlas$gene %in% biomarker_genes & !is.na(peptideatlas$norm_PSMs_per_100K)])), 
               aes(x = value), fill = "#E69F00", alpha = 0.3) +
  labs(
    title = "Distribution of Protein Intensities (Gene-Deduplicated)",
    subtitle = "All genes vs. Biomarker genes (using median values for duplicate genes)",
    x = "log10(PSMs per 100K)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )