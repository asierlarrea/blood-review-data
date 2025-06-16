#!/usr/bin/env Rscript
# ============================================================================
# Biomarker Plasma Analysis (Refined)
# Description: Plots expression distributions for each plasma source, overlays lines for each biomarker gene
# Output: Plots in outputs/plots/biomarker_plasma_analysis
# ============================================================================

# Function to create a waterfall plot for a specific database
create_waterfall_plot <- function(data, database_name, intensity_col, gene_col = "gene") {
  plot_data <- data %>%
    mutate(
      log_intensity = log10(.data[[intensity_col]] + 1),
      z_score = (log_intensity - mean(log_intensity, na.rm = TRUE)) / sd(log_intensity, na.rm = TRUE),
      gene = .data[[gene_col]],
      is_highlight = gene %in% highlight_biomarkers,
      is_biomarker = gene %in% biomarker_genes,
      rank = rank(z_score, ties.method = "first")
    ) %>%
    arrange(z_score)
  
  # Calculate mean and SD for annotations
  stats <- plot_data %>%
    summarise(
      biomarker_mean = mean(z_score[is_biomarker], na.rm = TRUE),
      all_mean = mean(z_score, na.rm = TRUE),
      n_proteins = n(),
      n_biomarkers = sum(is_biomarker)
    )
  
  # Create annotation text
  annotation_text <- sprintf(
    "Total proteins: %d\nBiomarkers: %d (%.1f%%)",
    stats$n_proteins,
    stats$n_biomarkers,
    100 * stats$n_biomarkers / stats$n_proteins
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = rank, y = z_score)) +
    # Add background for the plot
    theme_bw() +
    # Add background shading for different z-score regions
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2, ymax = 2,
             fill = "grey95", alpha = 0.5) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1, ymax = 1,
             fill = "grey90", alpha = 0.5) +
    # Add horizontal lines for z-score reference
    geom_hline(yintercept = c(-2, -1, 0, 1, 2), 
               color = "grey80", linetype = "dashed") +
    # Add the main bars
    geom_col(aes(fill = is_biomarker), width = 1) +
    # Add points and labels for highlighted biomarkers
    geom_point(data = filter(plot_data, is_highlight),
              color = "#FF4B4B", size = 2) +
    geom_text_repel(
      data = filter(plot_data, is_highlight),
      aes(label = gene),
      direction = "both",
      nudge_y = 0.5,
      size = 3,
      color = "#FF4B4B",
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.5,
      segment.color = "#FF4B4B",
      segment.alpha = 0.5
    ) +
    # Add mean lines
    geom_hline(yintercept = stats$biomarker_mean, 
               color = "#E69F00", 
               linetype = "dashed") +
    geom_hline(yintercept = stats$all_mean, 
               color = "#69b3a2", 
               linetype = "dashed") +
    # Customize colors and labels
    scale_fill_manual(
      values = c("TRUE" = "#E69F00", "FALSE" = "#69b3a2"),
      name = "Protein Type",
      labels = c("All proteins", "Biomarkers")
    ) +
    # Add titles and labels
    labs(
      title = paste(database_name, "Protein Abundance Distribution"),
      subtitle = "Z-score normalized abundances with highlighted key biomarkers",
      x = "Proteins (ranked by abundance)",
      y = "Z-Score"
    ) +
    # Add annotation text
    annotate(
      "text",
      x = max(plot_data$rank) * 0.95,
      y = max(plot_data$z_score) * 0.8,
      label = annotation_text,
      hjust = 1,
      size = 3,
      color = "grey30"
    ) +
    # Customize theme
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

# Load utilities and set up output directories
source("scripts/utilities/load_packages.R")
source("scripts/config/analysis_config.R")
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "ggthemes", "patchwork", "UpSetR", "ggrepel")
load_packages(required_packages)

# Parse command line arguments (simple approach)
args <- commandArgs(trailingOnly = TRUE)
force_mapping <- "--force-mapping" %in% args

# Parse formats argument
formats_arg <- NULL
if (length(args) > 0) {
  formats_idx <- which(args == "--formats")
  if (length(formats_idx) > 0 && formats_idx < length(args)) {
    formats_arg <- args[formats_idx + 1]
  }
}

# Set plot formats if provided
if (!is.null(formats_arg)) {
  set_plot_formats(formats_arg)
}

# Load gene mapping utility
source("scripts/data_processing/simple_id_mapping.R")

# Load gene deduplication utility
source("scripts/utilities/gene_deduplication.R")

# Load quantile normalization utility
source("scripts/utilities/quantile_normalization_functions.R")

# Set output directory
output_dir <- get_output_path("03_biomarker_plasma_analysis", subdir = "plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read biomarker gene list
biomarker_file <- "/Users/yperez/work/blood-review-data/data/metadata/biomarkers_list.csv"
biomarkers <- read_csv(biomarker_file, show_col_types = FALSE)
biomarker_genes <- unique(biomarkers$gene_name)

# Define biomarkers to highlight
highlight_biomarkers <- c("APOE", "APOB", "CRP", "SAA1", "TNF", "IL1B", "IFNG", "IL10", "ALB", "INS", "NPPB")

# Function to calculate z-score normalization within each database
calculate_zscore_normalization_biomarker <- function(data, intensity_col) {
  data %>%
    mutate(
      log_intensity = log10(.data[[intensity_col]] + 1),
      z_score = (log_intensity - mean(log_intensity, na.rm = TRUE)) / sd(log_intensity, na.rm = TRUE)
    )
}

# Helper function to create distribution and violin plots for each database (with z-score normalization)
create_database_plots <- function(data, gene_col, intensity_col, database_name, biomarker_genes, highlight_biomarkers) {
  # Filter out NA and infinite values and apply all normalization methods
  plot_data <- data %>%
    filter(!is.na(.data[[gene_col]]), !is.na(.data[[intensity_col]])) %>%
    calculate_zscore_normalization_biomarker(intensity_col) %>%
    mutate(
      is_biomarker = .data[[gene_col]] %in% biomarker_genes,
      is_highlight = .data[[gene_col]] %in% highlight_biomarkers
    )
  
  # Z-score distribution plot with highlighted biomarkers
  dist_plot_zscore <- ggplot(plot_data) +
    geom_density(aes(x = z_score, fill = is_biomarker), alpha = 0.5) +
    geom_vline(data = filter(plot_data, is_highlight), 
               aes(xintercept = z_score), 
               linetype = "dashed", 
               color = "#FF4B4B",
               alpha = 0.6) +
    geom_text_repel(
      data = filter(plot_data, is_highlight),
      aes(x = z_score, y = 0, label = .data[[gene_col]]),
      direction = "both",
      nudge_y = 0.1,
      size = 3,
      color = "#FF4B4B",
      max.overlaps = Inf
    ) +
    scale_fill_manual(values = c("FALSE" = "#69b3a2", "TRUE" = "#E69F00"),
                     labels = c("All proteins", "Biomarkers")) +
    labs(
      title = paste(database_name, "Protein Distribution (Z-Score)"),
      subtitle = sprintf("Z-score normalized with highlighted key biomarkers"),
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
  
  # Z-score violin plot with highlighted biomarkers
  violin_plot_zscore <- ggplot(plot_data, aes(x = is_biomarker, y = z_score, fill = is_biomarker)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    geom_point(data = filter(plot_data, is_highlight),
               aes(x = TRUE, y = z_score),
               color = "#FF4B4B",
               size = 2) +
    geom_text_repel(
      data = filter(plot_data, is_highlight),
      aes(x = TRUE, y = z_score, label = .data[[gene_col]]),
      direction = "y",
      nudge_x = 0.2,
      size = 3,
      color = "#FF4B4B",
      max.overlaps = Inf
    ) +
    scale_fill_manual(values = c("FALSE" = "#69b3a2", "TRUE" = "#E69F00"),
                     labels = c("All proteins", "Biomarkers")) +
    scale_x_discrete(labels = c("FALSE" = "All proteins", "TRUE" = "Biomarkers")) +
    labs(
      title = paste(database_name, "Abundance Distribution (Z-Score)"),
      subtitle = "Z-score normalized with highlighted key biomarkers",
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

peptideatlas_plots <- create_database_plots(peptideatlas, "gene", "norm_PSMs_per_100K", "PeptideAtlas", biomarker_genes, highlight_biomarkers)
peptideatlas_biomarkers <- unique(peptideatlas$gene[peptideatlas$gene %in% biomarker_genes])

# 2. HPA MS
message("Processing HPA MS...")
hpa_ms_raw <- read_csv("data/raw/hpa/hpa_ms.csv", show_col_types = FALSE, skip = 1)

# Deduplicate genes before analysis
hpa_ms <- deduplicate_genes(hpa_ms_raw, "Gene", "Concentration", aggregation_method = "median")

hpa_ms_plots <- create_database_plots(hpa_ms, "Gene", "Concentration", "HPA MS", biomarker_genes, highlight_biomarkers)
hpa_ms_biomarkers <- unique(hpa_ms$Gene[hpa_ms$Gene %in% biomarker_genes])

# 3. HPA PEA
message("Processing HPA PEA...")
hpa_pea_raw <- read_csv("data/raw/hpa/hpa_pea.csv", show_col_types = FALSE)

# Deduplicate genes before analysis
hpa_pea <- deduplicate_genes(hpa_pea_raw, "Gene", "median_npx", aggregation_method = "median")

hpa_pea_plots <- create_database_plots(hpa_pea, "Gene", "median_npx", "HPA PEA", biomarker_genes, highlight_biomarkers)
hpa_pea_biomarkers <- unique(hpa_pea$Gene[hpa_pea$Gene %in% biomarker_genes])

# 4. HPA Immunoassay
message("Processing HPA Immunoassay...")
hpa_imm_raw <- read_csv("data/raw/hpa/hpa_immunoassay_plasma.csv", show_col_types = FALSE)

# Deduplicate genes before analysis
hpa_imm <- deduplicate_genes(hpa_imm_raw, "Gene", "Concentration", aggregation_method = "median")

hpa_imm_plots <- create_database_plots(hpa_imm, "Gene", "Concentration", "HPA Immunoassay", biomarker_genes, highlight_biomarkers)
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

gpmdb_plots <- create_database_plots(gpmdb, "gene", "total", "GPMDB", biomarker_genes, highlight_biomarkers)
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

paxdb_plots <- create_database_plots(paxdb, "gene", "abundance", "PAXDB", biomarker_genes, highlight_biomarkers)
paxdb_biomarkers <- unique(paxdb$gene[paxdb$gene %in% biomarker_genes])

# Create waterfall plots for each database
message("[Biomarker Analysis] Creating waterfall plots...")

waterfall_plots <- list(
  peptideatlas = create_waterfall_plot(peptideatlas, "PeptideAtlas", "norm_PSMs_per_100K"),
  hpa_ms = create_waterfall_plot(hpa_ms, "HPA MS", "Concentration", "Gene"),
  hpa_pea = create_waterfall_plot(hpa_pea, "HPA PEA", "median_npx", "Gene"),
  hpa_imm = create_waterfall_plot(hpa_imm, "HPA Immunoassay", "Concentration", "Gene"),
  gpmdb = create_waterfall_plot(gpmdb, "GPMDB", "total"),
  paxdb = create_waterfall_plot(paxdb, "PAXDB", "abundance")
)

# Save individual waterfall plots
# for (db_name in names(waterfall_plots)) {
#   ggsave(
#     file.path(output_dir, paste0("waterfall_", db_name, ".png")),
#     waterfall_plots[[db_name]],
#     width = 12,
#     height = 6,
#     dpi = 300,
#     bg = "white"
#   )
# }

# Select only the databases we want to show in the combined plot
selected_waterfall_plots <- list(
  peptideatlas = waterfall_plots$peptideatlas,
  hpa_pea = waterfall_plots$hpa_pea,
  hpa_imm = waterfall_plots$hpa_imm,
  gpmdb = waterfall_plots$gpmdb,
  paxdb = waterfall_plots$paxdb,
  hpa_ms = waterfall_plots$hpa_ms
)

# Combine selected waterfall plots into one figure
combined_waterfall <- ggpubr::ggarrange(
  plotlist = selected_waterfall_plots,
  ncol = 2,
  nrow = 3,
  labels = "AUTO",
  font.label = list(size = 12, face = "bold"),
  common.legend = TRUE,
  legend = "bottom"
)

# Save combined waterfall plot
ggsave(
  file.path(output_dir, "02_combined_waterfall.png"),
  combined_waterfall,
  width = 20,
  height = 16,
  dpi = 300,
  bg = "white"
)

# Create final combined plot
message("[Biomarker Analysis] Creating final combined visualization...")

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

# Create a named list for UpSet plot
gene_lists <- list(
  PeptideAtlas = peptideatlas_biomarkers,
  "HPA MS" = hpa_ms_biomarkers,
  "HPA PEA" = hpa_pea_biomarkers,
  "HPA Immunoassay" = hpa_imm_biomarkers,
  GPMDB = gpmdb_biomarkers,
  PAXDB = paxdb_biomarkers
)

# Create UpSet plot
upset_plot <- upset(
  fromList(gene_lists),
  nsets = 6,
  sets = names(gene_lists),
  order.by = "freq",
  nintersects = 30,
  point.size = 3,
  line.size = 1,
  mainbar.y.label = "Biomarker Intersections",
  sets.x.label = "Biomarkers per Technology",
  text.scale = 1.2,
  mb.ratio = c(0.6, 0.4),
  keep.order = TRUE,
  main.bar.color = "#0072B2",
  sets.bar.color = "#0072B2"
)

# Save individual UpSet plot
png(file.path(output_dir, "03_biomarker_technology_overlap.png"), 
    width = 12, height = 8, units = "in", res = 300)
print(upset_plot)
dev.off()

# Save UpSet plot to a temporary file to read it back as a grob
temp_upset_file <- tempfile(fileext = ".png")
png(temp_upset_file, width = 10, height = 8, units = "in", res = 300, bg = "white")  # Reduced width for UpSet
print(upset_plot)
dev.off()

# Read the UpSet plot as a grob
upset_grob <- grid::rasterGrob(png::readPNG(temp_upset_file), interpolate = TRUE)

# Create boxplots for each database
create_boxplot <- function(data, database_name, intensity_col, gene_col = "gene") {
  plot_data <- data %>%
    mutate(
      log_intensity = log10(.data[[intensity_col]] + 1),
      z_score = (log_intensity - mean(log_intensity, na.rm = TRUE)) / sd(log_intensity, na.rm = TRUE),
      gene = .data[[gene_col]],
      is_highlight = gene %in% highlight_biomarkers,
      is_biomarker = gene %in% biomarker_genes
    )
  
  # Calculate statistics for annotation
  stats <- plot_data %>%
    group_by(is_biomarker) %>%
    summarise(
      n = n(),
      mean_z = mean(z_score, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create annotation text
  annotation_text <- sprintf(
    "n(biomarkers) = %d\nn(other) = %d",
    stats$n[stats$is_biomarker],
    stats$n[!stats$is_biomarker]
  )
  
  p <- ggplot(plot_data, aes(x = is_biomarker, y = z_score)) +
    # Add background
    theme_bw() +
    # Add reference lines
    geom_hline(yintercept = c(-2, -1, 0, 1, 2), 
               color = "grey80", linetype = "dashed") +
    # Add boxplot
    geom_boxplot(aes(fill = is_biomarker), width = 0.7, outlier.shape = NA) +
    # Add individual points for biomarkers
    geom_jitter(data = filter(plot_data, is_biomarker),
                width = 0.2, alpha = 0.4, size = 1) +
    # Add highlighted points
    geom_point(data = filter(plot_data, is_highlight),
               color = "#FF4B4B", size = 2) +
    # Add labels for highlighted points
    geom_text_repel(
      data = filter(plot_data, is_highlight),
      aes(label = gene),
      direction = "both",
      nudge_x = 0.3,
      size = 3,
      color = "#FF4B4B",
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.5,
      segment.color = "#FF4B4B",
      segment.alpha = 0.5
    ) +
    # Customize colors
    scale_fill_manual(
      values = c("TRUE" = "#E69F00", "FALSE" = "#69b3a2"),
      name = "Protein Type",
      labels = c("All proteins", "Biomarkers")
    ) +
    # Customize x-axis
    scale_x_discrete(labels = c("FALSE" = "All proteins", "TRUE" = "Biomarkers")) +
    # Add titles and labels
    labs(
      title = database_name,
      x = "",
      y = "Z-Score"
    ) +
    # Add annotation
    annotate(
      "text",
      x = 1.5,
      y = max(plot_data$z_score, na.rm = TRUE),
      label = annotation_text,
      hjust = 0,
      size = 3,
      color = "grey30"
    ) +
    # Customize theme
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

# Create boxplots for each database
boxplots <- list(
  peptideatlas = create_boxplot(peptideatlas, "PeptideAtlas", "norm_PSMs_per_100K"),
  hpa_ms = create_boxplot(hpa_ms, "HPA MS", "Concentration", "Gene"),
  hpa_pea = create_boxplot(hpa_pea, "HPA PEA", "median_npx", "Gene"),
  hpa_imm = create_boxplot(hpa_imm, "HPA Immunoassay", "Concentration", "Gene"),
  gpmdb = create_boxplot(gpmdb, "GPMDB", "total"),
  paxdb = create_boxplot(paxdb, "PAXDB", "abundance")
)

# Combine boxplots into one figure with shared legend
combined_boxplots <- ggpubr::ggarrange(
  boxplots$peptideatlas + theme(legend.position = "none"),
  boxplots$hpa_pea + theme(legend.position = "none"),
  boxplots$hpa_imm + theme(legend.position = "none"),
  boxplots$gpmdb + theme(legend.position = "none"),
  ncol = 2,
  nrow = 2,
  labels = "AUTO",
  font.label = list(size = 12, face = "bold")
) %>%
  ggpubr::annotate_figure(
    top = ggpubr::text_grob("Protein Abundance Distribution by Database",
                    face = "bold", size = 14)
  )

# Add a shared legend at the bottom
legend_plot <- boxplots$peptideatlas + 
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title.position = "top", nrow = 1))
shared_legend <- ggpubr::get_legend(legend_plot)

# Combine plots with shared legend
combined_boxplots_with_legend <- ggpubr::ggarrange(
  combined_boxplots,
  shared_legend,
  ncol = 1,
  heights = c(1, 0.1)
)

# Create biomarker profile matrix visualization
message("[Biomarker Analysis] Creating biomarker profile matrix...")

# Function to get z-scores for a database
get_zscore_data <- function(data, intensity_col, gene_col = "gene") {
  data %>%
    mutate(
      log_intensity = log10(.data[[intensity_col]] + 1),
      z_score = (log_intensity - mean(log_intensity, na.rm = TRUE)) / sd(log_intensity, na.rm = TRUE),
      gene = .data[[gene_col]]
    ) %>%
    select(gene, z_score)
}

# Collect z-scores from all databases
zscore_data <- list(
  PeptideAtlas = get_zscore_data(peptideatlas, "norm_PSMs_per_100K"),
  "HPA MS" = get_zscore_data(hpa_ms, "Concentration", "Gene"),
  "HPA PEA" = get_zscore_data(hpa_pea, "median_npx", "Gene"),
  "HPA Immunoassay" = get_zscore_data(hpa_imm, "Concentration", "Gene"),
  GPMDB = get_zscore_data(gpmdb, "total"),
  PAXDB = get_zscore_data(paxdb, "abundance")
)

# Create a matrix of z-scores for biomarkers
biomarker_matrix <- matrix(NA, 
                          nrow = length(biomarker_genes), 
                          ncol = length(zscore_data),
                          dimnames = list(biomarker_genes, names(zscore_data)))

# Fill the matrix with z-scores
for (db in names(zscore_data)) {
  db_data <- zscore_data[[db]]
  biomarker_matrix[, db] <- db_data$z_score[match(biomarker_genes, db_data$gene)]
}

# Calculate detection frequency for each biomarker
detection_freq <- rowSums(!is.na(biomarker_matrix)) / ncol(biomarker_matrix) * 100

# Calculate mean z-score (when detected)
mean_zscore <- rowMeans(biomarker_matrix, na.rm = TRUE)

# Create data frame for plotting
plot_data <- as.data.frame(biomarker_matrix) %>%
  tibble::rownames_to_column("Biomarker") %>%
  tidyr::pivot_longer(-Biomarker, names_to = "Database", values_to = "Z_score") %>%
  mutate(
    Database = factor(Database, levels = names(zscore_data)),
    Biomarker = factor(Biomarker, levels = biomarker_genes[order(mean_zscore, decreasing = TRUE)]),
    Detection = ifelse(is.na(Z_score), "Not detected", "Detected"),
    Z_score = ifelse(is.na(Z_score), 0, Z_score),  # Replace NA with 0 for visualization
    Z_score_capped = pmin(pmax(Z_score, -3), 3)    # Cap z-scores for better visualization
  )

# Create the profile matrix plot
profile_plot <- ggplot(plot_data, aes(x = Database, y = Biomarker)) +
  # Add tiles colored by z-score
  geom_tile(aes(fill = Z_score_capped), color = "white", size = 0.5) +
  # Add detection indicators
  geom_point(data = filter(plot_data, Detection == "Not detected"),
             color = "grey80", size = 1) +
  # Customize colors
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-3, 3),
    name = "Z-score"
  ) +
  # Add detection frequency
  geom_text(
    data = unique(plot_data[c("Biomarker")]) %>%
      mutate(
        Database = "Detection\nFrequency",
        label = sprintf("%.0f%%", detection_freq[as.character(Biomarker)])
      ),
    aes(label = label),
    hjust = 0.5,
    size = 3
  ) +
  # Customize theme
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  # Add titles
  labs(
    title = "Biomarker Detection Profile Across Databases",
    subtitle = "Z-scores shown by color intensity. Grey dots indicate non-detection."
  )

# Create right panel with UpSet and boxplots
right_panel <- ggpubr::ggarrange(
  upset_grob,
  combined_boxplots_with_legend,
  ncol = 1,
  heights = c(1, 1.2),
  labels = c("(B)", "(C)"),
  font.label = list(size = 16, face = "bold")
)

# Adjust profile matrix plot for left panel
profile_plot_adjusted <- profile_plot +
  theme(
    plot.margin = margin(t = 35, r = 10, b = 10, l = 10),  # Add top margin to align with right panel
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Slightly smaller text
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 12)
  )

# Create final combined plot with side-by-side layout
final_combined_plot <- ggpubr::ggarrange(
  profile_plot_adjusted,
  right_panel,
  ncol = 2,
  widths = c(0.6, 0.4),  # Matrix takes 60% of width, right panel 40%
  labels = c("(A)", ""),
  font.label = list(size = 16, face = "bold")
)

# Save the final combined plot
ggsave(
  file.path(output_dir, "00_final_combined_analysis.png"),
  final_combined_plot,
  width = 20,  # Reduced total width
  height = 20,  # Keep height the same
  dpi = 300,
  bg = "white"
)

# Clean up temporary file
unlink(temp_upset_file)

# Also save individual plots for flexibility
ggsave(
  file.path(output_dir, "01_biomarker_profile_matrix.png"),
  profile_plot,
  width = 15,
  height = 20,
  dpi = 300,
  bg = "white"
)

message("\nAnalysis complete. Files saved to: ", output_dir)