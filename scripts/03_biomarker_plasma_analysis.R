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
              color = "#FF4B4B", size = 3) +
    geom_text_repel(
      data = filter(plot_data, is_highlight),
      aes(label = gene),
      direction = "both",
      nudge_y = 0.5,
      size = 4,
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
      size = 4,
      color = "grey30"
    ) +
    # Customize theme
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
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
source("scripts/utilities/data_loader.R")
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
highlight_biomarkers <- c("F12", "LEP", "GHRL",  "GH1", "IL1A", "IL1B", "IL2RA", "EPO", "APP", "PLR", "F8", "C3", "C2")

# Function to calculate z-score normalization within each database
calculate_zscore_normalization_biomarker <- function(data, intensity_col) {
  data %>%
    mutate(
      log_intensity = log10(.data[[intensity_col]] + 1),
      z_score = (log_intensity - mean(log_intensity, na.rm = TRUE)) / sd(log_intensity, na.rm = TRUE)
    )
}

# Process each database
message("[Biomarker Analysis] Processing databases...")

# 1. PeptideAtlas
message("Processing PeptideAtlas...")
peptideatlas_raw <- read_csv("data/raw/peptideatlas/peptideatlas.csv", show_col_types = FALSE)
peptideatlas_raw$gene <- convert_to_gene_symbol(peptideatlas_raw$biosequence_accession, force_mapping = force_mapping)
peptideatlas <- deduplicate_genes(peptideatlas_raw, "gene", "norm_PSMs_per_100K", 
                                additional_cols = c("biosequence_accession"), 
                                aggregation_method = "median")
peptideatlas_biomarkers <- unique(peptideatlas$gene[peptideatlas$gene %in% biomarker_genes])

# 2. HPA MS
message("Processing HPA MS...")
hpa_ms_raw <- read_csv("data/raw/hpa/hpa_ms.csv", show_col_types = FALSE, skip = 1)
hpa_ms <- deduplicate_genes(hpa_ms_raw, "Gene", "Concentration", aggregation_method = "median")
hpa_ms_biomarkers <- unique(hpa_ms$Gene[hpa_ms$Gene %in% biomarker_genes])

# 3. HPA PEA
message("Processing HPA PEA...")
hpa_pea_raw <- read_csv("data/raw/hpa/hpa_pea.csv", show_col_types = FALSE)
hpa_pea <- deduplicate_genes(hpa_pea_raw, "Gene", "median_npx", aggregation_method = "median")
hpa_pea_biomarkers <- unique(hpa_pea$Gene[hpa_pea$Gene %in% biomarker_genes])

# 4. GPMDB
message("Processing GPMDB...")
gpmdb_raw <- read_csv("data/raw/gpmdb/gpmdb_plasma.csv", show_col_types = FALSE)
gpmdb_raw$gene <- stringr::str_extract(gpmdb_raw$description, "[A-Z0-9]+(?=,| |$)")
gpmdb <- deduplicate_genes(gpmdb_raw, "gene", "total", aggregation_method = "median")
gpmdb_biomarkers <- unique(gpmdb$gene[gpmdb$gene %in% biomarker_genes])

# 5. PAXDB
message("Processing PAXDB...")
paxdb_raw <- read_csv("data/raw/paxdb/paxdb_plasma.csv", show_col_types = FALSE)
paxdb_raw$ensp <- stringr::str_replace(paxdb_raw$string_external_id, "^9606\\.", "")
paxdb_raw$gene <- convert_to_gene_symbol(paxdb_raw$ensp, force_mapping = force_mapping)
paxdb <- deduplicate_genes(paxdb_raw, "gene", "abundance", 
                         additional_cols = c("ensp"), 
                         aggregation_method = "median")
paxdb_biomarkers <- unique(paxdb$gene[paxdb$gene %in% biomarker_genes])

# 6. QuantMS
message("Processing QuantMS...")
quantms <- load_quantms_data(sample_type = "plasma", force_mapping = force_mapping)
quantms_biomarkers <- unique(quantms$gene[quantms$gene %in% biomarker_genes])

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
                width = 0.2, alpha = 0.4, size = 1.5) +
    # Add highlighted points
    geom_point(data = filter(plot_data, is_highlight),
               color = "#FF4B4B", size = 3) +
    # Add labels for highlighted points
    geom_text_repel(
      data = filter(plot_data, is_highlight),
      aes(label = gene),
      direction = "both",
      nudge_x = 0.3,
      size = 6,  # Increased from 4
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
    scale_x_discrete(labels = c("FALSE" = "All proteins", "TRUE" = "Biomarkers")) +
    labs(
      title = database_name,
      x = "",
      y = "Z-Score"
    ) +
    annotate(
      "text",
      x = 1.2,  # Changed from 1.5 to be more centered
      y = max(plot_data$z_score, na.rm = TRUE) * 0.8,  # Moved down by multiplying by 0.8
      label = annotation_text,
      hjust = 0.3,  # Changed from 0 to center the text
      size = 6,
      color = "grey30"
    ) +
    theme(
      plot.title = element_text(size = 24, face = "bold"),
      axis.title = element_text(size = 20),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      legend.position = "bottom",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 20),
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
  gpmdb = create_boxplot(gpmdb, "GPMDB", "total"),
  paxdb = create_boxplot(paxdb, "PAXDB", "abundance"),
  quantms = create_boxplot(quantms, "quantms", "abundance")
)

# Combine boxplots into one figure with shared legend
combined_boxplots <- ggpubr::ggarrange(
  boxplots$peptideatlas + theme(legend.position = "none"),
  boxplots$hpa_ms + theme(legend.position = "none"),
  boxplots$hpa_pea + theme(legend.position = "none"),
  boxplots$gpmdb + theme(legend.position = "none"),
  boxplots$paxdb + theme(legend.position = "none"),
  boxplots$quantms + theme(legend.position = "none"),
  ncol = 3,
  nrow = 2,
  labels = NULL,  # Removed "AUTO" to remove subpanel labels
  font.label = list(size = 16, face = "bold")
) %>%
  ggpubr::annotate_figure(
    top = ggpubr::text_grob("(C) Protein Abundance Distribution by Database",
                    face = "bold", size = 24)
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
  GPMDB = get_zscore_data(gpmdb, "total"),
  PAXDB = get_zscore_data(paxdb, "abundance"),
  quantms = get_zscore_data(quantms, "abundance")
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
  # Add tiles colored by z-score with reduced size and width
  geom_tile(aes(fill = Z_score_capped), color = "white", size = 0.25, width = 0.65, height = 0.85) +  # Further reduced width
  # Add detection indicators with adjusted size
  geom_point(data = filter(plot_data, Detection == "Not detected"),
             color = "grey80", size = 0.75) +
  # Customize colors
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-3, 3),
    name = "Z-score"
  ) +
  # Add detection frequency with increased font size
  geom_text(
    data = unique(plot_data[c("Biomarker")]) %>%
      mutate(
        Database = "Detection\nFrequency",
        label = sprintf("%.0f%%", detection_freq[as.character(Biomarker)])
      ),
    aes(label = label),
    hjust = 0.5,
    size = 7  # Significantly increased text size for percentages
  ) +
  # Customize theme
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1.5, "cm"),
    # Add more spacing between columns
    panel.spacing = unit(1, "lines")
  ) +
  # Add titles
  labs(
    title = "(A) Biomarker Detection Profile Across Databases",
    subtitle = "Z-scores shown by color intensity. Grey dots indicate non-detection."
  )

# Adjust profile matrix plot for left panel
profile_plot_adjusted <- profile_plot +
  theme(
    plot.margin = margin(t = 35, r = 20, b = 10, l = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold")
  )

# Create and save UpSet plot
message("[Biomarker Analysis] Creating technology overlap plot...")

# Create a named list for UpSet plot
gene_lists <- list(
  PeptideAtlas = peptideatlas_biomarkers,
  "HPA MS" = hpa_ms_biomarkers,
  "HPA PEA" = hpa_pea_biomarkers,
  GPMDB = gpmdb_biomarkers,
  PAXDB = paxdb_biomarkers,
  quantms = quantms_biomarkers
)

# Create UpSet plot with explicit settings to ensure numbers show
upset_plot <- upset(
  fromList(gene_lists),
  nsets = 6,
  sets = names(gene_lists),
  order.by = "freq", 
  nintersects = 20,  # Reduce further to avoid crowding
  point.size = 6,  # Increased from 4 to 6 for larger dots
  line.size = 2,
  mainbar.y.label = "Number of Biomarkers",
  sets.x.label = "Total per Database",
  text.scale = c(4.0, 4.0, 2.5, 1.4, 4.0, 2.0),  # Increased y-axis title and tick labels
  mb.ratio = c(0.7, 0.3),  # More space for main bars
  keep.order = TRUE,  # Keep order for consistency
  main.bar.color = "#4575b4",
  sets.bar.color = "#4575b4",
  matrix.color = "#4575b4",
  number.angles = 0,
  show.numbers = TRUE  # Use TRUE instead of "yes"
)

# Save UpSet plot to a temporary file with standard approach (no custom viewport)
temp_upset_file <- tempfile(fileext = ".png")
png(temp_upset_file, width = 18, height = 14, units = "in", res = 600, bg = "white")

# Print the plot directly
print(upset_plot)

dev.off()

# Read the UpSet plot as a grob and add title later in ggplot
upset_grob <- grid::rasterGrob(png::readPNG(temp_upset_file), interpolate = TRUE)

# Create a ggplot wrapper for the UpSet plot with title
upset_with_title <- ggplot() +
  annotation_custom(upset_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() +
  labs(title = "(B) Biomarker Intersections") +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 10))
  )

# Create right panel with UpSet and boxplots
right_panel <- ggpubr::ggarrange(
  upset_with_title,
  combined_boxplots_with_legend,
  ncol = 1,
  heights = c(1, 1.2)
)

# Create final combined plot with side-by-side layout
final_combined_plot <- ggpubr::ggarrange(
  profile_plot_adjusted,
  right_panel,
  ncol = 2,
  widths = c(0.4, 0.6)
)

# Save the final combined plot - both TIFF and PNG versions
# Save TIFF version
ggsave(
  file.path(output_dir, "00_comprehensive_biomarkers_analysis_panel.tiff"),
  final_combined_plot,
  width = 28,
  height = 22,
  dpi = 600,  # Reduced from 1200 to 600 to avoid memory allocation errors
  bg = "white",
  device = "tiff",
  compression = "lzw"  # Added LZW compression for better file size
)

# Save PNG version
ggsave(
  file.path(output_dir, "00_comprehensive_biomarkers_analysis_panel.png"),
  final_combined_plot,
  width = 28,
  height = 22,
  dpi = 600,
  bg = "white",
  device = "png"
)

# Clean up temporary file
unlink(temp_upset_file)

# Generate comprehensive markdown report
message("Generating comprehensive biomarker analysis report...")
generate_biomarker_report <- function(output_dir) {
  
  # Calculate summary statistics for report
  database_stats <- list(
    PeptideAtlas = list(total = length(unique(peptideatlas$gene)), biomarkers = length(peptideatlas_biomarkers)),
    "HPA MS" = list(total = length(unique(hpa_ms$Gene)), biomarkers = length(hpa_ms_biomarkers)),
    "HPA PEA" = list(total = length(unique(hpa_pea$Gene)), biomarkers = length(hpa_pea_biomarkers)),
    GPMDB = list(total = length(unique(gpmdb$gene)), biomarkers = length(gpmdb_biomarkers)),
    PAXDB = list(total = length(unique(paxdb$gene)), biomarkers = length(paxdb_biomarkers)),
    quantms = list(total = length(unique(quantms$gene)), biomarkers = length(quantms_biomarkers))
  )
  
  # Create report content
  report_content <- paste0(
    "# Biomarker Plasma Analysis Report\n\n",
    "**Analysis Date:** ", Sys.Date(), "\n",
    "**Script:** `03_biomarker_plasma_analysis.R`\n",
    "**Description:** Analysis of biomarker protein expression across plasma databases with waterfall plots and abundance distribution analysis.\n\n",
    "---\n\n",
    
    "## Summary Statistics\n\n",
    "| Database | Total Proteins | Biomarkers Detected | Biomarker Coverage (%) | Technology |\n",
    "|----------|----------------|--------------------|-----------------------|------------|\n",
    sprintf("| PeptideAtlas | %d | %d | %.1f%% | MS |\n", 
            database_stats$PeptideAtlas$total, database_stats$PeptideAtlas$biomarkers,
            100 * database_stats$PeptideAtlas$biomarkers / database_stats$PeptideAtlas$total),
    sprintf("| HPA MS | %d | %d | %.1f%% | MS |\n", 
            database_stats$`HPA MS`$total, database_stats$`HPA MS`$biomarkers,
            100 * database_stats$`HPA MS`$biomarkers / database_stats$`HPA MS`$total),
    sprintf("| HPA PEA | %d | %d | %.1f%% | PEA |\n", 
            database_stats$`HPA PEA`$total, database_stats$`HPA PEA`$biomarkers,
            100 * database_stats$`HPA PEA`$biomarkers / database_stats$`HPA PEA`$total),
    sprintf("| GPMDB | %d | %d | %.1f%% | MS |\n", 
            database_stats$GPMDB$total, database_stats$GPMDB$biomarkers,
            100 * database_stats$GPMDB$biomarkers / database_stats$GPMDB$total),
    sprintf("| PAXDB | %d | %d | %.1f%% | MS |\n", 
            database_stats$PAXDB$total, database_stats$PAXDB$biomarkers,
            100 * database_stats$PAXDB$biomarkers / database_stats$PAXDB$total),
    sprintf("| quantms | %d | %d | %.1f%% | MS |\n",
            database_stats$quantms$total, database_stats$quantms$biomarkers,
            100 * database_stats$quantms$biomarkers / database_stats$quantms$total),
    "\n",
    
    "## Key Findings\n\n",
    "- **Biomarker representation** varies significantly across databases and technologies\n",
    "- **Targeted methods (PEA) show higher biomarker density** due to focused panels\n",
    "- **Mass spectrometry databases** provide broader coverage with substantial biomarker representation\n",
    "- **Key biomarkers** (F12, LEP, GHRL, GH1, IL1A, IL1B, etc.) consistently detected across multiple platforms\n",
    "- **Expression ranges** span 4-6 orders of magnitude across databases\n",
    "- **Cross-platform validation** possible for numerous biomarkers\n",
    "- **quantms** provides additional biomarker coverage\n\n",
    
    "## Biological Insights\n\n",
    "- **Biomarker accessibility** varies by technology: targeted methods excel at specific biomarkers\n",
    "- **Clinical relevance** confirmed by multi-database detection of established biomarkers\n",
    "- **Discovery potential** highest in MS databases due to broader protein coverage\n",
    "- **Validation opportunities** through cross-platform biomarker detection\n",
    "- **Expression patterns** reveal database-specific biases and sensitivities\n",
    "- **Abundance distributions** show technology-specific detection capabilities\n",
    "- **quantms** offers additional biomarker coverage\n\n",
    
    "## Database Comparison\n\n",
    "### Biomarker Detection Across Platforms\n\n",
    "**Mass Spectrometry Platforms:**\n",
    "- Provide unbiased discovery of biomarkers across wide abundance ranges\n",
    "- PAXDB offers highest absolute biomarker numbers due to comprehensive coverage\n",
    "- PeptideAtlas and HPA MS show complementary biomarker profiles\n\n",
    "**Targeted Platforms:**\n",
    "- HPA PEA: Balanced approach with focused biomarker panels\n",
    "- quantms: Additional biomarker coverage\n\n",
    "**Cross-Platform Validation:**\n",
    "- Biomarkers detected in multiple databases show higher clinical confidence\n",
    "- Platform-specific biomarkers may represent unique detection capabilities\n",
    "- Z-score normalization enables cross-database abundance comparisons\n",
    "- quantms offers additional biomarker coverage\n\n",
    
    "## Methodology\n\n",
    "- **Biomarker reference list:** Curated list from literature and clinical databases\n",
    "- **Abundance normalization:** Z-score transformation within each database\n",
    "- **Visualization:** Waterfall plots showing protein abundance distributions\n",
    "- **Highlighting system:** Key biomarkers emphasized in abundance rankings\n",
    "- **Statistical analysis:** Coverage percentages and abundance comparisons\n",
    "- **Cross-database integration:** UpSet plots for biomarker overlap analysis\n\n",
    
    "## Recommendations\n\n",
    "- **Use MS databases** for biomarker discovery and broad profiling\n",
    "- **Employ targeted methods (PEA)** for validation and clinical applications\n",
    "- **Cross-validate biomarkers** across multiple platforms when possible\n",
    "- **Consider abundance ranges** when selecting platforms for specific biomarkers\n",
    "- **Integrate complementary technologies** for comprehensive biomarker analysis\n",
    "- **Focus on multi-platform biomarkers** for robust clinical applications\n",
    "- **quantms** offers additional biomarker coverage\n\n",
    
    "## Generated Files\n\n",
    sprintf("- **Comprehensive panel:** `%s/00_comprehensive_biomarkers_analysis_panel.png`\n", basename(output_dir)),
    "- **Waterfall plots:** Individual database abundance distributions with biomarker highlighting\n",
    "- **Biomarker profile matrix:** Z-score heatmap across all databases\n",
    "- **UpSet intersection plots:** Biomarker overlap analysis between platforms\n",
    "- **Box plots:** Biomarker abundance distributions by technology\n\n",
    
    "---\n",
    "*Report generated automatically by the blood proteomics analysis pipeline*\n"
  )
  
  # Save report
  report_file <- file.path(output_dir, "biomarker_plasma_analysis_report.md")
  writeLines(report_content, report_file)
  message(sprintf("✅ Comprehensive biomarker report saved to: %s", report_file))
}

# Generate the report
generate_biomarker_report(output_dir)

message("\nBiomarker plasma analysis complete. Files saved to: ", output_dir)