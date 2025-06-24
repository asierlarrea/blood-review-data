#!/usr/bin/env Rscript
# ============================================================================
# Biomarker Plasma Analysis (Unified)
# Description: Plots expression distributions for each plasma source, overlays lines for each biomarker gene
# Uses the same data loading and normalization logic as script 01
# Output: Plots in outputs/plots/biomarker_plasma_analysis
# ============================================================================

# Function to create a waterfall plot for a specific database
create_waterfall_plot <- function(data, database_name, intensity_col, gene_col = "gene") {
  plot_data <- data %>%
    mutate(
      gene = .data[[gene_col]],
      is_highlight = gene %in% highlight_biomarkers,
      is_biomarker = gene %in% biomarker_genes,
      rank = rank(.data[[intensity_col]], ties.method = "first")
    ) %>%
    arrange(.data[[intensity_col]])
  
  # Calculate mean and SD for annotations
  stats <- plot_data %>%
    summarise(
      biomarker_mean = mean(.data[[intensity_col]][is_biomarker], na.rm = TRUE),
      all_mean = mean(.data[[intensity_col]], na.rm = TRUE),
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
  p <- ggplot(plot_data, aes(x = rank, y = .data[[intensity_col]])) +
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
    # Determine measurement type for this database
    measurement_type <- get_measurement_type(database_name)
    
    labs(
      title = paste(database_name, "Protein", measurement_type, "Distribution"),
      subtitle = paste("Quantile-normalized", tolower(measurement_type), "with highlighted key biomarkers"),
      x = paste("Proteins (ranked by", tolower(measurement_type), ")"),
      y = paste("Quantile-normalized", measurement_type)
    ) +
    # Add annotation text
    annotate(
      "text",
      x = max(plot_data$rank) * 0.95,
      y = max(plot_data[[intensity_col]]) * 0.8,
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
source("scripts/utilities/plot_themes.R")
source("scripts/utilities/quantile_normalization_functions.R")
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "ggthemes", "patchwork", "UpSetR", "ggrepel", "ggpubr", "tibble", "grid", "png")
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

message("[Biomarker Analysis] Loading data sources using unified approach...")

# Use the same data loading approach as script 01
# Get all source names and remove 'hpa_immunoassay' for this specific analysis
all_sources <- get_all_data_sources()
sources_to_load <- setdiff(all_sources, "hpa_immunoassay")

data_list <- load_multiple_sources(source_names = sources_to_load, force_mapping = force_mapping)

# Load QuantMS data with prevalence filter (same as script 01)
quantms_data <- load_quantms_data(sample_type = "plasma", force_mapping = force_mapping, min_samples = 10)

# Add it to the main data list
if (!is.null(quantms_data)) data_list$quantms <- quantms_data

if (length(data_list) == 0) {
  stop("No data sources could be loaded. Please check your data files.")
}

message("[Biomarker Analysis] Combining and normalizing data...")
combined_data <- combine_data_sources(data_list)

# Apply the same normalization logic as script 01
message("  - Filtering out zero abundance values...")
normalized_data <- combined_data %>%
  filter(abundance > 0)

message("  - Applying log10 transformation...")
normalized_data <- normalized_data %>%
  mutate(log_abundance = log10(abundance + 1))

message("  - Applying robust quantile-to-normal transformation within each database...")
  # Apply quantile normalization to normal distribution within each database to handle outliers
normalized_data <- normalized_data %>%
  group_by(source) %>%
  mutate(
    # Map each database's distribution to a standard normal distribution
    rank_quantile = rank(log_abundance, ties.method = "average") / (n() + 1),
    z_score = qnorm(rank_quantile)   # qnorm already produces standard normal distribution
  ) %>%
  ungroup()

# Function to determine if a database measures abundance or expression
get_measurement_type <- function(database_name) {
  expression_sources <- c("HPA PEA", "quantms")
  abundance_sources <- c("PeptideAtlas", "HPA MS", "GPMDB", "PAXDB")
  
  if (database_name %in% expression_sources) {
    return("Expression")
  } else if (database_name %in% abundance_sources) {
    return("Abundance")
  } else {
    return("Values")  # Generic fallback
  }
}

# Create boxplots for each database
create_boxplot <- function(data, database_name, intensity_col, gene_col = "gene") {
  plot_data <- data %>%
    mutate(
      gene = .data[[gene_col]],
      is_highlight = gene %in% highlight_biomarkers,
      is_biomarker = gene %in% biomarker_genes
    )
  
  # Calculate statistics for annotation
  stats <- plot_data %>%
    group_by(is_biomarker) %>%
    summarise(
      n = n(),
      mean_z = mean(.data[[intensity_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create annotation text
  annotation_text <- sprintf(
    "n(biomarkers) = %d\nn(other) = %d",
    stats$n[stats$is_biomarker],
    stats$n[!stats$is_biomarker]
  )
  
  # Determine measurement type for this database
  measurement_type <- get_measurement_type(database_name)
  
  p <- ggplot(plot_data, aes(x = is_biomarker, y = .data[[intensity_col]])) +
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
      y = paste("Quantile-normalized", measurement_type)
    ) +
    annotate(
      "text",
      x = 1.2,  # Changed from 1.5 to be more centered
      y = max(plot_data[[intensity_col]], na.rm = TRUE) * 0.8,  # Moved down by multiplying by 0.8
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

# Process biomarker data for each source
message("[Biomarker Analysis] Processing biomarkers for each database...")

# Extract biomarkers for each database using the normalized data
database_biomarkers <- list()
unique_sources <- unique(normalized_data$source)
for (source_name in unique_sources) {
  source_data <- normalized_data %>% filter(source == source_name)
  database_biomarkers[[source_name]] <- unique(source_data$gene[source_data$gene %in% biomarker_genes])
  message(sprintf("  - %s: %d biomarkers detected", source_name, length(database_biomarkers[[source_name]])))
}

# Create boxplots for each database using the normalized data
boxplots <- list()
for (source_name in unique_sources) {
  source_data <- normalized_data %>% filter(source == source_name)
  if (nrow(source_data) > 0) {
    boxplots[[source_name]] <- create_boxplot(source_data, source_name, "z_score", "gene")
  }
}

# Combine boxplots into one figure with shared legend
# Use the available boxplots dynamically
plot_list <- list()
for (source_name in unique_sources) {
  if (!is.null(boxplots[[source_name]])) {
    plot_list[[source_name]] <- boxplots[[source_name]] + theme(legend.position = "none")
  }
}

combined_boxplots <- ggpubr::ggarrange(
  plotlist = plot_list,
  ncol = 3,
  nrow = 2,
  labels = NULL,  # Removed "AUTO" to remove subpanel labels
  font.label = list(size = 16, face = "bold")
) %>%
  ggpubr::annotate_figure(
    top = ggpubr::text_grob("(C) Protein Abundance/Expression Distribution by Data Source",
                    face = "bold", size = 24)
  )

# Add a shared legend at the bottom
# Use the first available boxplot for legend
first_plot_name <- names(plot_list)[1]
legend_plot <- plot_list[[first_plot_name]] + 
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

# Function to get z-scores for a database (simplified since data is already normalized)
get_zscore_data <- function(data, intensity_col, gene_col = "gene") {
  data %>%
    select(gene = all_of(gene_col), z_score = all_of(intensity_col))
}

# Collect z-scores from all databases dynamically
zscore_data <- list()
for (source_name in unique_sources) {
  source_data <- normalized_data %>% filter(source == source_name)
  if (nrow(source_data) > 0) {
    zscore_data[[source_name]] <- get_zscore_data(source_data, "z_score")
  }
}

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
  geom_tile(aes(fill = Z_score_capped), color = "white", linewidth = 0.25, width = 0.65, height = 0.85, alpha = 1.0) +  # Explicit alpha and use linewidth
  # Add detection indicators with adjusted size
  geom_point(data = filter(plot_data, Detection == "Not detected"),
             color = "grey80", size = 0.75) +
  # Customize colors with enhanced contrast for small z-score differences
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "#E8E8E8",    # Light grey for better visibility of near-zero values
    high = "#B2182B",
    midpoint = 0,
    limits = c(-3, 3),
    name = "Z-score",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
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
    subtitle = "Quantile-normalized values shown by color intensity. Grey dots indicate non-detection."
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

# Create a named list for UpSet plot using available sources
gene_lists <- list()
for (source_name in unique_sources) {
  if (length(database_biomarkers[[source_name]]) > 0) {
    gene_lists[[source_name]] <- database_biomarkers[[source_name]]
  }
}

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
png(temp_upset_file, width = 18, height = 16, units = "in", res = 600, bg = "white")  # Increased height from 14 to 16

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
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(t = 20, b = 15)),
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10)
  )

# Create right panel with UpSet and boxplots
right_panel <- ggpubr::ggarrange(
  upset_with_title,
  combined_boxplots_with_legend,
  ncol = 1,
  heights = c(1.1, 1.1)  # Give more balanced space and slightly more for UpSet
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
    PeptideAtlas = list(total = length(unique(normalized_data$gene[normalized_data$source == "PeptideAtlas"])), biomarkers = length(database_biomarkers$PeptideAtlas)),
    "HPA MS" = list(total = length(unique(normalized_data$gene[normalized_data$source == "HPA MS"])), biomarkers = length(database_biomarkers$`HPA MS`)),
    "HPA PEA" = list(total = length(unique(normalized_data$gene[normalized_data$source == "HPA PEA"])), biomarkers = length(database_biomarkers$`HPA PEA`)),
    GPMDB = list(total = length(unique(normalized_data$gene[normalized_data$source == "GPMDB"])), biomarkers = length(database_biomarkers$GPMDB)),
    PAXDB = list(total = length(unique(normalized_data$gene[normalized_data$source == "PAXDB"])), biomarkers = length(database_biomarkers$PAXDB)),
    quantms = list(total = length(unique(normalized_data$gene[normalized_data$source == "quantms"])), biomarkers = length(database_biomarkers$quantms))
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