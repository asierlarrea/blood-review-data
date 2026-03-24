#!/usr/bin/env Rscript
# ============================================================================
# PLASMA PROTEIN ANALYSIS - DIVIDED PLOTS
# ============================================================================
# Description: Creates separate plots for panels A/B/C and panel D
# This allows testing panel D width independently
# ============================================================================

# Load utilities and configuration
source("scripts/utilities/load_packages.R")
source("scripts/config/analysis_config.R")
source("scripts/utilities/data_loader.R")
source("scripts/utilities/plot_themes.R")

# Set up environment
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "patchwork", "ggpubr", "UpSetR", "tibble", "ggupset", "purrr")
load_packages(required_packages)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
force_mapping <- "--force-mapping" %in% args

message(paste(rep("=", 60), collapse = ""))
message("PLASMA PROTEIN ANALYSIS - DIVIDED PLOTS")
message(paste(rep("=", 60), collapse = ""))
message("Script started at: ", Sys.time())

# Wrap main execution in try-catch
tryCatch({
  
# Step 1: Load all data sources
message("\n[STEP 1] Loading data sources...")
all_sources <- get_all_data_sources()
sources_to_load <- setdiff(all_sources, "hpa_immunoassay")

data_list <- load_multiple_sources(source_names = sources_to_load, force_mapping = force_mapping)

# Load QuantMS data with prevalence filter
quantms_data <- load_quantms_data(sample_type = "plasma", force_mapping = force_mapping)

# Add it to the main data list
if (!is.null(quantms_data)) data_list$quantms <- quantms_data

if (length(data_list) == 0) {
  stop("No data sources could be loaded. Please check your data files.")
}

# Step 2: Combine and normalize data
message("\n[STEP 2] Combining and normalizing data...")
combined_data <- combine_data_sources(data_list)

# Apply normalization
message("  - Filtering out zero abundance values...")
normalized_data <- combined_data %>%
  filter(abundance > 0)

message("  - Applying log10 transformation (where needed)...")
normalized_data <- normalized_data %>%
  mutate(
    log_abundance = ifelse(source == "HPA PEA",
                          abundance,              # Keep linear for HPA PEA (samples_above_lod)
                          log10(abundance + 1))   # Apply log to all other sources
  )

message("  - Applying quantile-to-normal normalization...")
normalized_data <- normalized_data %>%
  group_by(source) %>%
  mutate(
    rank_quantile = rank(log_abundance) / (n() + 1),
    z_score = qnorm(rank_quantile)
  ) %>%
  ungroup()

# Set up output directory
plot_dir <- get_output_path("01_plasma_protein_analysis", subdir = "plots")

#' Create UpSet plot for integration into comprehensive panel using ggupset
#'
create_upset_plot_for_panel <- function(gene_lists, set_colors) {
  
  # Convert gene lists to a tidy format for ggupset
  upset_data <- tibble()
  
  for (db_name in names(gene_lists)) {
    if (length(gene_lists[[db_name]]) > 0) {
      db_data <- tibble(
        gene = gene_lists[[db_name]],
        database = db_name
      )
      upset_data <- bind_rows(upset_data, db_data)
    }
  }
  
  # Create the data structure for ggupset
  upset_data_wide <- upset_data %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = database, values_from = present, values_fill = 0) %>%
    rowwise() %>%
    mutate(
      databases = list(names(gene_lists)[c_across(-gene) == 1])
    ) %>%
    ungroup() %>%
    select(gene, databases)
  
  # Count frequencies and filter out those below 100
  intersection_counts <- upset_data_wide %>%
    count(databases) %>%
    filter(n >= 100) %>%
    arrange(desc(n))
  
  # Create UpSet plot using ggupset with single blue color
  panel_plot <- upset_data_wide %>%
    filter(map_chr(databases, paste, collapse = ",") %in% 
           map_chr(intersection_counts$databases, paste, collapse = ",")) %>%
    ggplot(aes(x = databases)) +
    geom_bar(fill = "#4575b4", alpha = 0.85) +
    geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.3, size = 8) +
    scale_x_upset(order_by = "freq") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_text(size = 20, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      strip.text = element_text(size = 22, face = "bold"),
      strip.background = element_blank()
    ) +
    labs(
      title = "(A) Intersection between different data sources",
      x = "Data Source Combinations",
      y = "Number of Proteins"
    )
  
  return(panel_plot)
}

# === CREATE PANELS A, B, C ===
message("\n[STEP 3] Creating panels A, B, C...")

# === SECTION 1: OVERLAP ANALYSIS ===

# Panel A: UpSet Plot for Protein Overlap (Enhanced) - renamed from Panel B
upset_data_panel <- normalized_data %>%
  select(gene, source, technology) %>%
  distinct() %>%
  filter(!is.na(gene), gene != "")

gene_lists_panel <- split(upset_data_panel$gene, upset_data_panel$source)
gene_lists_panel <- lapply(gene_lists_panel, function(x) unique(x[!is.na(x) & x != ""]))

tech_mapping_panel <- upset_data_panel %>%
  select(source, technology) %>%
  distinct() %>%
  deframe()

tech_colors_panel <- get_plot_colors("technology")
set_colors_panel <- tech_colors_panel[tech_mapping_panel[names(gene_lists_panel)]]

panel_A <- create_upset_plot_for_panel(gene_lists_panel, set_colors_panel) +
  theme(
    plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
  ) +
  labs(
    title = "(A) Intersection between different data sources",
    x = "Data Source Combinations",
    y = "Number of Proteins"
  )

# === SECTION 2: CROSS-DATABASE ANALYSIS ===

# Panel B: Cross-Database Dot Plot (Enhanced) - renamed from Panel D
peptideatlas_data <- normalized_data %>%
  filter(source == "PeptideAtlas") %>%
  arrange(z_score) %>%
  mutate(order = row_number()) %>%
  select(gene, order)

dot_plot_data <- normalized_data %>%
  inner_join(peptideatlas_data, by = "gene") %>%
  filter(source %in% c("PeptideAtlas", "HPA MS", "PAXDB", "GPMDB")) %>%
  group_by(gene, source) %>%
  summarise(
    z_score = median(z_score, na.rm = TRUE),
    order = first(order),
    .groups = "drop"
  )

db_colors <- get_plot_colors("databases")
db_colors["PeptideAtlas"] <- "#1a1a1a"  # Darker reference color

panel_B <- ggplot(dot_plot_data, aes(x = order, y = z_score, color = source)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = db_colors, name = "Data Source") +
  scale_x_continuous(labels = scales::comma) +
  theme_blood_proteomics() +
  theme(
    plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
    plot.subtitle = element_text(size = 20),  # Added styling for subtitle
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    panel.grid.major.y = element_line(color = "grey90", size = 0.3)
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 6, alpha = 0.9))) +
  labs(
    title = "(B) Protein abundance/detection frequency correlation with PeptideAtlas",
    x = NULL,  # Removed x-axis title
    y = "Quantile-normalized Values"
  )

# === SECTION 3: QUANTMS SAMPLE DISTRIBUTION ANALYSIS ===

# Panel C: QuantMS Sample Count Distribution - renamed from Panel E
message("  - Creating Panel C: QuantMS sample distribution...")

# Get the QuantMS data with quantile-normalized z-scores
quantms_sample_data <- normalized_data %>%
  filter(source == "quantms") %>%
  select(gene, z_score, abundance, log_abundance)

# Calculate sample counts directly from raw data
message("  - Calculating sample counts from raw QuantMS data...")
quantms_dir <- file.path("data/raw/quantms/plasma")
quantms_files <- list.files(quantms_dir, pattern = "\\.csv$", full.names = TRUE)

# Count samples per gene across all files
sample_counts <- map_dfr(quantms_files, ~{
  read_csv(.x, show_col_types = FALSE) %>%
    select(protein_accession = ProteinName, Ibaq = IbaqNorm, SampleID) %>%
    filter(!is.na(Ibaq), Ibaq > 0)
}) %>%
  # Map to genes (use existing mapping)
  mutate(gene = convert_to_gene_symbol(protein_accession, force_mapping = FALSE)) %>%
  filter(!is.na(gene), gene != "") %>%
  group_by(gene) %>%
  summarise(sample_count = n_distinct(SampleID), .groups = "drop") %>%
  filter(sample_count >= 10) # Only genes present in 10+ samples (reduces sample count bias)

# Merge sample count information with abundance data
quantms_sample_data <- quantms_sample_data %>%
  left_join(sample_counts, by = "gene") %>%
  filter(!is.na(sample_count))

# Get PeptideAtlas data for comparison
peptideatlas_genes <- normalized_data %>%
  filter(source == "PeptideAtlas") %>%
  arrange(desc(log_abundance)) %>%
  mutate(peptideatlas_rank = row_number())

# Define gene categories for highlighting

# 1) Top 10 abundant genes in PeptideAtlas that are also in quantms
top_peptideatlas_in_quantms <- peptideatlas_genes %>%
  filter(gene %in% quantms_sample_data$gene) %>%
  slice_head(n = 10) %>%
  pull(gene)

# 2) Bottom 10 low-abundant genes in PeptideAtlas that are also in quantms  
low_peptideatlas_in_quantms <- peptideatlas_genes %>%
  filter(gene %in% quantms_sample_data$gene) %>%
  slice_tail(n = 10) %>%
  pull(gene)

# 3) Top 10 quantms genes not in PeptideAtlas - prioritize by z-score, then sample count
quantms_only_high <- quantms_sample_data %>%
  filter(!gene %in% peptideatlas_genes$gene) %>%
  arrange(desc(z_score), desc(sample_count)) %>%
  slice_head(n = 10) %>%
  pull(gene)

# 4) Bottom 10 low-abundant quantms genes not in PeptideAtlas
quantms_only_low <- quantms_sample_data %>%
  filter(!gene %in% peptideatlas_genes$gene) %>%
  arrange(z_score) %>%
  slice_head(n = 10) %>%
  pull(gene)

# 5) Well-known plasma proteins from literature
plasma_literature_genes <- c("ALB", "APOA1", "APOA2", "APOB", "APOE", "FGB", "FGG", "FGA", "HP", "TF")
plasma_lit_in_data <- intersect(plasma_literature_genes, quantms_sample_data$gene)

# Create sample count bins and gene categories
quantms_sample_data <- quantms_sample_data %>%
  mutate(
    sample_group = case_when(
      sample_count >= 10 & sample_count <= 20 ~ "10-20",
      sample_count >= 21 & sample_count <= 39 ~ "21-39", 
      sample_count >= 40 ~ "40+",
      TRUE ~ "Other"
    ),
    sample_group = factor(sample_group, levels = c("10-20", "21-39", "40+")),
             gene_category = case_when(
         gene %in% top_peptideatlas_in_quantms ~ "High PeptideAtlas (shared)",
         gene %in% low_peptideatlas_in_quantms ~ "Low PeptideAtlas (shared)", 
         gene %in% quantms_only_high ~ "High quantms (unique)",
         gene %in% quantms_only_low ~ "Low quantms (unique)",
         gene %in% plasma_lit_in_data ~ "Plasma literature",
         TRUE ~ "Other"
       ),
       gene_category = factor(gene_category, levels = c("High PeptideAtlas (shared)", "Low PeptideAtlas (shared)", 
                                                       "High quantms (unique)", "Low quantms (unique)",
                                                       "Plasma literature", "Other"))
          ) %>%
  filter(!is.na(sample_group), sample_group != "Other")

# Print summary of highlighted genes
message(sprintf("  - High PeptideAtlas (shared): %s", paste(top_peptideatlas_in_quantms, collapse = ", ")))
message(sprintf("  - Low PeptideAtlas (shared): %s", paste(low_peptideatlas_in_quantms, collapse = ", ")))
message(sprintf("  - High quantms (unique): %s", paste(quantms_only_high, collapse = ", ")))
message(sprintf("  - Low quantms (unique): %s", paste(quantms_only_low, collapse = ", ")))
message(sprintf("  - Plasma literature: %s", paste(plasma_lit_in_data, collapse = ", ")))

# Define colors and shapes for gene categories
category_colors <- c(
  "High PeptideAtlas (shared)" = "#d73027",      # Red
  "Low PeptideAtlas (shared)" = "#4575b4",       # Blue  
  "High quantms (unique)" = "#85218e",           # Orange
  "Low quantms (unique)" = "#1a7922",            # Light Blue
  "Plasma literature" = "#a2aa09",               # Purple
  "Other" = "grey40"                             # Grey
)

category_shapes <- c(
  "High PeptideAtlas (shared)" = 16,             # Circle
  "Low PeptideAtlas (shared)" = 17,              # Triangle
  "High quantms (unique)" = 15,                  # Square
  "Low quantms (unique)" = 18,                   # Diamond
  "Plasma literature" = 8,                       # Star
  "Other" = 20                                   # Small circle
)

# Create violin plot with highlighted gene categories
panel_C <- ggplot(quantms_sample_data, aes(x = sample_group, y = z_score)) +
  geom_violin(alpha = 0.6, scale = "width", trim = TRUE, width = 0.7, fill = "lightgrey", color = "grey60") +
  # Background points for all other genes
  geom_jitter(data = filter(quantms_sample_data, gene_category == "Other"),
              alpha = 0.3, size = 2.5, width = 0.15, color = "grey50") +
  # Highlighted points for special gene categories
  geom_jitter(data = filter(quantms_sample_data, gene_category != "Other"),
              aes(color = gene_category, shape = gene_category), 
              alpha = 0.9, size = 4.0, width = 0.15, seed = 42,
              show.legend = c(color = FALSE, shape = TRUE)) +
  # Add gene labels with dashed lines for highlighted genes
  ggrepel::geom_label_repel(
    data = filter(quantms_sample_data, gene_category != "Other"),
    aes(label = gene, color = gene_category),
    size = 5,
    alpha = 0.9,
    fontface = "bold",
    fill = "white",
    label.padding = 0.3,
    label.r = 0.15,
    box.padding = 0.4,
    point.padding = 0.6,
    segment.linetype = "dashed",
    segment.size = 0.7,
    segment.alpha = 0.8,
    max.overlaps = Inf,
    force = 2,
    seed = 42
  ) +
  geom_boxplot(width = 0.1, alpha = 0.9, outlier.size = 0, outlier.alpha = 0, 
               show.legend = FALSE, color = "black", fill = "white") +
  stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red", alpha = 0.8) +
  scale_color_manual(values = category_colors, guide = "none") +
  scale_shape_manual(values = category_shapes, name = "Gene Category") +
  theme_blood_proteomics() +
  theme(
    plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 22, face = "plain"),
    axis.text.x = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size = 22, face = "plain"),
    legend.position = "bottom",
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 22),
    legend.key.size = unit(2.0, "lines"),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    panel.grid.major.y = element_line(color = "grey90", size = 0.3)
  ) +
  guides(
    shape = guide_legend(
      nrow = 2,
      override.aes = list(
        size = 6, 
        alpha = 1, 
        color = c(
          "High PeptideAtlas (shared)" = "#d73027",
          "Low PeptideAtlas (shared)" = "#4575b4", 
          "High quantms (unique)" = "#85218e",
          "Low quantms (unique)" = "#1a7922",
          "Plasma literature" = "#a2aa09"
        )
      )
    )
  ) +
  labs(
    title = "(C) QuantMS protein abundance by sample presence",
    x = "Number of Samples", 
    y = "Quantile-normalized Abundance"
  )

# Combine panels A, B, C
message("  - Combining panels A, B, C...")
panels_ABC <- panel_A / panel_B / panel_C +
  plot_layout(heights = c(1, 1, 1.2)) &
  theme(legend.position = "right")

# Verify panels were created
message(sprintf("  - Panel A class: %s", class(panel_A)[1]))
message(sprintf("  - Panel B class: %s", class(panel_B)[1]))
message(sprintf("  - Panel C class: %s", class(panel_C)[1]))
message(sprintf("  - Combined panels class: %s", class(panels_ABC)[1]))

# Save panels A, B, C - disable RStudio graphics device to force file output
message("\n[STEP 4] Saving panels A, B, C...")
abc_png_path <- file.path(plot_dir, "00_comprehensive_plasma_analysis_panels_ABC.png")
abc_tiff_path <- file.path(plot_dir, "00_comprehensive_plasma_analysis_panels_ABC.tiff")

# Close any RStudio graphics device
if (exists(".rs.getActiveDevice")) {
  tryCatch({
    while (dev.cur() != 1) dev.off()
  }, error = function(e) {})
}

# Save PNG - use png() device directly and ensure plot is rendered
tryCatch({
  # Open PNG device
  png(abc_png_path, width = 12, height = 20, units = "in", res = 450, bg = "white", type = "cairo")
  # Render the plot - try both print and plot for patchwork
  if (inherits(panels_ABC, "patchwork")) {
    print(panels_ABC)
  } else {
    plot(panels_ABC)
  }
  # Close device
  dev.off()
  # Verify file was created and has content
  if (file.exists(abc_png_path) && file.info(abc_png_path)$size > 0) {
    message(sprintf("  ✓ PNG saved: %s (size: %d bytes)", abc_png_path, file.info(abc_png_path)$size))
  } else {
    message(sprintf("  ✗ PNG file is empty or not created: %s", abc_png_path))
  }
}, error = function(e) {
  message(sprintf("  ✗ Error saving PNG: %s", e$message))
  if (dev.cur() != 1) dev.off()
})

# Save TIFF
tryCatch({
  tiff(abc_tiff_path, width = 12, height = 20, units = "in", res = 450, compression = "lzw", bg = "white", type = "cairo")
  # Render the plot - try both print and plot for patchwork
  if (inherits(panels_ABC, "patchwork")) {
    print(panels_ABC)
  } else {
    plot(panels_ABC)
  }
  dev.off()
  # Verify file was created and has content
  if (file.exists(abc_tiff_path) && file.info(abc_tiff_path)$size > 0) {
    message(sprintf("  ✓ TIFF saved: %s (size: %d bytes)", abc_tiff_path, file.info(abc_tiff_path)$size))
  } else {
    message(sprintf("  ✗ TIFF file is empty or not created: %s", abc_tiff_path))
  }
}, error = function(e) {
  message(sprintf("  ✗ Error saving TIFF: %s", e$message))
  if (dev.cur() != 1) dev.off()
})

# === CREATE PANEL D ===
message("\n[STEP 5] Creating Panel D...")

# Load biomarker gene list
biomarker_file <- "data/metadata/biomarkers_list.csv"
biomarkers <- read_csv(biomarker_file, show_col_types = FALSE)
biomarker_genes <- unique(biomarkers$gene_name)

# Get unique sources
unique_sources <- unique(normalized_data$source)

# Collect z-scores from all databases
zscore_data <- list()
for (source_name in unique_sources) {
  source_data <- normalized_data %>% filter(source == source_name)
  if (nrow(source_data) > 0) {
    zscore_data[[source_name]] <- source_data %>%
      select(gene, z_score)
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

# Create data frame for plotting - include "Detection Frequency" as a Database level
database_levels <- c(names(zscore_data), "Detection\nFrequency")

plot_data <- as.data.frame(biomarker_matrix) %>%
  tibble::rownames_to_column("Biomarker") %>%
  tidyr::pivot_longer(-Biomarker, names_to = "Database", values_to = "Z_score") %>%
  mutate(
    Database = factor(Database, levels = database_levels),
    Biomarker = factor(Biomarker, levels = biomarker_genes[order(mean_zscore, decreasing = TRUE)]),
    Detection = ifelse(is.na(Z_score), "Not detected", "Detected"),
    Z_score = ifelse(is.na(Z_score), 0, Z_score),
    Z_score_capped = pmin(pmax(Z_score, -3), 3)
  )

# Create Panel D
panel_D <- ggplot(plot_data, aes(x = Database, y = Biomarker)) +
  geom_tile(data = filter(plot_data, Database != "Detection\nFrequency"),
            aes(fill = Z_score_capped), color = "white", linewidth = 0.25, width = 0.65, height = 0.85, alpha = 1.0) +
  geom_point(data = filter(plot_data, Detection == "Not detected", Database != "Detection\nFrequency"),
             color = "grey80", size = 0.75) +
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "#E8E8E8",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-3, 3),
    name = "Z-score",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  geom_text(
    data = unique(plot_data[c("Biomarker")]) %>%
      mutate(
        Database = "Detection\nFrequency",
        label = sprintf("%.0f%%", detection_freq[as.character(Biomarker)])
      ),
    aes(label = label),
    hjust = 0.5,
    size = 7
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(1.5, "cm"),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(20, 10, 20, 20)
  ) +
  labs(
    title = "(D) Biomarker Detection Profile Across Databases"
  )

# Save Panel D with 12 inches width only
message("\n[STEP 6] Saving Panel D (12 inches width)...")
message(sprintf("  - Panel D class: %s", class(panel_D)[1]))

filename <- "00_comprehensive_plasma_analysis_panel_D_width_12.png"
filepath <- file.path(plot_dir, filename)
tryCatch({
  ggsave(filepath, plot = panel_D, width = 12, height = 20, dpi = 450, device = "png", bg = "white")
  message(sprintf("  ✓ Saved Panel D with width = 12 inches: %s", filepath))
}, error = function(e) {
  message(sprintf("  ✗ Error saving Panel D: %s", e$message))
})

message("\n✓ Script execution completed!")
message(sprintf("  - Panels A, B, C: %s", file.path(plot_dir, "00_comprehensive_plasma_analysis_panels_ABC.png")))
message(sprintf("  - Panel D (multiple widths): %s", plot_dir))
message("Script finished at: ", Sys.time())

}, error = function(e) {
  message("\n✗ ERROR OCCURRED:")
  message(sprintf("  Error message: %s", e$message))
  message(sprintf("  Error at: %s", deparse(e$call)))
  traceback()
  quit(status = 1)
})

