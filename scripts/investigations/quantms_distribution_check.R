#!/usr/bin/env Rscript
# ============================================================================
# INVESTIGATION: QuantMS Abundance Distribution
# ============================================================================
# Description: This script investigates the per-protein abundance distributions
# in the QuantMS dataset to inform the best aggregation strategy.
# ============================================================================

# --- 1. Preamble: Load Packages and Utilities ---

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(stringr)
})

# Load project utilities
source("scripts/data_processing/simple_id_mapping.R")
source("scripts/utilities/plot_themes.R")
source("scripts/config/analysis_config.R")

# --- 2. Data Loading and Preparation ---

message("Loading and preparing QuantMS plasma data...")

# Path to the data
quantms_dir <- file.path(PROJECT_CONFIG$directories$data_raw, "quantms/plasma")

# Read and combine all QuantMS plasma files
quantms_files <- list.files(quantms_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(quantms_files) == 0) {
  stop("No QuantMS plasma data files found.")
}
raw_data <- lapply(quantms_files, read_csv, show_col_types = FALSE) %>%
  bind_rows()

# Map protein accessions to gene symbols
message("Mapping protein IDs to gene symbols...")
raw_data$gene <- convert_to_gene_symbol(raw_data$ProteinName, force_mapping = FALSE)

# Filter for proteins that mapped successfully to a single gene
mapped_data <- raw_data %>%
  filter(sapply(gene, function(x) length(x) == 1 && !is.na(x))) %>%
  mutate(gene = unlist(gene)) %>%
  filter(!is.na(Ibaq) & Ibaq > 0)

message(sprintf("Data prepared: %d measurements for %d unique genes.", nrow(mapped_data), n_distinct(mapped_data$gene)))

# --- 3. Identify Candidate Proteins for Visualization ---

message("Identifying candidate proteins for plotting...")

# Calculate summary statistics for each gene
gene_stats <- mapped_data %>%
  group_by(gene) %>%
  summarise(
    n_samples = n(),
    median_ibaq = median(Ibaq),
    mean_ibaq = mean(Ibaq),
    sd_ibaq = sd(Ibaq),
    variance_ibaq = var(Ibaq),
    .groups = "drop"
  ) %>%
  # Calculate coefficient of variation, a normalized measure of dispersion
  mutate(cv = sd_ibaq / mean_ibaq) %>%
  arrange(desc(variance_ibaq))

# Select interesting candidates:
# 1. Highest variance protein (with at least 10 observations)
high_var_gene <- gene_stats %>% filter(n_samples >= 10) %>% top_n(1, variance_ibaq) %>% pull(gene)

# 2. Low variance protein (with a similar number of observations to the high var one)
n_obs_high_var <- gene_stats %>% filter(gene == high_var_gene) %>% pull(n_samples)
low_var_gene <- gene_stats %>% 
  filter(n_samples >= (n_obs_high_var - 5) & n_samples <= (n_obs_high_var + 5)) %>%
  arrange(variance_ibaq) %>%
  head(1) %>%
  pull(gene)

# 3. A known biomarker for context
biomarker_gene <- "APP"
if (!biomarker_gene %in% gene_stats$gene) {
    biomarker_gene <- gene_stats$gene[10] # Fallback if APP not found
}

# Combine the candidates
candidate_genes <- c(high_var_gene, low_var_gene, biomarker_gene)
message(paste("Selected candidate genes:", paste(candidate_genes, collapse = ", ")))

# --- 4. Create Visualization ---

message("Generating diagnostic plots...")

# Function to create a distribution plot for a single gene
create_dist_plot <- function(gene_symbol) {
  
  plot_data <- mapped_data %>% filter(gene == gene_symbol)
  
  if (nrow(plot_data) == 0) return(NULL)
  
  # Calculate aggregation metrics
  agg_metrics <- plot_data %>%
    summarise(
      Median = median(Ibaq),
      Mean = mean(Ibaq),
      Max = max(Ibaq)
    )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = gene, y = Ibaq)) +
    geom_jitter(width = 0.1, alpha = 0.3, color = "grey40", shape = 16) +
    geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA, fill = "#2E86AB") +
    
    # Add lines for different aggregation methods
    geom_hline(aes(yintercept = agg_metrics$Median, color = "Median"), linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = agg_metrics$Mean, color = "Mean"), linetype = "dashed", size = 1) +
    geom_hline(aes(yintercept = agg_metrics$Max, color = "Max"), linetype = "dashed", size = 1) +
    
    scale_y_log10(labels = scales::label_scientific()) +
    scale_color_manual(
      name = "Aggregation",
      values = c("Median" = "#E15759", "Mean" = "#F28E2B", "Max" = "#4E79A7"),
      breaks = c("Median", "Mean", "Max")
    ) +
    theme_blood_proteomics() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right"
    ) +
    labs(
      title = paste("Distribution for:", gene_symbol),
      subtitle = paste(nrow(plot_data), "observations"),
      x = NULL,
      y = "iBAQ (log10 scale)"
    )
  
  return(p)
}

# Generate plots for all candidates
plot_list <- lapply(candidate_genes, create_dist_plot)

# Combine the plots
combined_plot <- wrap_plots(plot_list, ncol = 3, guides = "collect") +
  plot_annotation(
    title = "QuantMS iBAQ Distribution for Select Proteins",
    subtitle = "Comparing Median, Mean, and Max Aggregation Strategies",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 14, hjust = 0.5))
  )

# --- 5. Save Output ---

output_dir <- "outputs/plots/investigations"
output_file <- file.path(output_dir, "quantms_aggregation_investigation.png")

message(paste("Saving plot to:", output_file))
ggsave(output_file, combined_plot, width = 18, height = 7, dpi = 300, bg = "white")

message("✅ Investigation complete.") 