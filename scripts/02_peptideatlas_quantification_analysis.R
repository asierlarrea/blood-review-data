#!/usr/bin/env Rscript
# ============================================================================
# PeptideAtlas Quantification Methods Comparison Analysis
# Description: Compare n_observations vs norm_PSMs_per_100K quantification methods
# Output: Plots in outputs/plots/02_peptideatlas_quantification_analysis
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

# Load gene deduplication utility
source("scripts/utilities/gene_deduplication.R")

# Set output directories
output_dir <- get_output_path("peptideatlas_quantification")
plot_dir <- get_output_path("02_peptideatlas_quantification_analysis", subdir = "plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

message("PeptideAtlas Quantification Methods Comparison Analysis")
message("=====================================================")

# Load the complete dataset with new gene mapping
# Try to load processed file first, if not available, create it
processed_file <- "data/processed/peptideatlas_mapped_genes.csv"

if (file.exists(processed_file)) {
  cat("Loading processed dataset...\n")
  peptideatlas_data <- read_csv(processed_file, show_col_types = FALSE)
} else {
  cat("Processed file not found. Loading and mapping original data...\n")
  # Load original data
  peptideatlas_original <- read_csv("data/raw/peptideatlas/peptideatlas.csv", show_col_types = FALSE)
  
  # Map genes
  cat("Mapping genes (this may take a moment)...\n")
  peptideatlas_data <- peptideatlas_original %>%
    mutate(gene = convert_to_gene_symbol(biosequence_accession, force_mapping = force_mapping))
  
  # Create processed directory if it doesn't exist
  if (!dir.exists("data/processed")) {
    dir.create("data/processed", recursive = TRUE)
  }
  
  # Save processed data
  write_csv(peptideatlas_data, processed_file)
  cat("Saved processed data to:", processed_file, "\n")
}

# Convert quantification values to numeric
peptideatlas_data$n_observations <- as.numeric(peptideatlas_data$n_observations)
peptideatlas_data$norm_PSMs_per_100K <- as.numeric(peptideatlas_data$norm_PSMs_per_100K)

# Remove rows with missing quantification data and deduplicate genes
peptideatlas_clean_raw <- peptideatlas_data %>%
  filter(!is.na(n_observations) & !is.na(norm_PSMs_per_100K) & 
         n_observations > 0 & norm_PSMs_per_100K > 0)

# Deduplicate genes for both quantification methods
message("Deduplicating genes for n_observations...")
peptideatlas_n_obs <- deduplicate_genes(peptideatlas_clean_raw, "gene", "n_observations", 
                                       additional_cols = c("biosequence_accession"), 
                                       aggregation_method = "median")

message("Deduplicating genes for norm_PSMs_per_100K...")
peptideatlas_psms <- deduplicate_genes(peptideatlas_clean_raw, "gene", "norm_PSMs_per_100K", 
                                     additional_cols = c("biosequence_accession"), 
                                     aggregation_method = "median")

# Merge the deduplicated data
peptideatlas_clean <- peptideatlas_n_obs %>%
  select(gene, n_observations, biosequence_accession) %>%
  inner_join(
    peptideatlas_psms %>% select(gene, norm_PSMs_per_100K),
    by = "gene"
  )

cat(sprintf("Loaded %d unique genes with quantification data (after deduplication)\n", nrow(peptideatlas_clean)))

# Create correlation plot comparing the two quantification methods
cat("Creating correlation plot...\n")

correlation_coef <- cor(log10(peptideatlas_clean$n_observations), 
                       log10(peptideatlas_clean$norm_PSMs_per_100K))

p1 <- peptideatlas_clean %>%
  ggplot(aes(x = log10(n_observations), y = log10(norm_PSMs_per_100K))) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  labs(
    title = "Correlation between PeptideAtlas Quantification Methods",
    subtitle = sprintf("Pearson correlation: r = %.4f", correlation_coef),
    x = "Log10(Number of Observations)",
    y = "Log10(Normalized PSMs per 100K)",
    caption = "Each point represents a protein in PeptideAtlas"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold")
  )

# Create distribution comparison plot
cat("Creating distribution comparison plot...\n")

# Prepare data for distribution plot
dist_data <- peptideatlas_clean %>%
  select(n_observations, norm_PSMs_per_100K) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = case_when(
      metric == "n_observations" ~ "Number of Observations",
      metric == "norm_PSMs_per_100K" ~ "Normalized PSMs per 100K"
    ),
    log_value = log10(value)
  )

# Create histogram distribution plot
p2a <- dist_data %>%
  ggplot(aes(x = log_value, fill = metric)) +
  geom_histogram(alpha = 0.7, bins = 50) +
  facet_wrap(~metric, ncol = 1,
             labeller = labeller(metric = c("Number of Observations" = "(a) Number of Observations",
                                           "Normalized PSMs per 100K" = "(b) Normalized PSMs per 100K"))) +
  scale_fill_manual(values = c("Number of Observations" = "#E69F00", 
                              "Normalized PSMs per 100K" = "#56B4E9")) +
  labs(
    title = "PeptideAtlas Quantification Methods Distribution",
    subtitle = "Comparison of log-transformed quantification values",
    x = "Log10(Value + offset)",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12)
  )

# Create violin plot with boxplot overlay
p2b <- dist_data %>%
  ggplot(aes(x = metric, y = log_value, fill = metric)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("Number of Observations" = "#E69F00", 
                              "Normalized PSMs per 100K" = "#56B4E9")) +
  labs(
    title = "Abundance Distribution Comparison",
    subtitle = "Distribution of quantification values (log scale)",
    x = "Quantification Method",
    y = "Log10(Value + offset)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Combine distribution plots
p2 <- p2a | p2b + 
  plot_annotation(tag_levels = list(c('(a)', '(b)')))

# Combine all plots into comprehensive view
cat("Creating comprehensive plot...\n")

# Create the comprehensive plot with abundance distribution on left, distribution top right, correlation bottom right
comprehensive_plot <- p2b | (p2a / p1) +
  plot_annotation(
    title = "Comprehensive PeptideAtlas Quantification Methods Analysis",
    subtitle = "Distribution comparison and abundance analysis",
    tag_levels = list(c('(a)', '(b)', '(c)'))
  ) & 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

# Save plots with high resolution PNG format
ggsave(file.path(plot_dir, "correlation_plot.png"), p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(plot_dir, "distribution_plot.png"), p2, width = 15, height = 8, dpi = 300)
ggsave(file.path(plot_dir, "comprehensive_plot.png"), comprehensive_plot, width = 20, height = 12, dpi = 300)

# Print summary statistics
cat("\n=== PeptideAtlas Quantification Analysis Summary ===\n")
cat(sprintf("Total proteins analyzed: %d\n", nrow(peptideatlas_clean)))
cat(sprintf("Correlation between methods: r = %.4f\n", correlation_coef))

message("\nRECOMMENDATION: Use norm_PSMs_per_100K for normalized, quantitative analysis")

# Save data files
message("Saving comparison data...")
write_csv(peptideatlas_clean, file.path(output_dir, "peptideatlas_clean.csv"))

# Create summary report
sink(file.path(output_dir, "quantification_analysis_summary.txt"))
cat("PEPTIDEATLAS QUANTIFICATION METHODS COMPARISON\n")
cat("=============================================\n\n")

cat("DATA SUMMARY:\n")
cat("=============\n")
cat(sprintf("Total proteins with n_observations: %d\n", nrow(peptideatlas_clean)))
cat(sprintf("Total proteins with norm_PSMs_per_100K: %d\n", nrow(peptideatlas_clean)))
cat("\n")

cat("QUANTIFICATION STATISTICS:\n")
cat("==========================\n")
cat("n_observations:\n")
cat(sprintf("  Range: %.6f - %.6f\n", min(peptideatlas_clean$n_observations), max(peptideatlas_clean$n_observations)))
cat(sprintf("  Median: %.6f\n", median(peptideatlas_clean$n_observations)))
cat(sprintf("  Mean: %.6f\n", mean(peptideatlas_clean$n_observations)))
cat(sprintf("  Standard Deviation: %.6f\n", sd(peptideatlas_clean$n_observations)))
cat("\n")

cat("norm_PSMs_per_100K:\n")
cat(sprintf("  Range: %.6f - %.6f\n", min(peptideatlas_clean$norm_PSMs_per_100K), max(peptideatlas_clean$norm_PSMs_per_100K)))
cat(sprintf("  Median: %.6f\n", median(peptideatlas_clean$norm_PSMs_per_100K)))
cat(sprintf("  Mean: %.6f\n", mean(peptideatlas_clean$norm_PSMs_per_100K)))
cat(sprintf("  Standard Deviation: %.6f\n", sd(peptideatlas_clean$norm_PSMs_per_100K)))
cat("\n")

cat("CORRELATION ANALYSIS:\n")
cat("====================\n")
cat(sprintf("Pearson correlation: r = %.4f\n", correlation_coef))
cat("\n")

cat("RECOMMENDATION:\n")
cat("===============\n")
cat("norm_PSMs_per_100K is recommended because:\n")
cat("1. Normalized per 100K PSMs - accounts for run-to-run variation\n")
cat("2. Better abundance distribution for quantitative analysis\n")
cat("3. More comparable to concentration-based measures\n")
cat("4. Standard metric used in PeptideAtlas publications\n")
cat("5. Enables better cross-study comparisons\n")
cat("\n")

cat("GENERATED FILES:\n")
cat("================\n")
cat("• Correlation analysis: correlation_plot.png\n")
cat("• Distribution comparison: distribution_plot.png\n")
cat("• Comprehensive analysis: comprehensive_plot.png\n")
cat("• Data files: peptideatlas_*.csv\n")
cat("• Summary report: quantification_analysis_summary.txt\n")

sink()

message("\nPeptideAtlas quantification comparison completed!")
message(sprintf("Results saved to: %s", output_dir))
message(sprintf("Plots saved to: %s", plot_dir))
message("\nSUMMARY:")
message(sprintf("- n_observations proteins: %d", nrow(peptideatlas_clean)))
message(sprintf("- norm_PSMs_per_100K proteins: %d", nrow(peptideatlas_clean)))
message("\nRECOMMENDATION: Use norm_PSMs_per_100K for normalized, quantitative analysis") 