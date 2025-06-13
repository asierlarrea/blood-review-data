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

# Read PeptideAtlas data
message("Reading and processing PeptideAtlas data...")
peptideatlas <- read_csv("data/raw/peptideatlas/peptideatlas.csv", show_col_types = FALSE)

message(sprintf("Total rows in PeptideAtlas: %d", nrow(peptideatlas)))
message("Column names:", paste(colnames(peptideatlas), collapse = ", "))

# Map accessions to gene symbols using modern approach
peptideatlas$gene <- convert_to_gene_symbol(peptideatlas$biosequence_accession, force_mapping = force_mapping)

# Create datasets for both quantification methods
message("Processing quantification methods...")

# Method 1: n_observations
n_obs_data <- peptideatlas %>%
  filter(!is.na(gene), !is.na(n_observations), n_observations > 0) %>%
  select(gene, value = n_observations) %>%
  mutate(
    log_value = log10(value + 1),
    method = "n_observations"
  )

# Method 2: norm_PSMs_per_100K  
norm_psm_data <- peptideatlas %>%
  filter(!is.na(gene), !is.na(norm_PSMs_per_100K), norm_PSMs_per_100K > 0) %>%
  select(gene, value = norm_PSMs_per_100K) %>%
  mutate(
    log_value = log10(value + 0.001),  # Small offset for log transformation
    method = "norm_PSMs_per_100K"
  )

# Combine for comparison
all_data <- rbind(n_obs_data, norm_psm_data)

# Find genes present in both methods
genes_both <- intersect(n_obs_data$gene, norm_psm_data$gene)

message(sprintf("Genes with n_observations data: %d", nrow(n_obs_data)))
message(sprintf("Genes with norm_PSMs_per_100K data: %d", nrow(norm_psm_data)))
message(sprintf("Genes present in both methods: %d", length(genes_both)))

# Create paired comparison data
paired_data <- NULL
if(length(genes_both) > 0) {
  paired_data <- data.frame(
    gene = genes_both,
    n_observations = sapply(genes_both, function(g) n_obs_data$value[n_obs_data$gene == g][1]),
    norm_PSMs_per_100K = sapply(genes_both, function(g) norm_psm_data$value[norm_psm_data$gene == g][1]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      log_n_observations = log10(n_observations + 1),
      log_norm_PSMs_per_100K = log10(norm_PSMs_per_100K + 0.001)
    )
}

# Calculate summary statistics
stats_summary <- list(
  n_obs_count = nrow(n_obs_data),
  norm_psm_count = nrow(norm_psm_data),
  paired_count = length(genes_both),
  n_obs_stats = list(
    min = min(n_obs_data$value),
    median = median(n_obs_data$value),
    mean = mean(n_obs_data$value),
    max = max(n_obs_data$value),
    sd = sd(n_obs_data$value)
  ),
  norm_psm_stats = list(
    min = min(norm_psm_data$value),
    median = median(norm_psm_data$value),
    mean = mean(norm_psm_data$value),
    max = max(norm_psm_data$value),
    sd = sd(norm_psm_data$value)
  )
)

# Calculate correlations if paired data exists
if(!is.null(paired_data) && nrow(paired_data) > 0) {
  stats_summary$correlation_raw <- cor(paired_data$n_observations, paired_data$norm_PSMs_per_100K, use = "complete.obs")
  stats_summary$correlation_log <- cor(paired_data$log_n_observations, paired_data$log_norm_PSMs_per_100K, use = "complete.obs")
}

# Create visualizations
message("Creating visualizations...")

# 1. Distribution comparison plot
p1 <- ggplot(all_data, aes(x = log_value, fill = method)) +
  geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
  scale_fill_manual(values = c("n_observations" = "#2E86AB", "norm_PSMs_per_100K" = "#A23B72")) +
  facet_wrap(~method, scales = "free", ncol = 1, labeller = labeller(method = c(
    "n_observations" = "n_observations",
    "norm_PSMs_per_100K" = "norm_PSMs_per_100K"
  ))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "PeptideAtlas Quantification Methods Distribution",
    subtitle = "Comparison of log-transformed quantification values",
    x = "Log10(Value + offset)",
    y = "Frequency"
  )

ggsave(file.path(plot_dir, "quantification_methods_distribution.png"), p1, 
       width = 10, height = 8, dpi = 300, bg = "white")

# 2. Dynamic range comparison (box plot)
p2 <- ggplot(all_data, aes(x = method, y = log_value, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  geom_violin(alpha = 0.3) +
  scale_fill_manual(values = c("n_observations" = "#2E86AB", "norm_PSMs_per_100K" = "#A23B72")) +
  scale_x_discrete(labels = c("n_observations" = "n_observations", "norm_PSMs_per_100K" = "norm_PSMs_per_100K")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Dynamic Range Comparison",
    subtitle = "Distribution of quantification values (log scale)",
    x = "Quantification Method",
    y = "Log10(Value + offset)"
  )

ggsave(file.path(plot_dir, "quantification_methods_boxplot.png"), p2, 
       width = 10, height = 8, dpi = 300, bg = "white")

# 3. Correlation plots (if paired data exists)
if(!is.null(paired_data) && nrow(paired_data) > 0) {
  # Raw values correlation
  p3a <- ggplot(paired_data, aes(x = n_observations, y = norm_PSMs_per_100K)) +
    geom_point(alpha = 0.6, color = "#2E86AB") +
    geom_smooth(method = "lm", color = "#A23B72", linewidth = 1) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    ) +
    labs(
      title = "Raw Values Correlation",
      subtitle = sprintf("Pearson r = %.3f", stats_summary$correlation_raw),
      x = "n_observations",
      y = "norm_PSMs_per_100K"
    )
  
  # Log-transformed values correlation
  p3b <- ggplot(paired_data, aes(x = log_n_observations, y = log_norm_PSMs_per_100K)) +
    geom_point(alpha = 0.6, color = "#2E86AB") +
    geom_smooth(method = "lm", color = "#A23B72", linewidth = 1) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    ) +
    labs(
      title = "Log-transformed Values Correlation",
      subtitle = sprintf("Pearson r = %.3f", stats_summary$correlation_log),
      x = "Log10(n_observations + 1)",
      y = "Log10(norm_PSMs_per_100K + 0.001)"
    )
  
  # Combined correlation plot
  p3 <- p3a | p3b
  p3 <- p3 + plot_annotation(
    title = "PeptideAtlas Quantification Methods Correlation",
    subtitle = "Relationship between n_observations and norm_PSMs_per_100K",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )
  
  ggsave(file.path(plot_dir, "quantification_methods_correlation.png"), p3, 
         width = 16, height = 8, dpi = 300, bg = "white")
}

# 4. Create comprehensive combined plot
if(!is.null(paired_data) && nrow(paired_data) > 0) {
  # Create a 3-panel comprehensive plot with correlation
  combined_plot <- (p1 | p2) / p3
  combined_plot <- combined_plot + 
    plot_annotation(
      title = "Comprehensive PeptideAtlas Quantification Methods Analysis",
      subtitle = "Distribution comparison, dynamic range analysis, and correlation",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5)
      )
    )
} else {
  # Fallback to 2-panel plot if no correlation data
  combined_plot <- p1 | p2
  combined_plot <- combined_plot + 
    plot_annotation(
      title = "Comprehensive PeptideAtlas Quantification Methods Analysis",
      subtitle = "Distribution comparison and dynamic range analysis",
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5)
      )
    )
}

ggsave(file.path(plot_dir, "peptideatlas_quantification_comprehensive.png"), combined_plot, 
       width = 16, height = 12, dpi = 300, bg = "white")

message("Plots saved to:", plot_dir)

# Save data files
message("Saving comparison data...")
write_csv(n_obs_data, file.path(output_dir, "peptideatlas_n_observations.csv"))
write_csv(norm_psm_data, file.path(output_dir, "peptideatlas_norm_PSMs_per_100K.csv"))

if(!is.null(paired_data) && nrow(paired_data) > 0) {
  write_csv(paired_data, file.path(output_dir, "peptideatlas_paired_comparison.csv"))
}

# Save summary statistics
write_csv(
  data.frame(
    Metric = c("n_observations_count", "norm_PSMs_per_100K_count", "paired_count"),
    Value = c(stats_summary$n_obs_count, stats_summary$norm_psm_count, stats_summary$paired_count)
  ),
  file.path(output_dir, "quantification_summary_counts.csv")
)

# Create summary report
sink(file.path(output_dir, "quantification_analysis_summary.txt"))
cat("PEPTIDEATLAS QUANTIFICATION METHODS COMPARISON\n")
cat("=============================================\n\n")

cat("DATA SUMMARY:\n")
cat("=============\n")
cat(sprintf("Total proteins with n_observations: %d\n", stats_summary$n_obs_count))
cat(sprintf("Total proteins with norm_PSMs_per_100K: %d\n", stats_summary$norm_psm_count))
cat(sprintf("Proteins present in both methods: %d\n", stats_summary$paired_count))
cat("\n")

cat("QUANTIFICATION STATISTICS:\n")
cat("==========================\n")
cat("n_observations:\n")
cat(sprintf("  Range: %d - %d\n", stats_summary$n_obs_stats$min, stats_summary$n_obs_stats$max))
cat(sprintf("  Median: %.0f\n", stats_summary$n_obs_stats$median))
cat(sprintf("  Mean: %.3f\n", stats_summary$n_obs_stats$mean))
cat(sprintf("  Standard Deviation: %.3f\n", stats_summary$n_obs_stats$sd))
cat("\n")

cat("norm_PSMs_per_100K:\n")
cat(sprintf("  Range: %.6f - %.6f\n", stats_summary$norm_psm_stats$min, stats_summary$norm_psm_stats$max))
cat(sprintf("  Median: %.6f\n", stats_summary$norm_psm_stats$median))
cat(sprintf("  Mean: %.6f\n", stats_summary$norm_psm_stats$mean))
cat(sprintf("  Standard Deviation: %.6f\n", stats_summary$norm_psm_stats$sd))
cat("\n")

if(!is.null(stats_summary$correlation_raw)) {
  cat("CORRELATION ANALYSIS:\n")
  cat("====================\n")
  cat(sprintf("Pearson correlation (raw values): %.4f\n", stats_summary$correlation_raw))
  cat(sprintf("Pearson correlation (log values): %.4f\n", stats_summary$correlation_log))
  cat("\n")
}

cat("RECOMMENDATION:\n")
cat("===============\n")
cat("norm_PSMs_per_100K is recommended because:\n")
cat("1. Normalized per 100K PSMs - accounts for run-to-run variation\n")
cat("2. Better dynamic range for quantitative analysis\n")
cat("3. More comparable to concentration-based measures\n")
cat("4. Standard metric used in PeptideAtlas publications\n")
cat("5. Enables better cross-study comparisons\n")
cat("\n")

cat("GENERATED FILES:\n")
cat("================\n")
cat("• Distribution comparison: quantification_methods_distribution.png\n")
cat("• Dynamic range comparison: quantification_methods_boxplot.png\n")
if(!is.null(paired_data) && nrow(paired_data) > 0) {
  cat("• Correlation analysis: quantification_methods_correlation.png\n")
}
cat("• Comprehensive analysis: peptideatlas_quantification_comprehensive.png\n")
cat("• Data files: peptideatlas_*.csv\n")
cat("• Summary report: quantification_analysis_summary.txt\n")

sink()

message("\nPeptideAtlas quantification comparison completed!")
message(sprintf("Results saved to: %s", output_dir))
message(sprintf("Plots saved to: %s", plot_dir))
message("\nSUMMARY:")
message(sprintf("- n_observations proteins: %d", stats_summary$n_obs_count))
message(sprintf("- norm_PSMs_per_100K proteins: %d", stats_summary$norm_psm_count))
if(!is.null(stats_summary$correlation_log)) {
  message(sprintf("- Correlation between methods: %.4f", stats_summary$correlation_log))
}
message("\nRECOMMENDATION: Use norm_PSMs_per_100K for normalized, quantitative analysis") 