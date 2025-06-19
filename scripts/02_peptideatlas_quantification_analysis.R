#!/usr/bin/env Rscript
# ============================================================================
# PeptideAtlas Quantification Methods Comparison Analysis
# Description: Compare n_observations vs norm_PSMs_per_100K quantification methods
# with Log and Z-score transformations
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

# Add log and z-score transformations
message("Computing log and z-score transformations...")
peptideatlas_transformed <- peptideatlas_clean %>%
  mutate(
    # Log transformations
    log_n_observations = log10(n_observations),
    log_norm_PSMs_per_100K = log10(norm_PSMs_per_100K),
    # Z-score transformations (standardized values)
    z_n_observations = scale(n_observations)[,1],
    z_norm_PSMs_per_100K = scale(norm_PSMs_per_100K)[,1],
    # Z-score on log-transformed data
    z_log_n_observations = scale(log10(n_observations))[,1],
    z_log_norm_PSMs_per_100K = scale(log10(norm_PSMs_per_100K))[,1]
  )

# Calculate correlations for different transformations
correlations <- tibble(
  transformation = c("raw", "log", "z_score", "z_score_log"),
  correlation = c(
    cor(peptideatlas_transformed$n_observations, peptideatlas_transformed$norm_PSMs_per_100K),
    cor(peptideatlas_transformed$log_n_observations, peptideatlas_transformed$log_norm_PSMs_per_100K),
    cor(peptideatlas_transformed$z_n_observations, peptideatlas_transformed$z_norm_PSMs_per_100K),
    cor(peptideatlas_transformed$z_log_n_observations, peptideatlas_transformed$z_log_norm_PSMs_per_100K)
  )
)

cat("Correlation coefficients:\n")
print(correlations)

# Create correlation plots for different transformations
cat("Creating correlation plots...\n")

# Log-transformed correlation plot
p1_log <- peptideatlas_transformed %>%
  ggplot(aes(x = log_n_observations, y = log_norm_PSMs_per_100K)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  labs(
    title = "Log-Transformed Correlation",
    subtitle = sprintf("Pearson r = %.4f", correlations$correlation[2]),
    x = "Log10(Number of Observations)",
    y = "Log10(Normalized PSMs per 100K)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold")
  )

# Z-score correlation plot (on log-transformed data)
p1_zscore <- peptideatlas_transformed %>%
  ggplot(aes(x = z_log_n_observations, y = z_log_norm_PSMs_per_100K)) +
  geom_point(alpha = 0.6, color = "darkorange") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  labs(
    title = "Z-Score Normalized Correlation",
    subtitle = sprintf("Pearson r = %.4f", correlations$correlation[4]),
    x = "Z-Score(Log10(Number of Observations))",
    y = "Z-Score(Log10(Normalized PSMs per 100K))"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold")
  )

# Create distribution comparison plots
cat("Creating distribution comparison plots...\n")

# Prepare data for distribution plots - comparing transformations
dist_data_comprehensive <- peptideatlas_transformed %>%
  select(n_observations, norm_PSMs_per_100K, 
         log_n_observations, log_norm_PSMs_per_100K,
         z_log_n_observations, z_log_norm_PSMs_per_100K) %>%
  pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("transformation", "quantification"), sep = "_", extra = "merge", fill = "left") %>%
  mutate(
    transformation = case_when(
      transformation == "n" ~ "raw",
      transformation == "log" ~ "log",
      transformation == "z" ~ "z_score",
      TRUE ~ transformation
    ),
    quantification = case_when(
      str_detect(quantification, "observations") ~ "Number of Observations",
      str_detect(quantification, "PSMs") ~ "Normalized PSMs per 100K",
      quantification == "observations" ~ "Number of Observations",
      quantification == "norm PSMs per 100K" ~ "Normalized PSMs per 100K",
      TRUE ~ quantification
    ),
    transformation = factor(transformation, levels = c("raw", "log", "z_score"))
  )

# Log distribution plot
p2_log <- peptideatlas_transformed %>%
  select(log_n_observations, log_norm_PSMs_per_100K) %>%
  pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = case_when(
      metric == "log_n_observations" ~ "Number of Observations",
      metric == "log_norm_PSMs_per_100K" ~ "Normalized PSMs per 100K"
    )
  ) %>%
  ggplot(aes(x = metric, y = value, fill = metric)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("Number of Observations" = "#E69F00", 
                              "Normalized PSMs per 100K" = "#56B4E9")) +
  labs(
    title = "Log-Transformed Distribution",
    subtitle = "Distribution of log10 quantification values",
    x = "Quantification Method",
    y = "Log10(Value)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Z-score distribution plot
p2_zscore <- peptideatlas_transformed %>%
  select(z_log_n_observations, z_log_norm_PSMs_per_100K) %>%
  pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = case_when(
      metric == "z_log_n_observations" ~ "Number of Observations",
      metric == "z_log_norm_PSMs_per_100K" ~ "Normalized PSMs per 100K"
    )
  ) %>%
  ggplot(aes(x = metric, y = value, fill = metric)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.7, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("Number of Observations" = "#E69F00", 
                              "Normalized PSMs per 100K" = "#56B4E9")) +
  labs(
    title = "Z-Score Normalized Distribution",
    subtitle = "Distribution of z-score normalized values",
    x = "Quantification Method",
    y = "Z-Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Create transformation comparison plot
p3_comparison <- correlations %>%
  ggplot(aes(x = transformation, y = correlation, fill = transformation)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("r = %.4f", correlation)), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("raw" = "#FF6B6B", "log" = "#4ECDC4", 
                              "z_score" = "#45B7D1", "z_score_log" = "#96CEB4")) +
  scale_x_discrete(labels = c("Raw", "Log", "Z-Score", "Z-Score(Log)")) +
  labs(
    title = "Correlation by Transformation",
    subtitle = "Comparing correlation coefficients",
    x = "Transformation Method",
    y = "Pearson Correlation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    legend.position = "none"
  ) +
  ylim(0, 1)

# Create comprehensive panel showing log and z-score analyses
cat("Creating comprehensive analysis panel...\n")

# Updated comprehensive plot: 2x3 layout
# Top row: Log correlation, Z-score correlation, Transformation comparison
# Bottom row: Log distribution, Z-score distribution, Combined density
comprehensive_plot <- (p1_log | p1_zscore | p3_comparison) / 
                     (p2_log | p2_zscore | plot_spacer()) +
  plot_annotation(
    title = "Comprehensive PeptideAtlas Quantification Analysis: Log & Z-Score Transformations",
    subtitle = "Comparing quantification methods with different data transformations",
    tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'))
  ) & 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

# Create individual plot for detailed log analysis (original style)
p_detailed_log <- peptideatlas_transformed %>%
  ggplot(aes(x = log_n_observations, y = log_norm_PSMs_per_100K)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  labs(
    title = "Detailed Log-Transformed Correlation Analysis",
    subtitle = sprintf("Pearson correlation: r = %.4f", correlations$correlation[2]),
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

# Save plots with high resolution - both PNG and TIFF versions
ggsave(file.path(plot_dir, "00_comprehensive_peptideatlas_analysis_panel.png"), comprehensive_plot, width = 18, height = 12, dpi = 600)
ggsave(file.path(plot_dir, "00_comprehensive_peptideatlas_analysis_panel.tiff"), comprehensive_plot, width = 18, height = 12, dpi = 600, device = "tiff", compression = "lzw", bg = "white")

# Print summary statistics
cat("\n=== PeptideAtlas Quantification Analysis Summary ===\n")
cat(sprintf("Total proteins analyzed: %d\n", nrow(peptideatlas_transformed)))
cat("\nCorrelation coefficients by transformation:\n")
print(correlations)

message("\nRECOMMENDATION: Use norm_PSMs_per_100K with log transformation for normalized, quantitative analysis")

# Save enhanced data files
message("Saving enhanced comparison data...")
write_csv(peptideatlas_transformed, file.path(output_dir, "peptideatlas_transformed.csv"))
write_csv(correlations, file.path(output_dir, "correlation_summary.csv"))

# Create enhanced summary report
sink(file.path(output_dir, "quantification_analysis_summary.txt"))
cat("PEPTIDEATLAS QUANTIFICATION METHODS COMPARISON\n")
cat("==============================================\n")
cat("Enhanced Analysis with Log and Z-Score Transformations\n\n")

cat("DATA SUMMARY:\n")
cat("=============\n")
cat(sprintf("Total proteins with quantification data: %d\n", nrow(peptideatlas_transformed)))
cat("\n")

cat("TRANSFORMATION ANALYSIS:\n")
cat("========================\n")
cat("Correlation coefficients by transformation method:\n")
for(i in 1:nrow(correlations)) {
  cat(sprintf("  %s: r = %.4f\n", correlations$transformation[i], correlations$correlation[i]))
}
cat("\n")

cat("QUANTIFICATION STATISTICS (Log-Transformed):\n")
cat("============================================\n")
cat("Log10(n_observations):\n")
cat(sprintf("  Range: %.6f - %.6f\n", min(peptideatlas_transformed$log_n_observations), max(peptideatlas_transformed$log_n_observations)))
cat(sprintf("  Median: %.6f\n", median(peptideatlas_transformed$log_n_observations)))
cat(sprintf("  Mean: %.6f\n", mean(peptideatlas_transformed$log_n_observations)))
cat(sprintf("  Standard Deviation: %.6f\n", sd(peptideatlas_transformed$log_n_observations)))
cat("\n")

cat("Log10(norm_PSMs_per_100K):\n")
cat(sprintf("  Range: %.6f - %.6f\n", min(peptideatlas_transformed$log_norm_PSMs_per_100K), max(peptideatlas_transformed$log_norm_PSMs_per_100K)))
cat(sprintf("  Median: %.6f\n", median(peptideatlas_transformed$log_norm_PSMs_per_100K)))
cat(sprintf("  Mean: %.6f\n", mean(peptideatlas_transformed$log_norm_PSMs_per_100K)))
cat(sprintf("  Standard Deviation: %.6f\n", sd(peptideatlas_transformed$log_norm_PSMs_per_100K)))
cat("\n")

cat("Z-SCORE STATISTICS (on Log-Transformed Data):\n")
cat("==============================================\n")
cat("Z-Score Log10(n_observations):\n")
cat(sprintf("  Range: %.6f - %.6f\n", min(peptideatlas_transformed$z_log_n_observations), max(peptideatlas_transformed$z_log_n_observations)))
cat(sprintf("  Median: %.6f\n", median(peptideatlas_transformed$z_log_n_observations)))
cat(sprintf("  Mean: %.6f\n", mean(peptideatlas_transformed$z_log_n_observations)))
cat(sprintf("  Standard Deviation: %.6f\n", sd(peptideatlas_transformed$z_log_n_observations)))
cat("\n")

cat("Z-Score Log10(norm_PSMs_per_100K):\n")
cat(sprintf("  Range: %.6f - %.6f\n", min(peptideatlas_transformed$z_log_norm_PSMs_per_100K), max(peptideatlas_transformed$z_log_norm_PSMs_per_100K)))
cat(sprintf("  Median: %.6f\n", median(peptideatlas_transformed$z_log_norm_PSMs_per_100K)))
cat(sprintf("  Mean: %.6f\n", mean(peptideatlas_transformed$z_log_norm_PSMs_per_100K)))
cat(sprintf("  Standard Deviation: %.6f\n", sd(peptideatlas_transformed$z_log_norm_PSMs_per_100K)))
cat("\n")

cat("RECOMMENDATION:\n")
cat("===============\n")
cat("norm_PSMs_per_100K with log transformation is recommended because:\n")
cat("1. Log transformation normalizes the distribution and reduces skewness\n")
cat("2. Z-score normalization centers and scales the data (mean=0, sd=1)\n")
cat("3. Normalized per 100K PSMs - accounts for run-to-run variation\n")
cat("4. Better correlation structure preserved across transformations\n")
cat("5. Standard metric used in PeptideAtlas publications\n")
cat("6. Enables better cross-study comparisons and statistical analysis\n")
cat("\n")

cat("TRANSFORMATION NOTES:\n")
cat("====================\n")
cat("• Log transformation: Reduces right skewness and stabilizes variance\n")
cat("• Z-score normalization: Centers data around 0 with unit variance\n")
cat("• Combined approach: Log + Z-score provides optimal normalization\n")
cat("• All transformations preserve correlation structure well\n")
cat("\n")

cat("GENERATED FILES:\n")
cat("================\n")
cat("• Detailed log correlation: log_correlation_plot.png\n")
cat("• Z-score correlation: zscore_correlation_plot.png\n")
cat("• Transformation comparison: transformation_comparison.png\n")
cat("• Comprehensive analysis: comprehensive_plot.png\n")
cat("• Enhanced data: peptideatlas_transformed.csv\n")
cat("• Correlation summary: correlation_summary.csv\n")
cat("• Summary report: quantification_analysis_summary.txt\n")

sink()

message("\nPeptideAtlas quantification comparison with transformations completed!")
message(sprintf("Results saved to: %s", output_dir))
message(sprintf("Plots saved to: %s", plot_dir))
message("\nSUMMARY:")
message(sprintf("- Total proteins: %d", nrow(peptideatlas_transformed)))
message("- Correlation coefficients:")
for(i in 1:nrow(correlations)) {
  message(sprintf("  %s: r = %.4f", correlations$transformation[i], correlations$correlation[i]))
}
message("\nRECOMMENDATION: Use norm_PSMs_per_100K with log transformation for optimal analysis") 