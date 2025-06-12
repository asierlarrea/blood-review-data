#!/usr/bin/env Rscript

# Load utilities and set up paths
source("scripts/utilities/load_packages.R")
ensure_output_dirs()


# PeptideAtlas Quantification Method Comparison
# Comparing n_observations vs norm_PSMs_per_100K
# Author: Data Analysis Pipeline
# Date: 2024

# Load required libraries
message("Loading required libraries...")
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(corrplot)
  library(gridExtra)
  library(RColorBrewer)
})

# Create output directory
output_dir <- "outputs/plots/PeptideAtlas_Comparison"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to extract gene symbols from descriptions
extract_gene_symbol <- function(description) {
  if(is.na(description) || description == "") {
    return(NA)
  }
  
  desc <- as.character(description)
  
  # Pattern 1: GN=GENE_NAME (UniProt format)
  if(grepl("GN=([A-Za-z0-9_-]+)", desc)) {
    gene <- gsub(".*GN=([A-Za-z0-9_-]+).*", "\\1", desc)
    return(toupper(gene))
  }
  
  # Pattern 2: \GName=GENE_NAME (PeptideAtlas format)
  if(grepl("\\\\GName=([A-Za-z0-9_-]+)", desc)) {
    gene <- gsub(".*\\\\GName=([A-Za-z0-9_-]+).*", "\\1", desc)
    return(toupper(gene))
  }
  
  # Pattern 3: Gene name at the beginning of description
  if(grepl("^([A-Z][A-Z0-9_-]+)\\s", desc)) {
    gene <- gsub("^([A-Z][A-Z0-9_-]+)\\s.*", "\\1", desc)
    if(nchar(gene) >= 2 && nchar(gene) <= 15) {
      return(toupper(gene))
    }
  }
  
  return(NA)
}

# Read and process PeptideAtlas data
message("Reading PeptideAtlas data...")
peptideatlas <- read.csv(get_data_path("PeptideAtlas.csv"), stringsAsFactors = FALSE)

message(paste("Total rows in PeptideAtlas:", nrow(peptideatlas)))
message("Column names:", paste(colnames(peptideatlas), collapse = ", "))

# Process data for both quantification methods
message("Processing data for quantification comparison...")

# Initialize vectors for both methods
valid_entries_n_obs <- logical(nrow(peptideatlas))
valid_entries_norm_psm <- logical(nrow(peptideatlas))
gene_symbols <- character(nrow(peptideatlas))
n_observations <- numeric(nrow(peptideatlas))
norm_psms <- numeric(nrow(peptideatlas))

for(i in 1:nrow(peptideatlas)) {
  gene_symbol <- NA
  
  # First try biosequence_gene_name column
  if(!is.na(peptideatlas$biosequence_gene_name[i]) && peptideatlas$biosequence_gene_name[i] != "") {
    gene_name <- as.character(peptideatlas$biosequence_gene_name[i])
    if(!grepl("[;,]", gene_name)) {  # Single gene only
      gene_symbol <- toupper(trimws(gene_name))
    }
  }
  
  # If not found, try extracting from description
  if(is.na(gene_symbol)) {
    gene_symbol <- extract_gene_symbol(peptideatlas$biosequence_description[i])
  }
  
  # Only keep valid gene symbols
  if(!is.na(gene_symbol) && 
     !grepl("^(UNIPROT_|GENE_|UNKNOWN_|ENSP_)", gene_symbol) &&
     nchar(gene_symbol) >= 2 && nchar(gene_symbol) <= 15) {
    
    gene_symbols[i] <- gene_symbol
    
    # Get n_observations value
    n_obs_value <- suppressWarnings(as.numeric(peptideatlas$n_observations[i]))
    if(length(n_obs_value) > 0 && !is.na(n_obs_value) && n_obs_value > 0) {
      n_observations[i] <- n_obs_value
      valid_entries_n_obs[i] <- TRUE
    }
    
    # Get norm_PSMs_per_100K value
    norm_psm_value <- suppressWarnings(as.numeric(peptideatlas$norm_PSMs_per_100K[i]))
    if(length(norm_psm_value) > 0 && !is.na(norm_psm_value) && norm_psm_value > 0) {
      norm_psms[i] <- norm_psm_value
      valid_entries_norm_psm[i] <- TRUE
    }
  }
}

# Create comparison datasets
n_obs_data <- data.frame(
  Gene = gene_symbols[valid_entries_n_obs],
  Value = n_observations[valid_entries_n_obs],
  LogValue = log10(n_observations[valid_entries_n_obs] + 1),
  Method = "n_observations",
  stringsAsFactors = FALSE
)

norm_psm_data <- data.frame(
  Gene = gene_symbols[valid_entries_norm_psm],
  Value = norm_psms[valid_entries_norm_psm],
  LogValue = log10(norm_psms[valid_entries_norm_psm] + 0.001),  # Small offset for log transformation
  Method = "norm_PSMs_per_100K",
  stringsAsFactors = FALSE
)

# Combine for comparison
all_data <- rbind(n_obs_data, norm_psm_data)

# Find genes present in both methods
genes_both <- intersect(n_obs_data$Gene, norm_psm_data$Gene)
message(paste("Genes with n_observations data:", nrow(n_obs_data)))
message(paste("Genes with norm_PSMs_per_100K data:", nrow(norm_psm_data)))
message(paste("Genes present in both methods:", length(genes_both)))

# Create paired comparison data
if(length(genes_both) > 0) {
  paired_data <- data.frame(
    Gene = genes_both,
    n_observations = sapply(genes_both, function(g) n_obs_data$Value[n_obs_data$Gene == g][1]),
    norm_PSMs_per_100K = sapply(genes_both, function(g) norm_psm_data$Value[norm_psm_data$Gene == g][1]),
    stringsAsFactors = FALSE
  )
  
  paired_data$log_n_observations <- log10(paired_data$n_observations + 1)
  paired_data$log_norm_PSMs_per_100K <- log10(paired_data$norm_PSMs_per_100K + 0.001)
}

# Summary statistics
message("\nSummary Statistics:")
message("n_observations:")
message(paste("  Range:", min(n_obs_data$Value), "-", max(n_obs_data$Value)))
message(paste("  Median:", median(n_obs_data$Value)))
message(paste("  Mean:", round(mean(n_obs_data$Value), 3)))

message("norm_PSMs_per_100K:")
message(paste("  Range:", min(norm_psm_data$Value), "-", max(norm_psm_data$Value)))
message(paste("  Median:", round(median(norm_psm_data$Value), 6)))
message(paste("  Mean:", round(mean(norm_psm_data$Value), 6)))

# Plot 1: Distribution comparison
message("Creating distribution comparison plot...")
tiff(file.path(output_dir, "quantification_methods_distribution.tiff"), 
     width = 14, height = 8, units = "in", res = 600)

p1 <- ggplot(all_data, aes(x = LogValue, fill = Method)) +
  geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
  scale_fill_manual(values = c("n_observations" = "blue", "norm_PSMs_per_100K" = "red")) +
  facet_wrap(~Method, scales = "free", ncol = 1) +
  theme_minimal() +
  labs(title = "PeptideAtlas Quantification Methods Comparison",
       subtitle = "Distribution of log-transformed values",
       x = "Log10(Value + offset)", y = "Frequency") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none")

print(p1)
dev.off()

# Plot 2: Correlation between methods (if paired data exists)
if(length(genes_both) > 0) {
  message("Creating correlation plot...")
  
  # Calculate correlation
  correlation <- cor(paired_data$n_observations, paired_data$norm_PSMs_per_100K, 
                    use = "complete.obs")
  log_correlation <- cor(paired_data$log_n_observations, paired_data$log_norm_PSMs_per_100K, 
                        use = "complete.obs")
  
  message(paste("Correlation (raw values):", round(correlation, 4)))
  message(paste("Correlation (log values):", round(log_correlation, 4)))
  
  tiff(file.path(output_dir, "quantification_methods_correlation.tiff"), 
       width = 12, height = 10, units = "in", res = 600)
  
  # Create correlation plots
  p2a <- ggplot(paired_data, aes(x = n_observations, y = norm_PSMs_per_100K)) +
    geom_point(alpha = 0.6, color = "darkblue") +
    geom_smooth(method = "lm", color = "red") +
    theme_minimal() +
    labs(title = "Raw Values Correlation",
         subtitle = paste("Pearson r =", round(correlation, 4)),
         x = "n_observations", y = "norm_PSMs_per_100K") +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  p2b <- ggplot(paired_data, aes(x = log_n_observations, y = log_norm_PSMs_per_100K)) +
    geom_point(alpha = 0.6, color = "darkgreen") +
    geom_smooth(method = "lm", color = "red") +
    theme_minimal() +
    labs(title = "Log-transformed Values Correlation",
         subtitle = paste("Pearson r =", round(log_correlation, 4)),
         x = "Log10(n_observations + 1)", y = "Log10(norm_PSMs_per_100K + 0.001)") +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  grid.arrange(p2a, p2b, ncol = 2, 
               top = "PeptideAtlas Quantification Methods: n_observations vs norm_PSMs_per_100K")
  
  dev.off()
}

# Plot 3: Dynamic range comparison
message("Creating dynamic range comparison...")
tiff(file.path(output_dir, "quantification_methods_boxplot.tiff"), 
     width = 12, height = 8, units = "in", res = 600)

p3 <- ggplot(all_data, aes(x = Method, y = LogValue, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  geom_violin(alpha = 0.3) +
  scale_fill_manual(values = c("n_observations" = "blue", "norm_PSMs_per_100K" = "red")) +
  theme_minimal() +
  labs(title = "Dynamic Range Comparison",
       subtitle = "Distribution of quantification values (log scale)",
       x = "Quantification Method", y = "Log10(Value + offset)") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)
dev.off()

# Plot 4: Top proteins comparison
message("Creating top proteins comparison...")
top_n_obs <- n_obs_data %>% 
  arrange(desc(Value)) %>% 
  head(20) %>%
  mutate(Rank = 1:20, Method = "n_observations")

top_norm_psm <- norm_psm_data %>% 
  arrange(desc(Value)) %>% 
  head(20) %>%
  mutate(Rank = 1:20, Method = "norm_PSMs_per_100K")

top_combined <- rbind(
  top_n_obs %>% select(Gene, Value, Rank, Method),
  top_norm_psm %>% select(Gene, Value, Rank, Method)
)

tiff(file.path(output_dir, "top_proteins_comparison.tiff"), 
     width = 14, height = 10, units = "in", res = 600)

p4 <- ggplot(top_combined, aes(x = Rank, y = Value, color = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("n_observations" = "blue", "norm_PSMs_per_100K" = "red")) +
  facet_wrap(~Method, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(title = "Top 20 Proteins by Quantification Method",
       subtitle = "Ranking comparison between methods",
       x = "Rank", y = "Quantification Value") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.position = "bottom")

print(p4)
dev.off()

# Save comparison data to CSV
message("Saving comparison data...")
write.csv(n_obs_data, file.path(output_dir, "peptideatlas_n_observations.csv"), row.names = FALSE)
write.csv(norm_psm_data, file.path(output_dir, "peptideatlas_norm_PSMs_per_100K.csv"), row.names = FALSE)

if(length(genes_both) > 0) {
  write.csv(paired_data, file.path(output_dir, "peptideatlas_paired_comparison.csv"), row.names = FALSE)
}

# Create summary report
sink(file.path(output_dir, "quantification_comparison_summary.txt"))
cat("PeptideAtlas Quantification Methods Comparison\n")
cat("============================================\n\n")

cat("DATA SUMMARY:\n")
cat(paste("Total proteins with n_observations:", nrow(n_obs_data), "\n"))
cat(paste("Total proteins with norm_PSMs_per_100K:", nrow(norm_psm_data), "\n"))
cat(paste("Proteins present in both methods:", length(genes_both), "\n\n"))

cat("QUANTIFICATION STATISTICS:\n")
cat("n_observations:\n")
cat(paste("  Range:", min(n_obs_data$Value), "-", max(n_obs_data$Value), "\n"))
cat(paste("  Median:", median(n_obs_data$Value), "\n"))
cat(paste("  Mean:", round(mean(n_obs_data$Value), 3), "\n"))
cat(paste("  Standard Deviation:", round(sd(n_obs_data$Value), 3), "\n\n"))

cat("norm_PSMs_per_100K:\n")
cat(paste("  Range:", min(norm_psm_data$Value), "-", max(norm_psm_data$Value), "\n"))
cat(paste("  Median:", round(median(norm_psm_data$Value), 6), "\n"))
cat(paste("  Mean:", round(mean(norm_psm_data$Value), 6), "\n"))
cat(paste("  Standard Deviation:", round(sd(norm_psm_data$Value), 6), "\n\n"))

if(length(genes_both) > 0) {
  cat("CORRELATION ANALYSIS:\n")
  cat(paste("Pearson correlation (raw values):", round(correlation, 4), "\n"))
  cat(paste("Pearson correlation (log values):", round(log_correlation, 4), "\n\n"))
}

cat("RECOMMENDATION:\n")
cat("norm_PSMs_per_100K is recommended because:\n")
cat("1. Normalized per 100K PSMs - accounts for run-to-run variation\n")
cat("2. Better dynamic range for quantitative analysis\n")
cat("3. More comparable to concentration-based measures\n")
cat("4. Standard metric used in PeptideAtlas publications\n")

sink()

message("\nPeptideAtlas quantification comparison completed!")
message(paste("Results saved to:", output_dir))
message("\nSUMMARY:")
message(paste("- n_observations proteins:", nrow(n_obs_data)))
message(paste("- norm_PSMs_per_100K proteins:", nrow(norm_psm_data)))
if(length(genes_both) > 0) {
  message(paste("- Correlation between methods:", round(log_correlation, 4)))
}
message("\nRECOMMENDATION: Use norm_PSMs_per_100K for normalized, quantitative analysis") 
