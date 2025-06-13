#!/usr/bin/env Rscript
#
# Analysis of Plasma Proteins from Multiple Data Sources
# Script: 01_plasma_protein_analysis.R
# 
# Purpose: Analyze the total number of genes/proteins quantified in plasma by different sources:
# - PeptideAtlas (UniProt accessions)
# - HPA (Human Protein Atlas) - different technologies: MS, PEA, Immunoassays (Gene names)
# - GPMDB (Ensembl Protein accessions)
# - PAXDB (Ensembl Protein accessions)
#
# Updated to use consistent gene name mapping for accurate comparisons

# Load utilities and set up output directories
source("scripts/utilities/load_packages.R")
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "ggthemes", "patchwork", "UpSetR", "ggupset")
load_packages(required_packages)

# Parse command line arguments (simple approach)
args <- commandArgs(trailingOnly = TRUE)
force_mapping <- "--force-mapping" %in% args

# Load gene mapping utility
source("scripts/data_processing/simple_id_mapping.R")

# Load gene deduplication utility
source("scripts/utilities/gene_deduplication.R")

# Set output directory
output_dir <- get_output_path("", subdir = "plasma_protein")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read and process data from each source
message("Reading and processing data from each source...")

# 1. PeptideAtlas
message("Processing PeptideAtlas data...")
peptideatlas_raw <- read_csv("data/raw/peptideatlas/peptideatlas.csv", show_col_types = FALSE)
peptideatlas_raw$gene <- convert_to_gene_symbol(peptideatlas_raw$biosequence_accession, force_mapping = force_mapping)
peptideatlas_raw <- peptideatlas_raw %>% filter(!is.na(norm_PSMs_per_100K))

# Deduplicate genes using median quantification values
peptideatlas <- deduplicate_genes(peptideatlas_raw, "gene", "norm_PSMs_per_100K", 
                                additional_cols = c("biosequence_accession"), 
                                aggregation_method = "median")

# 2. HPA MS
message("Processing HPA MS data...")
hpa_ms_raw <- read_csv("data/raw/hpa/hpa_ms.csv", show_col_types = FALSE, skip = 1)
hpa_ms_raw <- hpa_ms_raw %>% rename(gene = Gene, expr = Concentration)

# Deduplicate genes using median quantification values
hpa_ms <- deduplicate_genes(hpa_ms_raw, "gene", "expr", aggregation_method = "median")

# 3. HPA PEA
message("Processing HPA PEA data...")
hpa_pea_raw <- read_csv("data/raw/hpa/hpa_pea.csv", show_col_types = FALSE)
hpa_pea_raw <- hpa_pea_raw %>% rename(gene = Gene, expr = `Variation between individuals`)

# Deduplicate genes using median quantification values
hpa_pea <- deduplicate_genes(hpa_pea_raw, "gene", "expr", aggregation_method = "median")

# 4. HPA Immunoassay
message("Processing HPA Immunoassay data...")
hpa_imm_raw <- read_csv("data/raw/hpa/hpa_immunoassay_plasma.csv", show_col_types = FALSE)
hpa_imm_raw <- hpa_imm_raw %>% rename(gene = Gene, expr = Concentration)

# Deduplicate genes using median quantification values
hpa_imm <- deduplicate_genes(hpa_imm_raw, "gene", "expr", aggregation_method = "median")

# 5. GPMDB
message("Processing GPMDB data...")
gpmdb_raw <- read_csv("data/raw/gpmdb/gpmdb_plasma.csv", show_col_types = FALSE)
gpmdb_raw$gene <- stringr::str_extract(gpmdb_raw$description, "[A-Z0-9]+(?=,| |$)")

# Deduplicate genes using median quantification values (assuming 'total' column exists)
# If no quantification column, we'll filter for unique genes only
if ("total" %in% colnames(gpmdb_raw)) {
  gpmdb <- deduplicate_genes(gpmdb_raw, "gene", "total", aggregation_method = "median")
} else {
  # Just remove duplicate genes if no quantification column
  gpmdb <- gpmdb_raw %>% 
    filter(!is.na(gene) & gene != "") %>%
    distinct(gene, .keep_all = TRUE)
}

# 6. PAXDB
message("Processing PAXDB data...")
paxdb_raw <- read_csv("data/raw/paxdb/paxdb_plasma.csv", show_col_types = FALSE)
paxdb_raw$ensp <- stringr::str_replace(paxdb_raw$string_external_id, "^9606\\.", "")
paxdb_raw$gene <- convert_to_gene_symbol(paxdb_raw$ensp, force_mapping = force_mapping)

# Deduplicate genes using median quantification values (assuming 'abundance' column exists)
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

# Create summary statistics
message("Creating summary statistics...")
stats_summary <- list(
  peptideatlas = length(unique(peptideatlas$gene)),
  hpa_ms = length(unique(hpa_ms$gene)),
  hpa_pea = length(unique(hpa_pea$gene)),
  hpa_immunoassay = length(unique(hpa_imm$gene)),
  gpmdb = length(unique(gpmdb$gene)),
  paxdb = length(unique(paxdb$gene))
)

# Calculate total unique genes across sources
all_genes <- unique(c(
  peptideatlas$gene,
  hpa_ms$gene,
  hpa_pea$gene,
  hpa_imm$gene,
  gpmdb$gene,
  paxdb$gene
))
stats_summary$total_across_sources <- length(all_genes)

# Calculate MS technologies total (including PAXDB)
ms_genes <- unique(c(
  peptideatlas$gene,
  hpa_ms$gene,
  gpmdb$gene,
  paxdb$gene
))
stats_summary$ms_technologies <- length(ms_genes)

# Calculate HPA total (all technologies)
hpa_genes <- unique(c(
  hpa_ms$gene,
  hpa_pea$gene,
  hpa_imm$gene
))
stats_summary$hpa_total <- length(hpa_genes)

# Save summary statistics
write.csv(
  data.frame(
    Source = names(stats_summary),
    Gene_Count = unlist(stats_summary)
  ),
  file.path(output_dir, "plasma_protein_counts_summary.csv"),
  row.names = FALSE
)

# Create visualizations
message("Creating visualizations...")

# Set up plot output directory
plot_dir <- get_output_path("01_plasma_protein_analysis", subdir = "plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# 1. Bar plot by individual sources
source_data <- data.frame(
  Source = c("PeptideAtlas", "HPA MS", "HPA PEA", "HPA Immunoassay", "GPMDB", "PAXDB"),
  Count = c(stats_summary$peptideatlas, stats_summary$hpa_ms, stats_summary$hpa_pea, 
            stats_summary$hpa_immunoassay, stats_summary$gpmdb, stats_summary$paxdb),
  Technology = c("MS", "MS", "PEA", "Immunoassay", "MS", "MS"),
  Database = c("PeptideAtlas", "HPA", "HPA", "HPA", "GPMDB", "PAXDB")
)

p1 <- ggplot(source_data, aes(x = reorder(Source, Count), y = Count, fill = Technology)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("MS" = "#2E86AB", "PEA" = "#A23B72", 
                               "Immunoassay" = "#F18F01")) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Plasma Proteins Quantified by Data Source",
    subtitle = "Number of unique genes detected in each database",
    x = "Data Source",
    y = "Number of Genes",
    fill = "Technology"
  )

ggsave(file.path(plot_dir, "plasma_proteins_by_source.png"), p1, 
       width = 10, height = 6, dpi = 300, bg = "white")

# 2. Bar plot by technology grouping (comprehensive breakdown)
tech_data <- data.frame(
  Technology = c("Mass Spectrometry", "PEA", "Immunoassay", "Total Across Sources"),
  Count = c(stats_summary$ms_technologies, stats_summary$hpa_pea, stats_summary$hpa_immunoassay, 
            stats_summary$total_across_sources),
  Type = c("Technology", "Technology", "Technology", "Overall")
)

p2 <- ggplot(tech_data, aes(x = reorder(Technology, Count), y = Count, fill = Type)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("Technology" = "#2E86AB", "Overall" = "#C73E1D")) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Plasma Proteins by Technology Classification",
    subtitle = "MS includes PeptideAtlas, HPA MS, GPMDB, and PAXDB; PEA and Immunoassay from HPA",
    x = "Technology Category",
    y = "Number of Genes",
    fill = "Category Type"
  )

ggsave(file.path(plot_dir, "plasma_proteins_by_technology.png"), p2, 
       width = 10, height = 6, dpi = 300, bg = "white")

# 3. Create violin plot showing quantification value distributions
message("Creating quantification distribution violin plot...")

# Prepare combined data for violin plot
violin_data <- data.frame()

# Add PeptideAtlas data
if (nrow(peptideatlas) > 0) {
  violin_data <- rbind(violin_data, data.frame(
    Database = "PeptideAtlas",
    Technology = "MS",
    log_value = log10(peptideatlas$norm_PSMs_per_100K + 1),
    stringsAsFactors = FALSE
  ))
}

# Add HPA MS data
if (nrow(hpa_ms) > 0) {
  violin_data <- rbind(violin_data, data.frame(
    Database = "HPA MS",
    Technology = "MS", 
    log_value = log10(hpa_ms$expr + 1),
    stringsAsFactors = FALSE
  ))
}

# Add HPA PEA data
if (nrow(hpa_pea) > 0) {
  violin_data <- rbind(violin_data, data.frame(
    Database = "HPA PEA",
    Technology = "PEA",
    log_value = log10(hpa_pea$expr + 1),
    stringsAsFactors = FALSE
  ))
}

# Add HPA Immunoassay data
if (nrow(hpa_imm) > 0) {
  violin_data <- rbind(violin_data, data.frame(
    Database = "HPA Immunoassay", 
    Technology = "Immunoassay",
    log_value = log10(hpa_imm$expr + 1),
    stringsAsFactors = FALSE
  ))
}

# Add GPMDB data (if total column exists)
if (nrow(gpmdb) > 0 && "total" %in% colnames(gpmdb)) {
  violin_data <- rbind(violin_data, data.frame(
    Database = "GPMDB",
    Technology = "MS",
    log_value = log10(gpmdb$total + 1),
    stringsAsFactors = FALSE
  ))
}

# Add PAXDB data (if abundance column exists)
if (nrow(paxdb) > 0 && "abundance" %in% colnames(paxdb)) {
  violin_data <- rbind(violin_data, data.frame(
    Database = "PAXDB",
    Technology = "MS", 
    log_value = log10(paxdb$abundance + 1),
    stringsAsFactors = FALSE
  ))
}

# Create violin plot with improved aesthetics
p3 <- ggplot(violin_data, aes(x = reorder(Database, log_value, FUN = median), 
                              y = log_value, fill = Technology)) +
  # Violin plots with better scaling and smoothing
  geom_violin(alpha = 0.8, 
              scale = "area",           # Equal area scaling for better comparison
              trim = TRUE,              # Trim to data range
              adjust = 1.2,             # Slightly smoother curves
              draw_quantiles = c(0.25, 0.5, 0.75),  # Add quartile lines
              color = "#150b0b",          # White outline for better definition
              size = 0.3) +
  # Refined boxplot overlay with matching colors
  geom_boxplot(aes(fill = Technology),
               width = 0.15, 
               alpha = 0.6,             # More transparent to show violin underneath
               outlier.size = 1.2,      # Larger outlier points
               outlier.alpha = 0.7,     # More opaque outliers
               outlier.color = "gray20", # Dark color for outliers
               color = "gray30",        # Dark outline for box and whiskers
               size = 0.5) +
  # Add median points for emphasis
  stat_summary(fun = median, 
               geom = "point", 
               shape = 21, 
               size = 2.5, 
               fill = "white", 
               color = "gray20", 
               stroke = 1) +
  coord_flip() +
  scale_fill_manual(values = c("MS" = "#2E86AB", "PEA" = "#A23B72", 
                               "Immunoassay" = "#F18F01")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray90", fill = NA, size = 0.5)
  ) +
  labs(
    title = "Protein Abundance Distribution by Data Source",
    subtitle = "Log10-transformed quantification values with quartiles (after gene deduplication)",
    x = "Data Source",
    y = "Log10(Quantification Value + 1)",
    fill = "Technology"
  )

ggsave(file.path(plot_dir, "plasma_proteins_abundance_distributions.png"), p3, 
       width = 10, height = 6, dpi = 300, bg = "white")

# 4. Create UpSet plot showing gene overlaps between databases
message("Creating UpSet plot for gene overlaps...")

# Create gene presence/absence matrix
all_unique_genes <- unique(c(
  peptideatlas$gene,
  hpa_ms$gene,
  hpa_pea$gene,
  hpa_imm$gene,
  gpmdb$gene,
  paxdb$gene
))

# Remove any NA genes
all_unique_genes <- all_unique_genes[!is.na(all_unique_genes) & all_unique_genes != ""]

# Create binary matrix
upset_matrix <- data.frame(
  Gene = all_unique_genes,
  PeptideAtlas = as.integer(all_unique_genes %in% peptideatlas$gene),
  HPA_MS = as.integer(all_unique_genes %in% hpa_ms$gene),
  HPA_PEA = as.integer(all_unique_genes %in% hpa_pea$gene),
  HPA_Immunoassay = as.integer(all_unique_genes %in% hpa_imm$gene),
  GPMDB = as.integer(all_unique_genes %in% gpmdb$gene),
  PAXDB = as.integer(all_unique_genes %in% paxdb$gene),
  stringsAsFactors = FALSE
)

# Save the classic UpSet plot as PNG using UpSetR
upset_file <- file.path(plot_dir, "gene_overlaps_upset.png")
png(upset_file, width = 12, height = 8, units = "in", res = 300)

UpSetR::upset(
  upset_matrix,
  sets = c("PeptideAtlas", "HPA_MS", "HPA_PEA", "HPA_Immunoassay", "GPMDB", "PAXDB"),
  order.by = "freq",
  nsets = 6,
  nintersects = 25,
  point.size = 3.5,
  line.size = 1.2,
  mainbar.y.label = "Gene Intersection Size",
  sets.x.label = "Genes per Database",
  text.scale = c(1.3, 1.3, 1.2, 1.2, 1.4, 1.2),  # Increase text sizes
  mb.ratio = c(0.6, 0.4),
  keep.order = TRUE,
  main.bar.color = "#2E86AB",
  sets.bar.color = c("#2E86AB", "#2E86AB", "#A23B72", "#F18F01", "#2E86AB", "#2E86AB"),  # Match technology colors
  matrix.color = "gray30",
  shade.color = "gray90",
  set_size.show = TRUE,
  set_size.scale_max = max(colSums(upset_matrix[,-1])) * 1.1
)

dev.off()

# Create a ggplot2-compatible UpSet plot using ggupset
# Prepare data in the format ggupset expects (long format with list columns)
upset_data_long <- data.frame()

# Add genes from each database
databases <- list(
  "PeptideAtlas" = peptideatlas$gene[!is.na(peptideatlas$gene)],
  "HPA MS" = hpa_ms$gene[!is.na(hpa_ms$gene)],
  "HPA PEA" = hpa_pea$gene[!is.na(hpa_pea$gene)],
  "HPA Immunoassay" = hpa_imm$gene[!is.na(hpa_imm$gene)],
  "GPMDB" = gpmdb$gene[!is.na(gpmdb$gene)],
  "PAXDB" = paxdb$gene[!is.na(paxdb$gene)]
)

# Create a data frame where each row is a gene and columns indicate which databases contain it
all_genes_upset <- unique(unlist(databases))
upset_df <- data.frame(gene = all_genes_upset, stringsAsFactors = FALSE)

for(db_name in names(databases)) {
  upset_df[[db_name]] <- upset_df$gene %in% databases[[db_name]]
}

# Convert to the format ggupset needs (list column)
upset_df_tidy <- upset_df %>%
  pivot_longer(cols = -gene, names_to = "database", values_to = "present") %>%
  filter(present) %>%
  group_by(gene) %>%
  summarise(databases = list(database), .groups = "drop")

# Create UpSet plot using ggupset
p4 <- upset_df_tidy %>%
  ggplot(aes(x = databases)) +
  geom_bar(fill = "#2E86AB", alpha = 0.8) +
  geom_text(stat = "count", aes(label = after_stat(count)), 
            vjust = -0.5, size = 3) +
  scale_x_upset(order_by = "freq", n_intersections = 15) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray90", fill = NA, linewidth = 0.5)
  ) +
  labs(
    title = "Gene Database Intersections",
    subtitle = "UpSet plot showing overlap patterns between databases",
    x = "Database Combinations",
    y = "Number of Genes"
  )

ggsave(file.path(plot_dir, "gene_overlaps_upset_ggplot.png"), p4, 
       width = 10, height = 6, dpi = 300, bg = "white")

# 5. Create a comprehensive combined plot with all four analyses
# Layout: Gene counts by source (top left), Gene counts by technology (top right), 
#         Abundance distributions (bottom left), Gene intersections (bottom right)
combined_plot <- (p1 | p2) / (p3 | p4) + 
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Comprehensive Plasma Protein Quantification Analysis",
    subtitle = "Gene counts and technology classification (top), abundance distributions and intersections (bottom)",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  )

ggsave(file.path(plot_dir, "plasma_proteins_comprehensive.png"), combined_plot, 
       width = 20, height = 12, dpi = 300, bg = "white")

message("Plots saved to:", plot_dir)

# Create summary report
sink(file.path(output_dir, "analysis_summary.txt"))
cat("PLASMA PROTEIN ANALYSIS SUMMARY (WITH GENE DEDUPLICATION)\n")
cat("=========================================================\n\n")

cat("GENE DEDUPLICATION INFORMATION:\n")
cat("===============================\n")
cat("Multiple proteins mapping to the same gene have been deduplicated.\n")
cat("For genes with multiple protein entries, median quantification values are used.\n")
cat("This ensures accurate gene-level counts and quantification.\n\n")

cat("SUMMARY STATISTICS (AFTER DEDUPLICATION):\n")
cat("==========================================\n")
cat(sprintf("PeptideAtlas: %d unique genes\n", stats_summary$peptideatlas))
cat(sprintf("HPA MS: %d unique genes\n", stats_summary$hpa_ms))
cat(sprintf("HPA PEA: %d unique genes\n", stats_summary$hpa_pea))
cat(sprintf("HPA Immunoassay: %d unique genes\n", stats_summary$hpa_immunoassay))
cat(sprintf("GPMDB: %d unique genes\n", stats_summary$gpmdb))
cat(sprintf("PAXDB: %d unique genes\n", stats_summary$paxdb))
cat(sprintf("\nTotal unique genes across sources: %d\n", stats_summary$total_across_sources))
cat(sprintf("MS technologies total: %d unique genes\n", stats_summary$ms_technologies))
cat(sprintf("HPA total (all technologies): %d unique genes\n", stats_summary$hpa_total))
cat(rep("=", 60), "\n", sep = "")

cat("\nAnalysis completed successfully!\n")
cat("Generated files:\n")
cat("• Individual source plot: plasma_proteins_by_source.png\n")
cat("• Technology classification plot: plasma_proteins_by_technology.png\n")
cat("• Abundance distribution plot: plasma_proteins_abundance_distributions.png\n")
cat("• Gene overlap UpSet plot (classic): gene_overlaps_upset.png\n")
cat("• Gene overlap UpSet plot (ggplot2): gene_overlaps_upset_ggplot.png\n")
cat("• Comprehensive combined plot: plasma_proteins_comprehensive.png\n")
cat("• Data: outputs/plasma_protein_counts_summary.csv\n")
cat("• Report: outputs/analysis_summary.txt\n")
sink()

message("Analysis completed successfully!") 