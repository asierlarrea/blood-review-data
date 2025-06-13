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
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "ggthemes", "patchwork")
load_packages(required_packages)

# Parse command line arguments (simple approach)
args <- commandArgs(trailingOnly = TRUE)
force_mapping <- "--force-mapping" %in% args

# Load gene mapping utility
source("scripts/data_processing/simple_id_mapping.R")

# Set output directory
output_dir <- get_output_path("", subdir = "plasma_protein")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read and process data from each source
message("Reading and processing data from each source...")

# 1. PeptideAtlas
message("Processing PeptideAtlas data...")
peptideatlas <- read_csv("data/raw/peptideatlas/peptideatlas.csv", show_col_types = FALSE)
peptideatlas$gene <- convert_to_gene_symbol(peptideatlas$biosequence_accession, force_mapping = force_mapping)
peptideatlas <- peptideatlas %>% filter(!is.na(norm_PSMs_per_100K))

# 2. HPA MS
message("Processing HPA MS data...")
hpa_ms <- read_csv("data/raw/hpa/hpa_ms.csv", show_col_types = FALSE, skip = 1)
hpa_ms <- hpa_ms %>% rename(gene = Gene, expr = Concentration)

# 3. HPA PEA
message("Processing HPA PEA data...")
hpa_pea <- read_csv("data/raw/hpa/hpa_pea.csv", show_col_types = FALSE)
hpa_pea <- hpa_pea %>% rename(gene = Gene, expr = `Variation between individuals`)

# 4. HPA Immunoassay
message("Processing HPA Immunoassay data...")
hpa_imm <- read_csv("data/raw/hpa/hpa_immunoassay_plasma.csv", show_col_types = FALSE)
hpa_imm <- hpa_imm %>% rename(gene = Gene, expr = Concentration)

# 5. GPMDB
message("Processing GPMDB data...")
gpmdb <- read_csv("data/raw/gpmdb/gpmdb_plasma.csv", show_col_types = FALSE)
gpmdb$gene <- stringr::str_extract(gpmdb$description, "[A-Z0-9]+(?=,| |$)")

# 6. PAXDB
message("Processing PAXDB data...")
paxdb <- read_csv("data/raw/paxdb/paxdb_plasma.csv", show_col_types = FALSE)
paxdb$ensp <- stringr::str_replace(paxdb$string_external_id, "^9606\\.", "")
paxdb$gene <- convert_to_gene_symbol(paxdb$ensp, force_mapping = force_mapping)

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

# Calculate MS technologies total
ms_genes <- unique(c(
  peptideatlas$gene,
  hpa_ms$gene,
  gpmdb$gene
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
  Technology = c("MS", "MS", "PEA", "Immunoassay", "MS", "Expression"),
  Database = c("PeptideAtlas", "HPA", "HPA", "HPA", "GPMDB", "PAXDB")
)

p1 <- ggplot(source_data, aes(x = reorder(Source, Count), y = Count, fill = Technology)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("MS" = "#2E86AB", "PEA" = "#A23B72", 
                               "Immunoassay" = "#F18F01", "Expression" = "#C73E1D")) +
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

# 2. Bar plot by technology grouping
tech_data <- data.frame(
  Technology = c("Mass Spectrometry", "HPA (All Technologies)", "Total Across Sources"),
  Count = c(stats_summary$ms_technologies, stats_summary$hpa_total, stats_summary$total_across_sources),
  Type = c("Technology", "Database", "Overall")
)

p2 <- ggplot(tech_data, aes(x = reorder(Technology, Count), y = Count, fill = Type)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("Technology" = "#2E86AB", "Database" = "#A23B72", "Overall" = "#C73E1D")) +
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
    title = "Plasma Proteins by Technology and Overall Coverage",
    subtitle = "Aggregated counts across different analytical approaches",
    x = "Category",
    y = "Number of Genes",
    fill = "Category Type"
  )

ggsave(file.path(plot_dir, "plasma_proteins_by_technology.png"), p2, 
       width = 10, height = 6, dpi = 300, bg = "white")

# 3. Main databases comparison (excluding HPA subcategories)
main_data <- data.frame(
  Database = c("PeptideAtlas", "HPA (Combined)", "GPMDB", "PAXDB"),
  Count = c(stats_summary$peptideatlas, stats_summary$hpa_total, 
            stats_summary$gpmdb, stats_summary$paxdb),
  Description = c("MS-based proteomics", "Multi-technology platform", 
                  "MS-based proteomics", "Expression database")
)

p3 <- ggplot(main_data, aes(x = reorder(Database, Count), y = Count, fill = Description)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("MS-based proteomics" = "#2E86AB", 
                               "Multi-technology platform" = "#A23B72",
                               "Expression database" = "#F18F01")) +
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
    title = "Main Plasma Protein Databases Comparison",
    subtitle = "Unique genes detected in major databases and platforms",
    x = "Database/Platform",
    y = "Number of Genes",
    fill = "Database Type"
  )

ggsave(file.path(plot_dir, "plasma_proteins_main_databases.png"), p3, 
       width = 10, height = 6, dpi = 300, bg = "white")

# 4. Create a comprehensive combined plot
combined_plot <- (p1 / p2) | p3
combined_plot <- combined_plot + 
  plot_annotation(
    title = "Comprehensive Plasma Protein Quantification Analysis",
    subtitle = "Comparison across different data sources, technologies, and databases",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  )

ggsave(file.path(plot_dir, "plasma_proteins_comprehensive.png"), combined_plot, 
       width = 16, height = 12, dpi = 300, bg = "white")

message("Plots saved to:", plot_dir)

# Create summary report
sink(file.path(output_dir, "analysis_summary.txt"))
cat("PLASMA PROTEIN ANALYSIS SUMMARY\n")
cat("==============================\n\n")

cat("SUMMARY STATISTICS:\n")
cat("==================\n")
cat(sprintf("PeptideAtlas: %d genes\n", stats_summary$peptideatlas))
cat(sprintf("HPA MS: %d genes\n", stats_summary$hpa_ms))
cat(sprintf("HPA PEA: %d genes\n", stats_summary$hpa_pea))
cat(sprintf("HPA Immunoassay: %d genes\n", stats_summary$hpa_immunoassay))
cat(sprintf("GPMDB: %d genes\n", stats_summary$gpmdb))
cat(sprintf("PAXDB: %d genes\n", stats_summary$paxdb))
cat(sprintf("\nTotal genes across sources: %d\n", stats_summary$total_across_sources))
cat(sprintf("MS technologies total: %d genes\n", stats_summary$ms_technologies))
cat(sprintf("HPA total (all technologies): %d genes\n", stats_summary$hpa_total))
cat(rep("=", 60), "\n", sep = "")

cat("\nAnalysis completed successfully!\n")
cat("Generated files:\n")
cat("• Plots: plots/plasma_proteins_*.png\n")
cat("• Data: outputs/plasma_protein_counts_summary.csv\n")
cat("• Report: outputs/analysis_summary.txt\n")
sink()

message("Analysis completed successfully!") 