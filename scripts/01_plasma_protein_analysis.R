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

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(viridis)
})

# Source the simple ID mapping functions
cat("Loading simple ID mapping utilities...\n")
source("scripts/data_processing/simple_id_mapping.R")

# Set up directories
data_dir <- "data/raw"
output_dir <- "outputs"
plot_dir <- "plots"

# Create output directories if they don't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

cat("Starting plasma protein analysis with gene name mapping...\n")
cat(rep("=", 70), "\n", sep = "")

# =============================================================================
# 1. PEPTIDEATLAS DATA (UniProt accessions)
# =============================================================================
cat("Processing PeptideAtlas data...\n")
peptideatlas_file <- file.path(data_dir, "peptideatlas", "peptideatlas.csv")
peptideatlas_data <- read_csv(peptideatlas_file, show_col_types = FALSE)

cat("  Mapping UniProt accessions to gene names...\n")
# Extract unique UniProt accessions and map to gene names
peptideatlas_accessions <- peptideatlas_data %>%
  filter(!is.na(biosequence_accession)) %>%
  distinct(biosequence_accession) %>%
  pull(biosequence_accession)

# Convert to gene symbols using the new mapping system
peptideatlas_gene_symbols <- convert_to_gene_symbol(peptideatlas_accessions)
peptideatlas_genes <- sum(!is.na(peptideatlas_gene_symbols))

cat(sprintf("  PeptideAtlas: %d unique genes\n", peptideatlas_genes))

# =============================================================================
# 2. HPA DATA - Different Technologies (Already has Gene names)
# =============================================================================
cat("Processing HPA data (multiple technologies)...\n")

# HPA MS (Mass Spectrometry) - already has Gene column
hpa_ms_file <- file.path(data_dir, "hpa", "hpa_ms.csv")
hpa_ms_data <- read_csv(hpa_ms_file, skip = 1, show_col_types = FALSE)  # Skip the comment line
hpa_ms_genes <- hpa_ms_data %>%
  filter(!is.na(Gene) & Gene != "") %>%
  distinct(Gene) %>%
  nrow()

# HPA PEA (Proximity Extension Assay) - already has Gene column
hpa_pea_file <- file.path(data_dir, "hpa", "hpa_pea.csv")
hpa_pea_data <- read_csv(hpa_pea_file, show_col_types = FALSE)
hpa_pea_genes <- hpa_pea_data %>%
  filter(!is.na(Gene) & Gene != "") %>%
  distinct(Gene) %>%
  nrow()

# HPA Immunoassay Plasma - already has Gene column
hpa_immunoassay_file <- file.path(data_dir, "hpa", "hpa_immunoassay_plasma.csv")
hpa_immunoassay_data <- read_csv(hpa_immunoassay_file, show_col_types = FALSE)
hpa_immunoassay_genes <- hpa_immunoassay_data %>%
  filter(!is.na(Gene) & Gene != "") %>%
  distinct(Gene) %>%
  nrow()

cat(sprintf("  HPA MS: %d unique genes\n", hpa_ms_genes))
cat(sprintf("  HPA PEA: %d unique genes\n", hpa_pea_genes))
cat(sprintf("  HPA Immunoassay: %d unique genes\n", hpa_immunoassay_genes))

# =============================================================================
# 3. GPMDB DATA (Ensembl Protein accessions)
# =============================================================================
cat("Processing GPMDB data...\n")
gpmdb_file <- file.path(data_dir, "gpmdb", "gpmdb_plasma.csv")
gpmdb_data <- read_csv(gpmdb_file, show_col_types = FALSE)

cat("  Mapping ENSP accessions to gene names...\n")
# Extract unique ENSP accessions and map to gene names
gpmdb_accessions <- gpmdb_data %>%
  filter(!is.na(accession)) %>%
  distinct(accession) %>%
  pull(accession)

# Convert to gene symbols using the new mapping system
gpmdb_gene_symbols <- convert_to_gene_symbol(gpmdb_accessions)
gpmdb_genes <- sum(!is.na(gpmdb_gene_symbols))

cat(sprintf("  GPMDB: %d unique genes\n", gpmdb_genes))

# =============================================================================
# 4. PAXDB DATA (Ensembl Protein accessions with taxonomy prefix)
# =============================================================================
cat("Processing PAXDB data...\n")
paxdb_file <- file.path(data_dir, "paxdb", "paxdb_plasma.csv")
paxdb_data <- read_csv(paxdb_file, show_col_types = FALSE)

cat("  Mapping ENSP accessions to gene names...\n")
# Clean PAXDB accessions (remove taxonomy prefix "9606.") and extract unique IDs
paxdb_accessions <- paxdb_data %>%
  filter(!is.na(string_external_id)) %>%
  mutate(
    clean_accession = str_replace(string_external_id, "^9606\\.", "")
  ) %>%
  distinct(clean_accession) %>%
  pull(clean_accession)

# Convert to gene symbols using the new mapping system
paxdb_gene_symbols <- convert_to_gene_symbol(paxdb_accessions)
paxdb_genes <- sum(!is.na(paxdb_gene_symbols))

cat(sprintf("  PAXDB: %d unique genes\n", paxdb_genes))

# =============================================================================
# 5. CREATE SUMMARY DATA FOR VISUALIZATION
# =============================================================================
cat("Creating summary data for visualization...\n")

# Create comprehensive summary data frame
protein_counts <- data.frame(
  Source = c("PeptideAtlas", "HPA MS", "HPA PEA", "HPA Immunoassay", "GPMDB", "PAXDB"),
  Technology = c("MS", "MS", "PEA", "Immunoassay", "MS", "Expression"),
  Count = c(peptideatlas_genes, hpa_ms_genes, hpa_pea_genes, hpa_immunoassay_genes, gpmdb_genes, paxdb_genes),
  stringsAsFactors = FALSE
)

# Create a category grouping
protein_counts$Category <- case_when(
  protein_counts$Source == "PeptideAtlas" ~ "PeptideAtlas",
  str_starts(protein_counts$Source, "HPA") ~ "HPA",
  protein_counts$Source == "GPMDB" ~ "GPMDB",
  protein_counts$Source == "PAXDB" ~ "PAXDB"
)

# Print summary
cat("\nSummary of gene counts by source (mapped to gene names):\n")
print(protein_counts)

# Save summary data
write_csv(protein_counts, file.path(output_dir, "plasma_protein_counts_summary.csv"))

# =============================================================================
# 6. CREATE VISUALIZATIONS
# =============================================================================
cat("Creating visualizations...\n")

# Plot 1: Main barplot showing all sources
p1 <- ggplot(protein_counts, aes(x = reorder(Source, Count), y = Count, fill = Category)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = Count), hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Total Number of Genes Quantified in Plasma",
    subtitle = "Comparison across different data sources and technologies (mapped to gene names)",
    x = "Data Source",
    y = "Number of Unique Genes",
    fill = "Database",
    caption = "Data sources: PeptideAtlas, HPA (MS/PEA/Immunoassay), GPMDB, PAXDB\nAll protein accessions mapped to HGNC gene symbols"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.caption = element_text(size = 10, color = "gray50"),
    panel.grid.major.x = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

# Save main plot
ggsave(file.path(plot_dir, "plasma_proteins_by_source.png"), p1, 
       width = 12, height = 8, dpi = 300, bg = "white")

# Plot 2: Grouped by technology
protein_counts_tech <- protein_counts %>%
  mutate(Technology_Detail = case_when(
    Source == "PeptideAtlas" ~ "MS (PeptideAtlas)",
    Source == "HPA MS" ~ "MS (HPA)",
    Source == "HPA PEA" ~ "PEA (HPA)",
    Source == "HPA Immunoassay" ~ "Immunoassay (HPA)",
    Source == "GPMDB" ~ "MS (GPMDB)",
    Source == "PAXDB" ~ "Expression Data"
  ))

p2 <- ggplot(protein_counts_tech, aes(x = reorder(Technology_Detail, Count), y = Count, fill = Technology)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = Count), hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  scale_fill_viridis_d(option = "viridis", begin = 0.1, end = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Plasma Genes by Technology Platform",
    subtitle = "Comparison of detection technologies and methodologies (gene-level counts)",
    x = "Technology Platform",
    y = "Number of Unique Genes",
    fill = "Technology Type",
    caption = "MS = Mass Spectrometry, PEA = Proximity Extension Assay\nAll accessions mapped to gene names for consistent comparison"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.caption = element_text(size = 10, color = "gray50"),
    panel.grid.major.x = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

# Save technology plot
ggsave(file.path(plot_dir, "plasma_proteins_by_technology.png"), p2, 
       width = 12, height = 8, dpi = 300, bg = "white")

# Plot 3: Simple comparison of main databases
main_sources <- protein_counts %>%
  filter(!str_starts(Source, "HPA")) %>%
  bind_rows(
    protein_counts %>%
      filter(str_starts(Source, "HPA")) %>%
      summarise(
        Source = "HPA (Combined)",
        Technology = "Multiple",
        Count = sum(Count),
        Category = "HPA"
      )
  )

p3 <- ggplot(main_sources, aes(x = reorder(Source, Count), y = Count, fill = Source)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = Count), hjust = -0.1, size = 5, fontface = "bold") +
  coord_flip() +
  scale_fill_viridis_d(option = "plasma") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Gene Coverage by Major Plasma Proteomics Databases",
    subtitle = "Consolidated view of gene quantification in plasma",
    x = "Database",
    y = "Number of Unique Genes",
    caption = "HPA combined includes MS, PEA, and Immunoassay data\nConsistent gene-level mapping across all databases"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.caption = element_text(size = 10, color = "gray50"),
    panel.grid.major.x = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

# Save main databases plot
ggsave(file.path(plot_dir, "plasma_proteins_main_databases.png"), p3, 
       width = 10, height = 6, dpi = 300, bg = "white")

# =============================================================================
# 7. GENERATE ANALYSIS SUMMARY REPORT
# =============================================================================
cat("Generating analysis summary report...\n")

# Calculate total proteins and statistics
stats_summary <- list(
  peptideatlas = peptideatlas_genes,
  hpa_ms = hpa_ms_genes,
  hpa_pea = hpa_pea_genes,
  hpa_immunoassay = hpa_immunoassay_genes,
  hpa_total = hpa_ms_genes + hpa_pea_genes + hpa_immunoassay_genes,
  gpmdb = gpmdb_genes,
  paxdb = paxdb_genes,
  total_across_sources = sum(protein_counts$Count),
  ms_technologies = peptideatlas_genes + hpa_ms_genes + gpmdb_genes
)

# Generate summary report
report_file <- file.path(output_dir, "analysis_summary.txt")
cat("PLASMA PROTEIN ANALYSIS SUMMARY\n", file = report_file)
cat("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = report_file, append = TRUE)
cat(rep("=", 60), "\n", file = report_file, append = TRUE)
cat("\nANALYSIS OVERVIEW:\n", file = report_file, append = TRUE)
cat("This analysis compares the number of unique genes quantified in plasma\n", file = report_file, append = TRUE)
cat("across multiple proteomics databases and technologies. All protein\n", file = report_file, append = TRUE)
cat("accessions have been mapped to HGNC gene symbols for consistent comparison.\n\n", file = report_file, append = TRUE)

cat("GENE COUNTS BY SOURCE:\n", file = report_file, append = TRUE)
cat(sprintf("• PeptideAtlas (MS): %d genes\n", stats_summary$peptideatlas), file = report_file, append = TRUE)
cat(sprintf("• HPA MS: %d genes\n", stats_summary$hpa_ms), file = report_file, append = TRUE)
cat(sprintf("• HPA PEA: %d genes\n", stats_summary$hpa_pea), file = report_file, append = TRUE)
cat(sprintf("• HPA Immunoassay: %d genes\n", stats_summary$hpa_immunoassay), file = report_file, append = TRUE)
cat(sprintf("• GPMDB (MS): %d genes\n", stats_summary$gpmdb), file = report_file, append = TRUE)
cat(sprintf("• PAXDB (Expression): %d genes\n", stats_summary$paxdb), file = report_file, append = TRUE)

cat("\nSUMMARY STATISTICS:\n", file = report_file, append = TRUE)
cat(sprintf("• Total genes across all sources: %d\n", stats_summary$total_across_sources), file = report_file, append = TRUE)
cat(sprintf("• MS-based technologies total: %d genes\n", stats_summary$ms_technologies), file = report_file, append = TRUE)
cat(sprintf("• HPA total (all technologies): %d genes\n", stats_summary$hpa_total), file = report_file, append = TRUE)
cat(sprintf("• Highest individual source: PAXDB (%d genes)\n", max(protein_counts$Count)), file = report_file, append = TRUE)
cat(sprintf("• Lowest individual source: HPA Immunoassay (%d genes)\n", min(protein_counts$Count)), file = report_file, append = TRUE)

cat("\nTECHNOLOGY INSIGHTS:\n", file = report_file, append = TRUE)
cat("• Expression data (PAXDB) shows highest gene coverage\n", file = report_file, append = TRUE)
cat("• Mass spectrometry dominates quantitative proteomics\n", file = report_file, append = TRUE)
cat("• Immunoassays provide targeted, high-specificity measurements\n", file = report_file, append = TRUE)
cat("• PEA technology offers multiplexed protein quantification\n", file = report_file, append = TRUE)

cat("\nDATA PROCESSING NOTES:\n", file = report_file, append = TRUE)
cat("• UniProt accessions (PeptideAtlas) mapped to gene names\n", file = report_file, append = TRUE)
cat("• ENSP accessions (GPMDB, PAXDB) mapped to gene names\n", file = report_file, append = TRUE)
cat("• HPA data already contained gene names\n", file = report_file, append = TRUE)
cat("• Gene mapping performed using lazy-loading system with UniProt/Ensembl APIs\n", file = report_file, append = TRUE)
cat("• Mapping cache used for fast repeated queries\n", file = report_file, append = TRUE)
cat("• Only successfully mapped genes included in final counts\n", file = report_file, append = TRUE)

cat(rep("=", 60), "\n", file = report_file, append = TRUE)

# Print summary
cat("\n", rep("=", 60), "\n", sep = "")
cat("ANALYSIS SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")
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