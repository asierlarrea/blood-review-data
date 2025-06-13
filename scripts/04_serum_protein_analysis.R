#!/usr/bin/env Rscript
#
# Analysis of Serum Proteins from Multiple Data Sources
# Script: 04_serum_protein_analysis.R
# 
# Purpose: Analyze the total number of genes/proteins quantified in serum by different sources:
# - GPMDB (Ensembl Protein accessions)
# - PAXDB (Ensembl Protein accessions)  
# - HPA Immunoassay (Gene names)
#
# Updated to use consistent gene name mapping for accurate comparisons

# Load utilities and set up output directories
source("scripts/utilities/load_packages.R")
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "ggthemes", "patchwork", "UpSetR", "ggupset", "viridis")
load_packages(required_packages)

# Parse command line arguments (simple approach)
args <- commandArgs(trailingOnly = TRUE)
force_mapping <- "--force-mapping" %in% args

# Load gene mapping utility
source("scripts/data_processing/simple_id_mapping.R")

# Load gene deduplication utility
source("scripts/utilities/gene_deduplication.R")

# Set output directory
output_dir <- get_output_path("", subdir = "serum_protein")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read and process data from each source
message("Reading and processing serum data from each source...")

# 1. GPMDB Serum
message("Processing GPMDB serum data...")
gpmdb_serum_raw <- read_csv("data/raw/gpmdb/gpmdb_serum.csv", show_col_types = FALSE)
gpmdb_serum_raw$gene <- stringr::str_extract(gpmdb_serum_raw$description, "[A-Z0-9]+(?=,| |$)")

# Deduplicate genes using median quantification values
gpmdb_serum <- deduplicate_genes(gpmdb_serum_raw, "gene", "total", 
                               additional_cols = c("accession"), 
                               aggregation_method = "median")

# 2. PAXDB Serum
message("Processing PAXDB serum data...")
paxdb_serum_raw <- read_csv("data/raw/paxdb/paxdb_serum.csv", show_col_types = FALSE)
paxdb_serum_raw$ensp <- stringr::str_replace(paxdb_serum_raw$string_external_id, "^9606\\.", "")
paxdb_serum_raw$gene <- convert_to_gene_symbol(paxdb_serum_raw$ensp, force_mapping = force_mapping)

# Deduplicate genes using median quantification values
paxdb_serum <- deduplicate_genes(paxdb_serum_raw, "gene", "abundance", 
                               additional_cols = c("ensp"), 
                               aggregation_method = "median")

# 3. HPA Immunoassay Serum
message("Processing HPA Immunoassay serum data...")
hpa_serum_raw <- read_csv("data/raw/hpa/hpa_immunoassay_serum.csv", show_col_types = FALSE)
hpa_serum_raw <- hpa_serum_raw %>% rename(gene = Gene, expr = Concentration)

# Deduplicate genes using median quantification values
hpa_serum <- deduplicate_genes(hpa_serum_raw, "gene", "expr", 
                             additional_cols = c("ENSG_ID"), 
                             aggregation_method = "median")

# Create summary statistics
message("Creating summary statistics...")
stats_summary <- list(
  gpmdb_serum = length(unique(gpmdb_serum$gene)),
  paxdb_serum = length(unique(paxdb_serum$gene)),
  hpa_immunoassay_serum = length(unique(hpa_serum$gene))
)

# Calculate total unique genes across sources
all_genes <- unique(c(
  gpmdb_serum$gene,
  paxdb_serum$gene,
  hpa_serum$gene
))
stats_summary$total_across_sources <- length(all_genes)

# Calculate MS technologies total (GPMDB, PAXDB)
ms_genes <- unique(c(
  gpmdb_serum$gene,
  paxdb_serum$gene
))
stats_summary$ms_technologies <- length(ms_genes)

# Save summary statistics
write.csv(
  data.frame(
    Source = names(stats_summary),
    Gene_Count = unlist(stats_summary)
  ),
  file.path(output_dir, "serum_protein_counts_summary.csv"),
  row.names = FALSE
)

# Create visualizations
message("Creating visualizations...")

# Set up plot output directory
plot_dir <- get_output_path("04_serum_protein_analysis", subdir = "plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# 1. Bar plot by individual sources
source_data <- data.frame(
  Source = c("GPMDB", "PAXDB", "HPA Immunoassay"),
  Count = c(stats_summary$gpmdb_serum, stats_summary$paxdb_serum, 
            stats_summary$hpa_immunoassay_serum),
  Technology = c("MS", "MS", "Immunoassay"),
  Database = c("GPMDB", "PAXDB", "HPA")
)

p1 <- ggplot(source_data, aes(x = reorder(Source, Count), y = Count, fill = Technology)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("MS" = "#2E86AB", "Immunoassay" = "#F18F01")) +
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
    title = "Serum Proteins Quantified by Data Source",
    subtitle = "Number of unique genes detected in each database",
    x = "Data Source",
    y = "Number of Genes",
    fill = "Technology"
  )

ggsave(file.path(plot_dir, "serum_protein_counts_by_source.png"), 
       p1, width = 10, height = 6, dpi = 300, bg = "white")

# 2. Grouped bar plot by technology and total
summary_data <- data.frame(
  Category = c("GPMDB", "PAXDB", "HPA Immunoassay", "MS Technologies", "All Sources"),
  Count = c(stats_summary$gpmdb_serum, stats_summary$paxdb_serum, 
            stats_summary$hpa_immunoassay_serum, stats_summary$ms_technologies, 
            stats_summary$total_across_sources),
  Type = c("Individual", "Individual", "Individual", "Combined", "Combined")
)

p2 <- ggplot(summary_data, aes(x = reorder(Category, Count), y = Count, fill = Type)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("Individual" = "#4E79A7", "Combined" = "#E15759")) +
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
    title = "Serum Protein Detection: Individual vs Combined Sources",
    subtitle = "Comparison of gene counts across databases and technologies",
    x = "Database/Category",
    y = "Number of Genes",
    fill = "Category Type"
  )

ggsave(file.path(plot_dir, "serum_protein_counts_combined.png"), 
       p2, width = 10, height = 6, dpi = 300, bg = "white")

# 3. Venn diagram data preparation and UpSet plot
message("Creating overlap analysis...")

# Create gene lists for overlap analysis
gene_lists <- list(
  GPMDB = gpmdb_serum$gene[!is.na(gpmdb_serum$gene) & gpmdb_serum$gene != ""],
  PAXDB = paxdb_serum$gene[!is.na(paxdb_serum$gene) & paxdb_serum$gene != ""],
  HPA_Immunoassay = hpa_serum$gene[!is.na(hpa_serum$gene) & hpa_serum$gene != ""]
)

# Create UpSet plot
p3 <- UpSetR::upset(
  fromList(gene_lists),
  nsets = 3,
  sets = c("GPMDB", "PAXDB", "HPA_Immunoassay"),
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  text.scale = 1.2,
  point.size = 3,
  line.size = 1.2,
  mainbar.y.label = "Number of Genes",
  sets.x.label = "Total Genes per Database"
)

png(file.path(plot_dir, "serum_protein_overlap_upset.png"), 
    width = 12, height = 8, units = "in", res = 300, bg = "white")
print(p3)
dev.off()

# 4. Create quantification value distributions
message("Creating quantification distribution plots...")

# Prepare data for distribution plots
gpmdb_dist <- gpmdb_serum %>% 
  filter(!is.na(total) & total > 0) %>%
  mutate(Database = "GPMDB", log_value = log10(total))

paxdb_dist <- paxdb_serum %>% 
  filter(!is.na(abundance) & abundance > 0) %>%
  mutate(Database = "PAXDB", log_value = log10(abundance))

hpa_dist <- hpa_serum %>% 
  filter(!is.na(expr) & expr > 0) %>%
  mutate(Database = "HPA Immunoassay", log_value = log10(expr))

# Combined distribution data
dist_data <- bind_rows(
  gpmdb_dist %>% select(gene, Database, log_value),
  paxdb_dist %>% select(gene, Database, log_value),
  hpa_dist %>% select(gene, Database, log_value)
)

p4 <- ggplot(dist_data, aes(x = log_value, fill = Database)) +
  geom_histogram(alpha = 0.7, bins = 50) +
  facet_wrap(~Database, scales = "free", ncol = 1, 
             labeller = labeller(Database = c("GPMDB" = "(a) GPMDB", 
                                             "PAXDB" = "(b) PAXDB", 
                                             "HPA Immunoassay" = "(c) HPA Immunoassay"))) +
  scale_fill_viridis_d(option = "plasma") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "Distribution of Protein Quantification Values",
    subtitle = "Log10-transformed values across serum databases",
    x = "Log10(Quantification Value)",
    y = "Number of Proteins"
  )

ggsave(file.path(plot_dir, "serum_protein_quantification_distributions.png"), 
       p4, width = 10, height = 12, dpi = 300, bg = "white")

# 5. Density plots for direct comparison
p5 <- ggplot(dist_data, aes(x = log_value, fill = Database)) +
  geom_density(alpha = 0.6) +
  scale_fill_viridis_d(option = "plasma") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  ) +
  labs(
    title = "Protein Quantification Value Distributions Comparison",
    subtitle = "Density plots of log10-transformed values across serum databases",
    x = "Log10(Quantification Value)",
    y = "Density",
    fill = "Database"
  )

ggsave(file.path(plot_dir, "serum_protein_quantification_density.png"), 
       p5, width = 10, height = 6, dpi = 300, bg = "white")

# 6. Create overlap statistics table
message("Creating overlap statistics...")

# Calculate pairwise overlaps
overlap_stats <- data.frame(
  Comparison = character(),
  Overlap_Count = numeric(),
  Database1_Total = numeric(),
  Database2_Total = numeric(),
  Overlap_Percentage_DB1 = numeric(),
  Overlap_Percentage_DB2 = numeric()
)

# GPMDB vs PAXDB
gpmdb_paxdb_overlap <- length(intersect(gene_lists$GPMDB, gene_lists$PAXDB))
overlap_stats <- rbind(overlap_stats, data.frame(
  Comparison = "GPMDB vs PAXDB",
  Overlap_Count = gpmdb_paxdb_overlap,
  Database1_Total = length(gene_lists$GPMDB),
  Database2_Total = length(gene_lists$PAXDB),
  Overlap_Percentage_DB1 = round(gpmdb_paxdb_overlap / length(gene_lists$GPMDB) * 100, 2),
  Overlap_Percentage_DB2 = round(gpmdb_paxdb_overlap / length(gene_lists$PAXDB) * 100, 2)
))

# GPMDB vs HPA
gpmdb_hpa_overlap <- length(intersect(gene_lists$GPMDB, gene_lists$HPA_Immunoassay))
overlap_stats <- rbind(overlap_stats, data.frame(
  Comparison = "GPMDB vs HPA",
  Overlap_Count = gpmdb_hpa_overlap,
  Database1_Total = length(gene_lists$GPMDB),
  Database2_Total = length(gene_lists$HPA_Immunoassay),
  Overlap_Percentage_DB1 = round(gpmdb_hpa_overlap / length(gene_lists$GPMDB) * 100, 2),
  Overlap_Percentage_DB2 = round(gpmdb_hpa_overlap / length(gene_lists$HPA_Immunoassay) * 100, 2)
))

# PAXDB vs HPA
paxdb_hpa_overlap <- length(intersect(gene_lists$PAXDB, gene_lists$HPA_Immunoassay))
overlap_stats <- rbind(overlap_stats, data.frame(
  Comparison = "PAXDB vs HPA",
  Overlap_Count = paxdb_hpa_overlap,
  Database1_Total = length(gene_lists$PAXDB),
  Database2_Total = length(gene_lists$HPA_Immunoassay),
  Overlap_Percentage_DB1 = round(paxdb_hpa_overlap / length(gene_lists$PAXDB) * 100, 2),
  Overlap_Percentage_DB2 = round(paxdb_hpa_overlap / length(gene_lists$HPA_Immunoassay) * 100, 2)
))

# Three-way overlap
three_way_overlap <- length(Reduce(intersect, gene_lists))
overlap_stats <- rbind(overlap_stats, data.frame(
  Comparison = "All Three Databases",
  Overlap_Count = three_way_overlap,
  Database1_Total = NA,
  Database2_Total = NA,
  Overlap_Percentage_DB1 = NA,
  Overlap_Percentage_DB2 = NA
))

write.csv(overlap_stats, file.path(output_dir, "serum_protein_overlap_statistics.csv"), row.names = FALSE)

# 7. Create a comprehensive summary plot
p6 <- p1 + p2 + 
  plot_layout(ncol = 1) +
  plot_annotation(
    title = "Comprehensive Serum Protein Analysis",
    subtitle = "Comparison across GPMDB, PAXDB, and HPA databases",
    tag_levels = list(c('(a)', '(b)')),
    theme = theme(plot.title = element_text(size = 18, face = "bold"))
  )

ggsave(file.path(plot_dir, "serum_protein_comprehensive_summary.png"), 
       p6, width = 12, height = 12, dpi = 300, bg = "white")

# 8. Create an extended comprehensive plot with distributions
p7 <- p1 + p2 + p5 + 
  plot_layout(ncol = 1) +
  plot_annotation(
    title = "Comprehensive Serum Protein Analysis with Quantification Distributions",
    subtitle = "Complete comparison across GPMDB, PAXDB, and HPA databases",
    tag_levels = list(c('(a)', '(b)', '(c)')),
    theme = theme(plot.title = element_text(size = 18, face = "bold"))
  )

ggsave(file.path(plot_dir, "serum_protein_extended_comprehensive.png"), 
       p7, width = 12, height = 18, dpi = 300, bg = "white")

# Print summary information
message("\n=== SERUM PROTEIN ANALYSIS SUMMARY ===")
message(sprintf("GPMDB serum proteins: %s", scales::comma(stats_summary$gpmdb_serum)))
message(sprintf("PAXDB serum proteins: %s", scales::comma(stats_summary$paxdb_serum)))
message(sprintf("HPA Immunoassay serum proteins: %s", scales::comma(stats_summary$hpa_immunoassay_serum)))
message(sprintf("MS technologies total: %s", scales::comma(stats_summary$ms_technologies)))
message(sprintf("Total across all sources: %s", scales::comma(stats_summary$total_across_sources)))
message(sprintf("Three-way overlap: %s", scales::comma(three_way_overlap)))

message("\n=== OUTPUT FILES ===")
message(sprintf("Summary statistics: %s", file.path(output_dir, "serum_protein_counts_summary.csv")))
message(sprintf("Overlap statistics: %s", file.path(output_dir, "serum_protein_overlap_statistics.csv")))
message(sprintf("Plots directory: %s", plot_dir))

message("\n=== ANALYSIS COMPLETE ===") 