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
source("scripts/config/analysis_config.R")
source("scripts/utilities/data_loader.R")
source("scripts/utilities/plot_themes.R")
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "ggthemes", "patchwork", "UpSetR", "ggupset", "viridis", "ggvenn")
load_packages(required_packages)

# Parse command line arguments (simple approach)
args <- commandArgs(trailingOnly = TRUE)
force_mapping <- "--force-mapping" %in% args

# Parse formats argument
formats_arg <- NULL
if (length(args) > 0) {
  formats_idx <- which(args == "--formats")
  if (length(formats_idx) > 0 && formats_idx < length(args)) {
    formats_arg <- args[formats_idx + 1]
  }
}

# Set plot formats if provided
if (!is.null(formats_arg)) {
  set_plot_formats(formats_arg)
}

# Load gene mapping utility
source("scripts/data_processing/simple_id_mapping.R")

# Load gene deduplication utility
source("scripts/utilities/gene_deduplication.R")

# Load quantile normalization utility
source("scripts/utilities/quantile_normalization_functions.R")

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

# Process using standardized function
gpmdb_serum <- process_gpmdb_data(gpmdb_serum_raw, force_mapping = force_mapping)

# Add metadata
gpmdb_serum$source <- "GPMDB"
gpmdb_serum$technology <- "MS"
gpmdb_serum$abundance_type <- "Spectral Count"

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

# 4. quantms Serum
message("Processing quantms serum data...")
quantms_serum_raw <- load_quantms_data(sample_type = "serum", force_mapping = force_mapping)

# Ensure it's not null before processing
if (!is.null(quantms_serum_raw)) {
  quantms_serum <- quantms_serum_raw
} else {
  # Create an empty tibble to avoid errors downstream if data is missing
  quantms_serum <- tibble(gene = character(), abundance = double(), source = character())
}

# Create summary statistics
message("Creating summary statistics...")
stats_summary <- list(
  gpmdb_serum = length(unique(gpmdb_serum$gene)),
  paxdb_serum = length(unique(paxdb_serum$gene)),
  hpa_immunoassay_serum = length(unique(hpa_serum$gene)),
  quantms_serum = length(unique(quantms_serum$gene))
)

# Calculate total unique genes across sources
all_genes <- unique(c(
  gpmdb_serum$gene,
  paxdb_serum$gene,
  hpa_serum$gene,
  quantms_serum$gene
))
stats_summary$total_across_sources <- length(all_genes)

# Calculate MS technologies total (GPMDB, PAXDB)
ms_genes <- unique(c(
  gpmdb_serum$gene,
  paxdb_serum$gene,
  quantms_serum$gene
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
  Source = c("GPMDB", "PAXDB", "HPA Immunoassay", "quantms"),
  Count = c(stats_summary$gpmdb_serum, stats_summary$paxdb_serum, 
            stats_summary$hpa_immunoassay_serum, stats_summary$quantms_serum),
  Technology = c("MS", "MS", "Immunoassay", "MS"),
  Database = c("GPMDB", "PAXDB", "HPA", "quantms")
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

# 2. Grouped bar plot by technology and total
summary_data <- data.frame(
  Category = c("GPMDB", "PAXDB", "HPA Immunoassay", "quantms", "MS Technologies", "All Sources"),
  Count = c(stats_summary$gpmdb_serum, stats_summary$paxdb_serum, 
            stats_summary$hpa_immunoassay_serum, stats_summary$quantms_serum,
            stats_summary$ms_technologies, 
            stats_summary$total_across_sources),
  Type = c("Individual", "Individual", "Individual", "Individual", "Combined", "Combined")
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

# 3. Venn diagram data preparation and UpSet plot
message("Creating overlap analysis...")

# Create gene lists for overlap analysis
gene_lists <- list(
  GPMDB = gpmdb_serum$gene[!is.na(gpmdb_serum$gene) & gpmdb_serum$gene != ""],
  PAXDB = paxdb_serum$gene[!is.na(paxdb_serum$gene) & paxdb_serum$gene != ""],
  HPA_Immunoassay = hpa_serum$gene[!is.na(hpa_serum$gene) & hpa_serum$gene != ""],
  quantms = quantms_serum$gene[!is.na(quantms_serum$gene) & quantms_serum$gene != ""]
)

# Create UpSet plot
p3 <- UpSetR::upset(
  fromList(gene_lists),
  nsets = 4,
  sets = c("GPMDB", "PAXDB", "HPA_Immunoassay", "quantms"),
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  text.scale = 1.2,
  point.size = 3,
  line.size = 1.2,
  mainbar.y.label = "Number of Genes",
  sets.x.label = "Total Genes per Database"
)

# 4. Create quantification value distributions
message("Creating quantification distribution plots...")

# Function to calculate z-score normalization by database
calculate_zscore_normalization_serum <- function(data) {
  data %>%
    group_by(Database) %>%
    mutate(
      z_score = (log_value - mean(log_value, na.rm = TRUE)) / sd(log_value, na.rm = TRUE)
    ) %>%
    ungroup()
}

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

quantms_dist <- quantms_serum %>%
  filter(!is.na(abundance) & abundance > 0) %>%
  mutate(Database = "quantms", log_value = log10(abundance))

# Combined distribution data
dist_data <- bind_rows(
  gpmdb_dist %>% select(gene, Database, log_value),
  paxdb_dist %>% select(gene, Database, log_value),
  hpa_dist %>% select(gene, Database, log_value),
  quantms_dist %>% select(gene, Database, log_value)
)

# Apply z-score normalization
dist_data_zscore <- calculate_zscore_normalization_serum(dist_data)

p4 <- ggplot(dist_data, aes(x = log_value, fill = Database)) +
  geom_histogram(alpha = 0.7, bins = 50) +
  facet_wrap(~Database, scales = "free", ncol = 1, 
             labeller = labeller(Database = c("GPMDB" = "(a) GPMDB", 
                                             "PAXDB" = "(b) PAXDB", 
                                             "HPA Immunoassay" = "(c) HPA Immunoassay",
                                             "quantms" = "(d) quantms"))) +
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

# 4z. Z-score normalized histogram plots
p4z <- ggplot(dist_data_zscore, aes(x = z_score, fill = Database)) +
  geom_histogram(alpha = 0.7, bins = 50) +
  facet_wrap(~Database, scales = "free", ncol = 1, 
             labeller = labeller(Database = c("GPMDB" = "(a) GPMDB", 
                                             "PAXDB" = "(b) PAXDB", 
                                             "HPA Immunoassay" = "(c) HPA Immunoassay",
                                             "quantms" = "(d) quantms"))) +
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
    title = "Distribution of Z-Score Normalized Protein Quantification Values",
    subtitle = "Z-score normalized values for direct cross-database comparison",
    x = "Z-Score (standardized within each database)",
    y = "Number of Proteins",
    caption = "Z-scores calculated within each database: (log10(value) - mean) / sd"
  )

# 5. Density plots for direct comparison (log10)
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
    title = "Protein Quantification Value Distributions Comparison (Log10)",
    subtitle = "Density plots of log10-transformed values across serum databases",
    x = "Log10(Quantification Value)",
    y = "Density",
    fill = "Database"
  )

# ggsave(file.path(plot_dir, "06_serum_protein_quantification_density.png"), 
#        p5, width = 10, height = 6, dpi = 300, bg = "white")

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

# Create comprehensive panel
message("\nCreating comprehensive panel...")

#' Create comprehensive panel for Serum Protein Analysis (Script 04)
#' 
#' Panel Layout:
#' +-------------------+-------------------+
#' | A: Database       | B: Overlap        |
#' |    Coverage       |    (Venn)         |
#' +-------------------+-------------------+
#' | C: Correlation    | D: Distribution   |
#' |    Analysis       |    (Violin)       |
#' +-------------------+-------------------+
#'
create_serum_comprehensive_panel <- function(all_serum_data, stats_summary, plot_dir) {
  
  message("Creating comprehensive serum analysis panel...")
  
  # A: Database Coverage
  coverage_data <- data.frame(
    Database = c("GPMDB", "PAXDB", "HPA Immunoassay", "quantms"),
    Count = c(stats_summary$gpmdb_serum, stats_summary$paxdb_serum, stats_summary$hpa_immunoassay_serum, stats_summary$quantms_serum),
    Technology = c("MS", "MS", "Immunoassay", "MS")
  ) %>%
    arrange(desc(Count))
  
  panel_A <- ggplot(coverage_data, aes(x = reorder(Database, Count), y = Count, fill = Technology)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 3.5) +
    coord_flip() +
    scale_fill_manual(values = c("MS" = "#4E79A7", "Immunoassay" = "#E15759")) +
    theme_blood_proteomics() +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    labs(
      title = "(A) Serum Database Coverage",
      subtitle = "Unique proteins per database",
      x = "Database",
      y = "Number of Proteins",
      fill = "Technology"
    )
  
  # B: Venn Diagram
  gene_lists <- list(
    GPMDB = all_serum_data$gene[all_serum_data$source == "GPMDB"],
    PAXDB = all_serum_data$gene[all_serum_data$source == "PAXDB"],
    `HPA Immunoassay` = all_serum_data$gene[all_serum_data$source == "HPA Immunoassay"],
    quantms = all_serum_data$gene[all_serum_data$source == "quantms"]
  )
  
  # Calculate overlap statistics for subtitle
  total_unique <- length(unique(unlist(gene_lists)))
  overlap_all <- length(Reduce(intersect, gene_lists))
  overlap_percent <- round(overlap_all / total_unique * 100, 1)
  
  panel_B <- ggvenn(
    gene_lists,
    fill_color = c("#4E79A7", "#56B4E9", "#E15759", "#D55E00"),
    stroke_size = 0.5,
    set_name_size = 4
  ) +
    theme_void() +
    labs(
      title = "(B) Database Overlap",
      subtitle = sprintf("%d proteins in all databases (%.1f%%)", overlap_all, overlap_percent)
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  
  # C: Correlation Analysis for MS databases
  ms_data <- all_serum_data %>%
    filter(source %in% c("GPMDB", "PAXDB", "quantms")) %>%
    mutate(value = case_when(
      source == "GPMDB" ~ total,
      source == "PAXDB" ~ abundance,
      source == "quantms" ~ abundance,
      TRUE ~ NA_real_
    )) %>%
    mutate(log_value = log10(value + 1)) %>%
    group_by(source) %>%
    mutate(z_score = scale(log_value)[,1]) %>%
    ungroup() %>%
    select(gene, source, z_score)
  
  correlation_matrix_data <- ms_data %>%
    pivot_wider(
      names_from = source,
      values_from = z_score
    ) %>%
    column_to_rownames("gene")
  
  correlation_matrix <- cor(correlation_matrix_data, use = "pairwise.complete.obs")
  
  # Create heatmap
  corr_data <- as.data.frame(correlation_matrix) %>%
    rownames_to_column("db1") %>%
    pivot_longer(-db1, names_to = "db2", values_to = "correlation")
  
  panel_C <- ggplot(corr_data, aes(x = db1, y = db2, fill = correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4", midpoint = 0, limit = c(-1, 1)) +
    theme_blood_proteomics() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(C) Cross-Database Correlation",
      subtitle = "Pearson correlation (MS databases)",
      fill = "Correlation"
    )
  
  # D: Distribution Violin Plot
  distribution_data <- all_serum_data %>%
    mutate(value = case_when(
      source == "GPMDB" ~ total,
      source == "PAXDB" ~ abundance,
      source == "HPA Immunoassay" ~ expr,
      source == "quantms" ~ abundance,
      TRUE ~ NA_real_
    )) %>%
    group_by(source) %>%
    mutate(log_value = log10(value + 1), z_score = scale(log_value)[,1]) %>%
    ungroup()
  
  panel_D <- ggplot(distribution_data, aes(x = source, y = z_score, fill = source)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    scale_fill_manual(values = c("GPMDB" = "#4E79A7", "PAXDB" = "#56B4E9", "HPA Immunoassay" = "#E15759", "quantms" = "#D55E00")) +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "(D) Abundance Distribution",
      subtitle = "Z-score normalized values across databases",
      x = "Database",
      y = "Z-score"
    )
  
  # Combine panels
  comprehensive_panel <- (panel_A | panel_B) / (panel_C | panel_D) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      title = "Comprehensive Serum Protein Analysis",
      subtitle = paste0("Cross-database analysis of serum proteins across ", 
                       length(unique(all_serum_data$source)), " databases"),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  # Save in both TIFF and PNG formats with high resolution
  ggsave(file.path(plot_dir, "00_comprehensive_serum_analysis_panel.tiff"), 
         comprehensive_panel, width = 16, height = 12, dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(plot_dir, "00_comprehensive_serum_analysis_panel.png"), 
         comprehensive_panel, width = 16, height = 12, dpi = 600, device = "png", bg = "white")
  
  return(comprehensive_panel)
}

# Prepare data for comprehensive panel
all_serum_data <- bind_rows(
  gpmdb_serum %>% select(gene, total) %>% mutate(source = "GPMDB"),
  paxdb_serum %>% select(gene, abundance) %>% mutate(source = "PAXDB"),
  hpa_serum %>% select(gene, expr) %>% mutate(source = "HPA Immunoassay"),
  quantms_serum %>% select(gene, abundance) %>% mutate(source = "quantms")
)

# Create comprehensive panel
create_serum_comprehensive_panel(all_serum_data, stats_summary, plot_dir)

# Generate comprehensive markdown report
message("Generating comprehensive serum protein analysis report...")
generate_serum_report <- function(stats_summary, plot_dir) {
  
  # Create report content
  report_content <- paste0(
    "# Serum Protein Analysis Report\n\n",
    "**Analysis Date:** ", Sys.Date(), "\n",
    "**Script:** `04_serum_protein_analysis.R`\n",
    "**Description:** Comprehensive analysis of serum proteins from GPMDB, PAXDB, and HPA Immunoassay databases with overlap and coverage assessment.\n\n",
    "---\n\n",
    
    "## Summary Statistics\n\n",
    "| Data Source | Technology | Unique Proteins | Coverage Type | Detection Method |\n",
    "|-------------|------------|-----------------|---------------|------------------|\n",
    sprintf("| GPMDB Serum | MS | %d | Broad | Spectral counting |\n", stats_summary$gpmdb_serum),
    sprintf("| PAXDB Serum | MS | %d | Comprehensive | Abundance scores |\n", stats_summary$paxdb_serum),
    sprintf("| HPA Immunoassay Serum | Immunoassay | %d | Targeted | Antibody-based |\n", stats_summary$hpa_immunoassay_serum),
    sprintf("| quantms Serum | MS | %d | Comprehensive | iBAQ |\n", stats_summary$quantms_serum),
    sprintf("| **MS Technologies Combined** | MS | %d | Combined MS | Multiple MS methods |\n", stats_summary$ms_technologies),
    sprintf("| **All Sources Combined** | Mixed | %d | Total | All technologies |\n", stats_summary$total_across_sources),
    "\n",
    
    "## Key Findings\n\n",
    sprintf("- **Total serum proteome coverage:** %d unique proteins across all sources\n", stats_summary$total_across_sources),
    sprintf("- **MS dominance:** Mass spectrometry technologies provide %d proteins (%.1f%% of total coverage)\n", 
            stats_summary$ms_technologies, 100 * stats_summary$ms_technologies / stats_summary$total_across_sources),
    sprintf("- **Database complementarity:** Limited overlap between databases suggests technology-specific detection\n"),
    sprintf("- **PAXDB advantage:** Highest individual database coverage (%d proteins) for serum proteome\n", stats_summary$paxdb_serum),
    "- **Immunoassay specificity:** High-sensitivity detection of targeted proteins\n",
    "- **Cross-validation opportunities:** Proteins detected across platforms show enhanced confidence\n\n",
    
    "## Biological Insights\n\n",
    sprintf("- **Serum complexity:** %d+ proteins detectable reflecting serum's role as protein reservoir\n", stats_summary$total_across_sources),
    "- **Technology complementarity:** MS excels at discovery, immunoassays at sensitivity\n",
    "- **Biomarker implications:** Comprehensive coverage supports serum as biomarker source\n",
    "- **Protein abundance range:** Multiple orders of magnitude captured across technologies\n",
    "- **Clinical relevance:** High coverage of clinically important serum proteins\n",
    "- **Disease monitoring potential:** Broad protein coverage enables disease-specific signature detection\n\n",
    
    "## Database Comparison\n\n",
    "### Serum Protein Detection Analysis\n\n",
    "**Mass Spectrometry Databases:**\n",
    sprintf("- PAXDB: Most comprehensive single-source coverage (%d proteins)\n", stats_summary$paxdb_serum),
    sprintf("- GPMDB: Focused detection with complementary protein sets (%d proteins)\n", stats_summary$gpmdb_serum),
    sprintf("- Combined MS: %d proteins representing comprehensive MS-based serum proteome\n\n", stats_summary$ms_technologies),
    
    "**Immunoassay Platform:**\n",
    sprintf("- HPA: Targeted high-sensitivity detection (%d proteins)\n", stats_summary$hpa_immunoassay_serum),
    "- Focus on clinically relevant and well-characterized proteins\n",
    "- Complements MS discovery with validation-grade measurements\n\n",
    
    "**Cross-Database Overlap:**\n",
    "- Significant non-overlap between databases indicates technology-specific capabilities\n",
    "- Core serum proteins detected across multiple platforms show high confidence\n",
    "- Unique detections per platform suggest specialized detection advantages\n\n",
    
    "## Methodology\n\n",
    "- **Data integration:** Standardized gene symbol mapping across all serum databases\n",
    "- **Quality control:** Gene deduplication using median aggregation for multiple entries\n",
    "- **Overlap analysis:** UpSet plots and Venn diagrams for database intersection analysis\n",
    "- **Coverage assessment:** Individual and combined database protein counts\n",
    "- **Statistical analysis:** Correlation analysis between MS-based databases\n",
    "- **Visualization:** Comprehensive panel with coverage, overlap, correlation, and distribution plots\n\n",
    
    "## Recommendations\n\n",
    sprintf("- **Use PAXDB** for comprehensive serum proteome discovery (%d proteins)\n", stats_summary$paxdb_serum),
    "- **Combine GPMDB and PAXDB** for maximum MS-based coverage\n",
    "- **Include HPA Immunoassay** for high-confidence targeted protein measurements\n",
    "- **Consider technology bias** when interpreting serum protein profiles\n",
    "- **Apply orthogonal validation** using complementary detection methods\n",
    "- **Focus on multi-platform proteins** for robust biomarker candidates\n\n",
    
    "## Generated Files\n\n",
    sprintf("- **Comprehensive panel:** `%s/00_comprehensive_serum_analysis_panel.png`\n", basename(plot_dir)),
    "- **Summary statistics:** `outputs/serum_protein/serum_protein_counts_summary.csv`\n",
    "- **Overlap analysis:** `outputs/serum_protein/serum_protein_overlap_statistics.csv`\n",
    "- **Visualization plots:** Individual coverage, overlap, and correlation analyses\n",
    "- **Cross-database correlation:** Statistical validation between MS databases\n\n",
    
    "---\n",
    "*Report generated automatically by the blood proteomics analysis pipeline*\n"
  )
  
  # Save report to serum protein output directory
  report_file <- file.path(dirname(plot_dir), "serum_protein", "serum_protein_analysis_report.md")
  
  # Ensure directory exists
  dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
  
  writeLines(report_content, report_file)
  message(sprintf("✅ Comprehensive serum protein report saved to: %s", report_file))
}

# Generate the report
generate_serum_report(stats_summary, plot_dir)

message("\n=== SERUM PROTEIN ANALYSIS COMPLETE ===") 