#!/usr/bin/env Rscript
#
# Analysis of Serum Proteins from Multiple Data Sources
# Script: 04_serum_protein_analysis.R
# 
# Purpose: Analyze the total number of genes/proteins quantified in serum by different sources:
# - GPMDB (Ensembl Protein accessions)
# - PAXDB (Ensembl Protein accessions)  
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

# ---- SET GLOBAL PLOT THEME ----
library(ggplot2)

my_plot_theme <- theme_minimal(base_size = 18) +
  theme(
    plot.title    = element_text(size = 26, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 22, hjust = 0.5),
    axis.title    = element_text(size = 26, face = "bold"),
    axis.text     = element_text(size = 20),
    legend.title  = element_text(size = 22, face = "bold"),
    legend.text   = element_text(size = 20),
    strip.text    = element_text(size = 20, face = "bold"),
    plot.tag      = element_text(size = 22, face = "bold")
  )

# Set as default for all ggplots
theme_set(my_plot_theme)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
force_mapping <- "--force-mapping" %in% args

formats_arg <- NULL
if (length(args) > 0) {
  formats_idx <- which(args == "--formats")
  if (length(formats_idx) > 0 && formats_idx < length(args)) {
    formats_arg <- args[formats_idx + 1]
  }
}

if (!is.null(formats_arg)) set_plot_formats(formats_arg)

# Load gene mapping and deduplication utilities
source("scripts/data_processing/simple_id_mapping.R")
source("scripts/utilities/gene_deduplication.R")
source("scripts/utilities/quantile_normalization_functions.R")

# Set output directories
output_dir <- get_output_path("", subdir = "serum_protein")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- READ AND PROCESS DATA ----
message("Reading and processing serum data from each source...")

# 1. GPMDB Serum
gpmdb_serum_raw <- read_csv("data/raw/gpmdb/gpmdb_serum.csv", show_col_types = FALSE)
gpmdb_serum <- process_gpmdb_data(gpmdb_serum_raw, force_mapping = force_mapping)
gpmdb_serum$source <- "GPMDB"
gpmdb_serum$technology <- "MS"
gpmdb_serum$abundance_type <- "Spectral Count"

# 2. PAXDB Serum
paxdb_serum_raw <- read_csv("data/raw/paxdb/paxdb_serum.csv", show_col_types = FALSE)
paxdb_serum_raw$ensp <- stringr::str_replace(paxdb_serum_raw$string_external_id, "^9606\\.", "")
paxdb_serum_raw$gene <- convert_to_gene_symbol(paxdb_serum_raw$ensp, force_mapping = force_mapping)

paxdb_serum <- deduplicate_genes(
  paxdb_serum_raw, 
  "gene", 
  "abundance", 
  additional_cols = c("ensp"), 
  aggregation_method = "median"
)

# ---- SUMMARY STATISTICS ----
stats_summary <- list(
  gpmdb_serum = length(unique(gpmdb_serum$gene)),
  paxdb_serum = length(unique(paxdb_serum$gene))
)

all_genes <- unique(c(gpmdb_serum$gene, paxdb_serum$gene))
stats_summary$total_across_sources <- length(all_genes)

ms_genes <- unique(c(gpmdb_serum$gene, paxdb_serum$gene))
stats_summary$ms_technologies <- length(ms_genes)

write.csv(
  data.frame(Source = names(stats_summary), Gene_Count = unlist(stats_summary)),
  file.path(output_dir, "serum_protein_counts_summary.csv"),
  row.names = FALSE
)

# ---- CREATE PLOTS ----
message("Creating visualizations...")

plot_dir <- get_output_path("04_serum_protein_analysis", subdir = "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 1. Bar plot by individual sources
source_data <- data.frame(
  Source = c("GPMDB", "PAXDB"),
  Count = c(stats_summary$gpmdb_serum, stats_summary$paxdb_serum),
  Technology = c("MS", "MS")
)

p1 <- ggplot(source_data, aes(x = reorder(Source, Count), y = Count, fill = Technology)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 5) +
  coord_flip() +
  scale_fill_manual(values = c("MS" = "#2E86AB", "Immunoassay" = "#F18F01")) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Serum Proteins Quantified by Data Source",
    subtitle = "Number of unique genes detected in each database",
    x = "Data Source",
    y = "Number of Genes",
    fill = "Technology"
  )

# 2. Grouped bar plot
summary_data <- data.frame(
  Category = c("GPMDB", "PAXDB", "MS Technologies", "All Sources"),
  Count = c(stats_summary$gpmdb_serum, stats_summary$paxdb_serum, stats_summary$ms_technologies, stats_summary$total_across_sources),
  Type = c("Individual", "Individual", "Combined", "Combined")
)

p2 <- ggplot(summary_data, aes(x = reorder(Category, Count), y = Count, fill = Type)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 5) +
  coord_flip() +
  scale_fill_manual(values = c("Individual" = "#4E79A7", "Combined" = "#E15759")) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Serum Protein Detection: Individual vs Combined Sources",
    subtitle = "Comparison of gene counts across databases and technologies",
    x = "Database/Category",
    y = "Number of Genes",
    fill = "Category Type"
  )

# 3. UpSet plot (overlap)
gene_lists <- list(
  GPMDB = gpmdb_serum$gene[!is.na(gpmdb_serum$gene) & gpmdb_serum$gene != ""],
  PAXDB = paxdb_serum$gene[!is.na(paxdb_serum$gene) & paxdb_serum$gene != ""]
)

p3 <- UpSetR::upset(
  fromList(gene_lists),
  nsets = 4,
  sets = c("GPMDB", "PAXDB"),
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  text.scale = 1.2,
  point.size = 3,
  line.size = 1.2,
  mainbar.y.label = "Number of Genes",
  sets.x.label = "Total Genes per Database"
)

# 4. Distribution plots
calculate_zscore_normalization_serum <- function(data) {
  data %>%
    group_by(Database) %>%
    mutate(z_score = (log_value - mean(log_value, na.rm = TRUE)) / sd(log_value, na.rm = TRUE)) %>%
    ungroup()
}

gpmdb_dist <- gpmdb_serum %>% filter(!is.na(total) & total > 0) %>%
  mutate(Database = "GPMDB", log_value = log10(total))

paxdb_dist <- paxdb_serum %>% filter(!is.na(abundance) & abundance > 0) %>%
  mutate(Database = "PAXDB", log_value = log10(abundance))

dist_data <- bind_rows(
  gpmdb_dist %>% select(gene, Database, log_value),
  paxdb_dist %>% select(gene, Database, log_value)
)

dist_data_zscore <- calculate_zscore_normalization_serum(dist_data)

p4 <- ggplot(dist_data, aes(x = log_value, fill = Database)) +
  geom_histogram(alpha = 0.7, bins = 50) +
  facet_wrap(~Database, scales = "free", ncol = 1,
             labeller = labeller(Database = c("GPMDB" = "(a) GPMDB", "PAXDB" = "(b) PAXDB"))) +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Distribution of Protein Quantification Values",
    subtitle = "Unique proteins per database",
    x = "Log10(Quantification Value)",
    y = "Number of Proteins"
  )

p4z <- ggplot(dist_data_zscore, aes(x = z_score, fill = Database)) +
  geom_histogram(alpha = 0.7, bins = 50) +
  facet_wrap(~Database, scales = "free", ncol = 1,
             labeller = labeller(Database = c("GPMDB" = "(a) GPMDB", "PAXDB" = "(b) PAXDB"))) +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Distribution of Z-Score Normalized Protein Quantification Values",
    subtitle = "Z-score normalized values across databases",
    x = "Z-Score (standardized within each database)",
    y = "Number of Proteins"
  )

# ---- COMPREHENSIVE PANEL FUNCTION ----
create_serum_comprehensive_panel <- function(all_serum_data, stats_summary, plot_dir) {
  message("Creating comprehensive serum analysis panel...")
  
  # Panel A: Database coverage
  coverage_data <- data.frame(
    Database = c("GPMDB", "PAXDB"),
    Count = c(stats_summary$gpmdb_serum, stats_summary$paxdb_serum),
    Technology = c("MS", "MS")
  ) %>% arrange(desc(Count))
  
  panel_A <- ggplot(coverage_data, aes(x = reorder(Database, Count), y = Count, fill = Technology)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = scales::comma(Count)), hjust = -0.05, size = 6) +  
    coord_flip(clip = "off") +
    scale_fill_manual(values = c("MS" = "#4E79A7", "Immunoassay" = "#E15759")) +
    labs(
      title = "(A) Serum Database Coverage",
      subtitle = "Unique proteins per database",
      x = "Database",
      y = "Number of Proteins",
      fill = "Technology"
    ) +
    theme(plot.margin = margin(t = 5, r = 25, b = 5, l = 5))
  
  
  library(ggplotify)  # allows conversion of ggvenn to ggplot for theming
  
  library(ggplotify)
  
  library(ggplotify)
  
  # Panel B: Venn
  gene_lists <- list(
    GPMDB = all_serum_data$gene[all_serum_data$source == "GPMDB"],
    PAXDB = all_serum_data$gene[all_serum_data$source == "PAXDB"]
  )
  
  total_unique <- length(unique(unlist(gene_lists)))
  overlap_all <- length(Reduce(intersect, gene_lists))
  overlap_percent <- round(overlap_all / total_unique * 100, 1)
  
  panel_B <- ggvenn(
    gene_lists,
    fill_color = c("#4E79A7", "#56B4E9"),
    stroke_size = 0.5,
    set_name_size = 8,  # smaller titles
    text_size = 6       # larger numbers inside circles
  ) %>% 
    as.ggplot() +
    labs(
      title = "(B) Database Overlap"
    ) +
    theme(
      plot.title = element_text(size = 26, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 22, hjust = 0.5)
    )
  
  
  
  
  # Panel C: Correlation heatmap
  ms_data <- all_serum_data %>%
    filter(source %in% c("GPMDB", "PAXDB")) %>%
    mutate(value = case_when(
      source == "GPMDB" ~ total,
      source == "PAXDB" ~ abundance,
      TRUE ~ NA_real_
    )) %>%
    mutate(log_value = log10(value + 1)) %>%
    group_by(source) %>% 
    mutate(z_score = scale(log_value)[,1]) %>% 
    ungroup() %>%
    select(gene, source, z_score)
  
  correlation_matrix_data <- ms_data %>% 
    pivot_wider(names_from = source, values_from = z_score) %>% 
    column_to_rownames("gene")
  
  correlation_matrix <- cor(correlation_matrix_data, use = "pairwise.complete.obs")
  
  corr_data <- as.data.frame(correlation_matrix) %>% 
    rownames_to_column("db1") %>% 
    pivot_longer(-db1, names_to = "db2", values_to = "correlation")
  
  panel_C <- ggplot(corr_data, aes(x = db1, y = db2, fill = correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 8) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4", midpoint = 0, limit = c(-1,1)) +
    labs(
      title = "(C) Cross-Database Correlation",
      subtitle = "Pearson correlation (MS databases)",
      fill = "Correlation"
    ) +
    theme_minimal(base_size = 18) +  # base font size
    theme(
      axis.text.x  = element_text(size = 20, angle = 45, hjust = 1, face = "bold"),
      axis.text.y  = element_text(size = 20, face = "bold"),
      axis.title   = element_text(size = 22, face = "bold"),
      plot.title   = element_text(size = 26, face = "bold", hjust = 0.5),
      plot.subtitle= element_text(size = 22, hjust = 0.5),
      legend.title = element_text(size = 22, face = "bold"),
      legend.text  = element_text(size = 20),
      plot.margin  = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  
  
  # Panel D: Violin
  distribution_data <- all_serum_data %>%
    mutate(value = case_when(source == "GPMDB" ~ total,
                             source == "PAXDB" ~ abundance,
                             TRUE ~ NA_real_)) %>%
    group_by(source) %>% mutate(log_value = log10(value + 1), z_score = scale(log_value)[,1]) %>% ungroup()
  
  panel_D <- ggplot(distribution_data, aes(x = source, y = z_score, fill = source)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    scale_fill_manual(values = c("GPMDB" = "#4E79A7", "PAXDB" = "#56B4E9")) +
    labs(
      title = "(D) Abundance Distribution",
      subtitle = "Z-score normalized values across databases",
      x = "Database",
      y = "Z-score"
    )
  
  # Combine panels
  comprehensive_panel <- (panel_A | panel_B) / (panel_C | panel_D) +
    plot_layout(heights = c(1,1)) +
    plot_annotation(
      title = "Comprehensive Serum Protein Analysis",
      subtitle = paste0("Cross-database analysis of serum proteins across ", length(unique(all_serum_data$source)), " databases")
    )
  
  # Save
  ggsave(file.path(plot_dir, "00_comprehensive_serum_analysis_panel.tiff"), comprehensive_panel, width = 16, height = 12, dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(plot_dir, "00_comprehensive_serum_analysis_panel.png"), comprehensive_panel, width = 16, height = 12, dpi = 600, device = "png", bg = "white")
  
  return(comprehensive_panel)
}

all_serum_data <- bind_rows(
  gpmdb_serum %>% select(gene, total) %>% mutate(source = "GPMDB"),
  paxdb_serum %>% select(gene, abundance) %>% mutate(source = "PAXDB")
)

create_serum_comprehensive_panel(all_serum_data, stats_summary, plot_dir)

# ---- GENERATE REPORT ----
generate_serum_report <- function(stats_summary, plot_dir) {
  report_content <- paste0(
    "# Serum Protein Analysis Report\n\n",
    "**Analysis Date:** ", Sys.Date(), "\n",
    "**Script:** `04_serum_protein_analysis.R`\n",
    "**Description:** Comprehensive analysis of serum proteins from GPMDB, PAXDB databases.\n\n",
    "---\n\n",
    "## Summary Statistics\n\n",
    "| Data Source | Technology | Unique Proteins |\n",
    "|-------------|------------|----------------|\n",
    sprintf("| GPMDB Serum | MS | %d |\n", stats_summary$gpmdb_serum),
    sprintf("| PAXDB Serum | MS | %d |\n", stats_summary$paxdb_serum),
    sprintf("| **MS Technologies Combined** | MS | %d |\n", stats_summary$ms_technologies),
    sprintf("| **All Sources Combined** | Mixed | %d |\n", stats_summary$total_across_sources),
    "\n"
  )
  
  report_file <- file.path(dirname(plot_dir), "serum_protein", "serum_protein_analysis_report.md")
  dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(report_content, report_file)
  message(sprintf("Comprehensive serum protein report saved to: %s", report_file))
}

generate_serum_report(stats_summary, plot_dir)
message("\n=== SERUM PROTEIN ANALYSIS COMPLETE ===")
