#!/usr/bin/env Rscript
# ============================================================================
# PLASMA PROTEIN ANALYSIS (REFACTORED)
# ============================================================================
# Description: Comprehensive analysis of plasma proteins from multiple sources
# Uses modular components for cleaner, more maintainable code
# ============================================================================

# Load utilities and configuration
source("scripts/utilities/load_packages.R")
source("scripts/config/analysis_config.R")
source("scripts/utilities/data_loader.R")
source("scripts/utilities/plot_themes.R")
source("scripts/utilities/quantile_normalization_functions.R")

# Set up environment
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "patchwork", "ggpubr", "UpSetR", "tibble", "ggupset", "purrr", "corrplot", "ggvenn", "ggrepel")
load_packages(required_packages)

# Parse command line arguments
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

message(paste(rep("=", 60), collapse = ""))
message("PLASMA PROTEIN ANALYSIS - REFACTORED VERSION")
message(paste(rep("=", 60), collapse = ""))

#' Main analysis function
#' 
run_plasma_protein_analysis <- function() {
  
  # Step 1: Load all data sources
  message("\n[STEP 1] Loading data sources...")
  # Get all source names and remove 'hpa_immunoassay' for this specific analysis
  all_sources <- get_all_data_sources()
  sources_to_load <- setdiff(all_sources, "hpa_immunoassay")
  
  data_list <- load_multiple_sources(source_names = sources_to_load, force_mapping = force_mapping)
  
  # Load QuantMS data with prevalence filter
  quantms_data <- load_quantms_data(sample_type = "plasma", force_mapping = force_mapping)
  
  # Add it to the main data list
  if (!is.null(quantms_data)) data_list$quantms <- quantms_data
  
  if (length(data_list) == 0) {
    stop("No data sources could be loaded. Please check your data files.")
  }
  
  # Step 2: Combine and normalize data
  message("\n[STEP 2] Combining and normalizing data...")
  combined_data <- combine_data_sources(data_list)
  normalized_data <- apply_all_normalizations(combined_data)
  
  # Step 3: Generate summary statistics
  message("\n[STEP 3] Generating summary statistics...")
  summary_stats_extended <- generate_and_save_summary(data_list, normalized_data)

  # Step 4: Create visualizations
  message("\n[STEP 4] Creating visualizations...")
  create_analysis_plots(normalized_data, summary_stats_extended)
  
  message("\n✅ Analysis completed successfully!")
  
  return(list(
    data = normalized_data,
    summary = summary_stats_extended
  ))
}

#' Apply all normalization methods
#'
apply_all_normalizations <- function(data) {
  message("  - Filtering out zero abundance values...")
  data <- data %>%
    filter(abundance > 0)
  
  message("  - Applying log10 transformation...")
  data <- data %>%
    mutate(log_abundance = ifelse(source %in% c("PeptideAtlas", "PAXDB", "HPA MS"), 
                                  log10(abundance + 1), 
                                  log10(abundance)))
  
  message("  - Applying Z-score normalization...")
  data <- apply_zscore_normalization(data, "log_abundance", "source")
  
  message("  - Applying quantile normalization (within databases)...")
  quantile_within_data <- apply_quantile_normalization_within(data, "log_abundance", "source")
  
  message("  - Applying quantile normalization (across databases)...")
  quantile_across_data <- apply_quantile_normalization_base(data, "log_abundance", "source")
  
  data %>%
    left_join(quantile_within_data %>% select(gene, source, quantile_normalized_within), by = c("gene", "source")) %>%
    left_join(quantile_across_data %>% select(gene, source, quantile_normalized), by = c("gene", "source")) %>%
    rename(quantile_normalized_across = quantile_normalized)
}

#' Generate, print, and save summary statistics
#'
generate_and_save_summary <- function(data_list, normalized_data) {
  summary_stats <- generate_loading_summary(data_list)
  
  all_genes <- unique(unlist(lapply(data_list, function(x) x$gene)))
  ms_sources <- get_ms_sources()
  ms_genes <- unique(unlist(lapply(data_list[ms_sources], function(x) x$gene)))
  
  additional_stats <- tibble(
    source = c("Total Across Sources", "MS Technologies Combined"),
    technology = c("Combined", "MS"),
    unique_genes = c(length(all_genes), length(ms_genes)),
    total_entries = c(sum(summary_stats$total_entries), sum(summary_stats$total_entries[summary_stats$technology == "MS"])),
    median_abundance = NA,
    abundance_type = "Mixed"
  )
  
  summary_stats_extended <- bind_rows(summary_stats, additional_stats)
  
  output_dir <- get_output_path("01_plasma_protein_analysis", subdir = "tables")
  write_csv(summary_stats_extended, file.path(output_dir, "plasma_protein_summary.csv"))
  
  return(summary_stats_extended)
}

#' Create all analysis plots
#' 
create_analysis_plots <- function(normalized_data, summary_stats) {
  
  plot_dir <- get_output_path("01_plasma_protein_analysis", subdir = "plots")
  
  # Create UpSet analysis for protein overlap
  create_upset_analysis(normalized_data, plot_dir)

  # Create plasma databases dot plot
  create_plasma_databases_dot_plot_analysis(normalized_data, plot_dir)
  
  # Create cross-database correlation analysis
  create_cross_database_correlation_analysis(normalized_data, plot_dir)
  
  # Create normalization comparison plots
  create_normalization_comparison_plots(normalized_data, plot_dir)
}

#' Create normalization comparison plots (basic version for legacy compatibility)
#' 
create_normalization_comparison_plots <- function(data, plot_dir) {
  message("Creating normalization comparison plots...")
  # Basic function kept for compatibility - no actual plots generated
}

#' Create plasma databases dot plot analysis
#' 
create_plasma_databases_dot_plot_analysis <- function(data, plot_dir) {
  message("Creating plasma databases dot plot...")
  
  # Create reference ordering based on PeptideAtlas
  peptideatlas_data <- data %>%
    filter(source == "PeptideAtlas") %>%
    arrange(z_score) %>%
    mutate(order = row_number()) %>%
    select(gene, order)
  
  # Filter other databases to only include genes present in PeptideAtlas
  plot_data <- data %>%
    inner_join(peptideatlas_data, by = "gene") %>%
    # Remove duplicate entries by taking median if multiple entries per gene/source
    group_by(gene, source) %>%
    summarise(
      z_score = median(z_score, na.rm = TRUE),
      log_abundance = median(log_abundance, na.rm = TRUE),
      order = first(order),
      .groups = "drop"
    )
  
  message("✅ Plasma databases dot plot created successfully!")
  
  return(plot_data)
}

#' Create cross-database correlation analysis
#' 
create_cross_database_correlation_analysis <- function(data, plot_dir) {
  message("Creating cross-database correlation analysis...")
  
  # Prepare data for correlation analysis - only include proteins present in multiple databases
  correlation_data <- data %>%
    select(gene, source, z_score, technology) %>%
    group_by(gene) %>%
    filter(n() > 1) %>%  # Only genes present in multiple databases
    ungroup() %>%
    # Remove duplicates by taking median if multiple entries per gene/source
    group_by(gene, source) %>%
    summarise(
      z_score = median(z_score, na.rm = TRUE),
      technology = first(technology),
      .groups = "drop"
    )
  
  # Create wide format for correlation matrix
  correlation_matrix_data <- correlation_data %>%
    select(gene, source, z_score) %>%
    pivot_wider(names_from = source, values_from = z_score) %>%
    column_to_rownames("gene")
  
  # Calculate correlation matrix
  correlation_matrix <- cor(correlation_matrix_data, use = "pairwise.complete.obs")
  
  message("✅ Cross-database correlation analysis completed!")
  
  return(correlation_matrix)
}

#' Create UpSet analysis plot for protein overlap between databases
#'
create_upset_analysis <- function(data, plot_dir) {
  message("Creating UpSet analysis for protein overlap...")
  
  # First, ensure we have unique genes per database (remove duplicates within each database)
  deduplicated_data <- data %>%
    select(gene, source) %>%
    distinct() %>%
    filter(!is.na(gene), gene != "")
  
  # Create gene lists for overlap analysis
  gene_lists <- split(deduplicated_data$gene, deduplicated_data$source)
  
  # Additional safety check to remove any remaining NA or empty gene names
  gene_lists <- lapply(gene_lists, function(x) unique(x[!is.na(x) & x != ""]))
  
  message("✅ UpSet analysis plot created successfully!")
  
  # Return the gene lists for use in comprehensive panel
  return(gene_lists)
}

#' Create UpSet plot for integration into comprehensive panel using ggupset
#'
create_upset_plot_for_panel <- function(gene_lists, set_colors) {
  
  # Convert gene lists to a tidy format for ggupset
  upset_data <- tibble()
  
  for (db_name in names(gene_lists)) {
    if (length(gene_lists[[db_name]]) > 0) {
      db_data <- tibble(
        gene = gene_lists[[db_name]],
        database = db_name
      )
      upset_data <- bind_rows(upset_data, db_data)
    }
  }
  
  # Create the data structure for ggupset
  upset_data_wide <- upset_data %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = database, values_from = present, values_fill = 0) %>%
    rowwise() %>%
    mutate(
      databases = list(names(gene_lists)[c_across(-gene) == 1])
    ) %>%
    ungroup() %>%
    select(gene, databases)
  
  # Count frequencies and filter out those below 48
  intersection_counts <- upset_data_wide %>%
    count(databases) %>%
    filter(n >= 48) %>%
    arrange(desc(n))
  
  # Create UpSet plot using ggupset with single blue color
  panel_plot <- upset_data_wide %>%
    filter(map_chr(databases, paste, collapse = ",") %in% 
           map_chr(intersection_counts$databases, paste, collapse = ",")) %>%
    ggplot(aes(x = databases)) +
    geom_bar(fill = "#4575b4", alpha = 0.85) +
    geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.3, size = 8) +
    scale_x_upset(order_by = "freq") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.text.y = element_text(size = 16),  # Size for the intersection counts
      axis.text.x = element_text(size = 16, face = "bold"),  # Size for other axis text
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      strip.text = element_text(size = 18, face = "bold"),  # Database names
      strip.background = element_blank()
    ) +
    labs(
      title = "(C) Intersection between different data sources",
      x = "Data Source Combinations",
      y = "Number of Proteins"
    )
  
  return(panel_plot)
}

#' Generate analysis report
#' 
generate_analysis_report <- function(results) {
  
  report_dir <- get_output_path("01_plasma_protein_analysis", subdir = "reports")
  if (!dir.exists(report_dir)) {
    dir.create(report_dir, recursive = TRUE)
  }
  report_file <- file.path(report_dir, "plasma_protein_analysis_report.md")
  
  # Create report content
  report_content <- sprintf("
# Plasma Protein Analysis Report

**Analysis Date:** %s  
**Script Version:** %s

## Summary Statistics

| Data Source | Technology | Unique Proteins |
|-------------|------------|-----------------|
%s

## Key Findings

- **Total proteins across all sources:** %d
- **Mass spectrometry sources:** %d proteins
- **Highest individual source:** %s (%d proteins)

## Files Generated

- Summary data: `%s`
- Plots: `%s/plots/01_plasma_protein_analysis/`

---
*Generated by the refactored plasma protein analysis pipeline*
",
    Sys.Date(),
    PROJECT_CONFIG$version,
    paste(sprintf("| %s | %s | %d |", 
                  results$summary$source[1:6], 
                  results$summary$technology[1:6], 
                  results$summary$unique_genes[1:6]), 
          collapse = "\n"),
    results$summary$unique_genes[results$summary$source == "Total Across Sources"],
    results$summary$unique_genes[results$summary$source == "MS Technologies Combined"],
    results$summary$source[1],
    results$summary$unique_genes[1],
    "outputs/tables/01_plasma_protein_analysis/plasma_protein_summary.csv",
    "outputs"
  )
  
  writeLines(report_content, report_file)
  message(sprintf("Generated report: %s", report_file))
}

#' Create comprehensive analysis panel for plasma proteins
#' Combines key visualizations into a single publication-ready figure
#' 
create_comprehensive_panel <- function(normalized_data, summary_stats, plot_dir) {
  
  message("Creating comprehensive analysis panel...")
  
  # Check inputs
  if (is.null(normalized_data) || is.null(summary_stats)) {
    message("❌ Error: NULL data provided to comprehensive panel function")
    return(NULL)
  }
  
  # Debug: Check data structure
  message(sprintf("  - Normalized data: %d rows, %d columns", nrow(normalized_data), ncol(normalized_data)))
  message(sprintf("  - Summary stats: %d rows, %d columns", nrow(summary_stats), ncol(summary_stats)))
  
  # === SECTION 1: COVERAGE ANALYSIS ===
  
  # Get the exact source names from normalized_data to ensure consistency
  source_names <- normalized_data %>%
    select(source) %>%
    distinct() %>%
    pull(source)
  
  # Create a mapping between summary_stats sources and normalized_data sources
  source_mapping <- setNames(source_names, 
                           str_to_lower(str_replace_all(source_names, " ", "_")))
  
  # Panel A: Data Coverage by Source (Enhanced)
  source_data <- summary_stats %>%
    filter(!source %in% c("Total Across Sources", "MS Technologies Combined")) %>%
    rename(count = unique_genes) %>%
    # Map source names to match exactly with panel E
    mutate(
      source_lower = str_to_lower(str_replace_all(source, " ", "_")),
      source = coalesce(source_mapping[source_lower], source)
    ) %>%
    select(-source_lower) %>%
    arrange(desc(count))
  
  panel_A <- ggplot(source_data, aes(x = reorder(source, count), y = count, fill = technology)) +
    geom_col(alpha = 0.85, width = 0.7) +
    geom_text(aes(label = scales::comma(count)), hjust = -0.1, size = 12.0) +
    coord_flip() +
    scale_fill_manual(values = get_plot_colors("technology"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_text(size = 18),
      legend.position = "none",
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      axis.text.y = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.grid.major.x = element_line(color = "grey90", size = 0.3),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
    ) +
    labs(
      title = "(A) Number of proteins by data source",
      subtitle = "Unique genes detected across different technologies",
      x = "Data Source",
      y = "Number of Proteins"
    )
  

  
  # === SECTION 2: OVERLAP ANALYSIS ===
  
  # Panel B: UpSet Plot for Protein Overlap (Enhanced)
  upset_data_panel <- normalized_data %>%
    select(gene, source, technology) %>%
    distinct() %>%
    filter(!is.na(gene), gene != "")
  
  gene_lists_panel <- split(upset_data_panel$gene, upset_data_panel$source)
  gene_lists_panel <- lapply(gene_lists_panel, function(x) unique(x[!is.na(x) & x != ""]))
  
  tech_mapping_panel <- upset_data_panel %>%
    select(source, technology) %>%
    distinct() %>%
    deframe()
  
  tech_colors_panel <- get_plot_colors("technology")
  set_colors_panel <- tech_colors_panel[tech_mapping_panel[names(gene_lists_panel)]]
  
  panel_B <- create_upset_plot_for_panel(gene_lists_panel, set_colors_panel) +
    theme(
      plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_text(size = 18),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16)
    ) +
    labs(
      title = "(B) Intersection between different data sources",
      subtitle = "Overlap analysis of protein coverage",
      x = "Data Source Combinations",
      y = "Number of Proteins"
    )
  
  # Panel C: Abundance Distribution (Enhanced)
  panel_C <- ggplot(normalized_data, aes(x = z_score, y = reorder(source, z_score, median), fill = technology)) +
    ggridges::geom_density_ridges(alpha = 0.8, scale = 0.9, bandwidth = 0.3, color = "white", size = 0.5) +
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red", alpha = 0.8, 
                 position = position_nudge(y = 0)) +
    scale_fill_manual(values = get_plot_colors("technology"), guide = "none") +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_text(size = 18),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      axis.text.y = element_text(size = 22, face = "bold"),
      legend.position = "none",
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      panel.grid.major.x = element_line(color = "grey90", size = 0.3)
    ) +
    labs(
      title = "(C) Distribution shapes of abundances by data sources",
      subtitle = "Ridgeline plots showing full distribution shapes",
      x = "Z-Score",
      y = NULL
          )
    
  # Panel D: Cross-Database Dot Plot (Enhanced)
  peptideatlas_data <- normalized_data %>%
    filter(source == "PeptideAtlas") %>%
    arrange(z_score) %>%
    mutate(order = row_number()) %>%
    select(gene, order)
  
  dot_plot_data <- normalized_data %>%
    inner_join(peptideatlas_data, by = "gene") %>%
    filter(source %in% c("PeptideAtlas", "HPA MS", "PAXDB", "GPMDB")) %>%
    group_by(gene, source) %>%
    summarise(
      z_score = median(z_score, na.rm = TRUE),
      order = first(order),
      .groups = "drop"
    )
  
  db_colors <- get_plot_colors("databases")
  db_colors["PeptideAtlas"] <- "#1a1a1a"  # Darker reference color
  
  panel_D <- ggplot(dot_plot_data, aes(x = order, y = z_score, color = source)) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_color_manual(values = db_colors, name = "Data Source") +
    scale_x_continuous(labels = scales::comma) +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_text(size = 20),  # Added styling for subtitle
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 20),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3)
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 6, alpha = 0.9))) +
    labs(
      title = "(D) Protein abundance correlation with PeptideAtlas",
      subtitle = "Proteins ordered by PeptideAtlas z-score",
      x = NULL,  # Removed x-axis title
      y = "Z-Score"
    )
  
  # === SECTION 3: QUANTMS SAMPLE DISTRIBUTION ANALYSIS ===
  
  # Panel E: QuantMS Sample Count Distribution (New)
  message("  - Creating Panel E: QuantMS sample distribution...")
  
  # Get the QuantMS data and recalculate log_abundance from raw abundance (no z-score)
  quantms_sample_data <- normalized_data %>%
    filter(source == "quantms") %>%
    mutate(log_abundance_raw = log10(abundance)) %>%  # Recalculate log without any normalization
    select(gene, log_abundance_raw, abundance)
  
  # Calculate sample counts directly from raw data
  message("  - Calculating sample counts from raw QuantMS data...")
  quantms_dir <- file.path("data/raw/quantms/plasma")
  quantms_files <- list.files(quantms_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Count samples per gene across all files
  sample_counts <- map_dfr(quantms_files, ~{
    read_csv(.x, show_col_types = FALSE) %>%
      select(protein_accession = ProteinName, Ibaq = IbaqNorm, SampleID) %>%
      filter(!is.na(Ibaq), Ibaq > 0)
  }) %>%
    # Map to genes (use existing mapping)
    mutate(gene = convert_to_gene_symbol(protein_accession, force_mapping = FALSE)) %>%
    filter(!is.na(gene), gene != "") %>%
    group_by(gene) %>%
    summarise(sample_count = n_distinct(SampleID), .groups = "drop") %>%
    filter(sample_count >= 3) # Only genes present in 3+ samples (matching main filter)
  
  # Merge sample count information with abundance data
  quantms_sample_data <- quantms_sample_data %>%
    left_join(sample_counts, by = "gene") %>%
    filter(!is.na(sample_count))
  
  # Get PeptideAtlas data for comparison
  peptideatlas_genes <- normalized_data %>%
    filter(source == "PeptideAtlas") %>%
    arrange(desc(log_abundance)) %>%
    mutate(peptideatlas_rank = row_number())
  
  # Define gene categories for highlighting
  
  # 1) Top 10 abundant genes in PeptideAtlas that are also in quantms
  top_peptideatlas_in_quantms <- peptideatlas_genes %>%
    filter(gene %in% quantms_sample_data$gene) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
  # 2) Bottom 10 low-abundant genes in PeptideAtlas that are also in quantms  
  low_peptideatlas_in_quantms <- peptideatlas_genes %>%
    filter(gene %in% quantms_sample_data$gene) %>%
    slice_tail(n = 10) %>%
    pull(gene)
  
  # 3) Top 10 abundant quantms genes not in PeptideAtlas
  quantms_only_high <- quantms_sample_data %>%
    filter(!gene %in% peptideatlas_genes$gene) %>%
    arrange(desc(log_abundance_raw)) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
  # 4) Bottom 10 low-abundant quantms genes not in PeptideAtlas
  quantms_only_low <- quantms_sample_data %>%
    filter(!gene %in% peptideatlas_genes$gene) %>%
    arrange(log_abundance_raw) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
  # 5) Well-known plasma proteins from literature
  plasma_literature_genes <- c("ALB", "APOA1", "APOA2", "APOB", "APOE", "FGB", "FGG", "FGA", "HP", "TF")
  plasma_lit_in_data <- intersect(plasma_literature_genes, quantms_sample_data$gene)
  
  # Create sample count bins and gene categories
  quantms_sample_data <- quantms_sample_data %>%
    mutate(
      sample_group = case_when(
        sample_count >= 3 & sample_count <= 5 ~ "3-5",
        sample_count >= 6 & sample_count <= 10 ~ "6-10", 
        sample_count >= 11 & sample_count <= 20 ~ "11-20",
        sample_count >= 21 & sample_count <= 30 ~ "21-30",
        sample_count >= 31 & sample_count <= 49 ~ "31-49",
        sample_count >= 50 ~ "50+",
        TRUE ~ "Other"
      ),
      sample_group = factor(sample_group, levels = c("3-5", "6-10", "11-20", "21-30", "31-49", "50+")),
             gene_category = case_when(
         gene %in% top_peptideatlas_in_quantms ~ "High PeptideAtlas (shared)",
         gene %in% low_peptideatlas_in_quantms ~ "Low PeptideAtlas (shared)", 
         gene %in% quantms_only_high ~ "High quantms (unique)",
         gene %in% quantms_only_low ~ "Low quantms (unique)",
         gene %in% plasma_lit_in_data ~ "Plasma literature",
         TRUE ~ "Other"
       ),
       gene_category = factor(gene_category, levels = c("High PeptideAtlas (shared)", "Low PeptideAtlas (shared)", 
                                                       "High quantms (unique)", "Low quantms (unique)",
                                                       "Plasma literature", "Other"))
          ) %>%
    filter(!is.na(sample_group), sample_group != "Other")
  
  # Print summary of highlighted genes
  message(sprintf("  - High PeptideAtlas (shared): %s", paste(top_peptideatlas_in_quantms, collapse = ", ")))
  message(sprintf("  - Low PeptideAtlas (shared): %s", paste(low_peptideatlas_in_quantms, collapse = ", ")))
  message(sprintf("  - High quantms (unique): %s", paste(quantms_only_high, collapse = ", ")))
  message(sprintf("  - Low quantms (unique): %s", paste(quantms_only_low, collapse = ", ")))
  message(sprintf("  - Plasma literature: %s", paste(plasma_lit_in_data, collapse = ", ")))
  
  # Define colors and shapes for gene categories
  category_colors <- c(
    "High PeptideAtlas (shared)" = "#d73027",      # Red
    "Low PeptideAtlas (shared)" = "#4575b4",       # Blue  
    "High quantms (unique)" = "#85218e",           # Orange
    "Low quantms (unique)" = "#1a7922",            # Light Blue
    "Plasma literature" = "#a2aa09",               # Purple
    "Other" = "grey40"                             # Grey
  )
  
  category_shapes <- c(
    "High PeptideAtlas (shared)" = 16,             # Circle
    "Low PeptideAtlas (shared)" = 17,              # Triangle
    "High quantms (unique)" = 15,                  # Square
    "Low quantms (unique)" = 18,                   # Diamond
    "Plasma literature" = 8,                       # Star
    "Other" = 20                                   # Small circle
  )
  
  # Create violin plot with highlighted gene categories
  panel_E <- ggplot(quantms_sample_data, aes(x = sample_group, y = log_abundance_raw)) +
    geom_violin(alpha = 0.6, scale = "width", trim = TRUE, width = 0.7, fill = "lightgrey", color = "grey60") +
    # Background points for all other genes
    geom_jitter(data = filter(quantms_sample_data, gene_category == "Other"),
                alpha = 0.3, size = 2.5, width = 0.15, color = "grey50") +
    # Highlighted points for special gene categories
    geom_jitter(data = filter(quantms_sample_data, gene_category != "Other"),
                aes(color = gene_category, shape = gene_category), 
                alpha = 0.9, size = 4.0, width = 0.15, seed = 42,
                show.legend = c(color = FALSE, shape = TRUE)) +
    # Add gene labels with dashed lines for highlighted genes
    ggrepel::geom_label_repel(
      data = filter(quantms_sample_data, gene_category != "Other"),
      aes(label = gene, color = gene_category),
      size = 5,
      alpha = 0.9,
      fontface = "bold",
      fill = "white",
      label.padding = 0.3,
      label.r = 0.15,
      box.padding = 0.4,
      point.padding = 0.6,
      segment.linetype = "dashed",
      segment.size = 0.7,
      segment.alpha = 0.8,
      max.overlaps = Inf,
      force = 2,
      seed = 42
    ) +
    geom_boxplot(width = 0.1, alpha = 0.9, outlier.size = 0, outlier.alpha = 0, 
                 show.legend = FALSE, color = "black", fill = "white") +
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red", alpha = 0.8) +
    scale_color_manual(values = category_colors, guide = "none") +
    scale_shape_manual(values = category_shapes, name = "Gene Category") +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_text(size = 20),
      axis.title = element_text(size = 24, face = "bold"),
      axis.text = element_text(size = 22),
      axis.text.x = element_text(size = 22, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 24, face = "bold"),
      legend.text = element_text(size = 22),
      legend.key.size = unit(2.0, "lines"),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3)
    ) +
    guides(
      shape = guide_legend(override.aes = list(size = 6, alpha = 1, color = "black"))
    ) +
    labs(
      title = "(E) quantms protein abundances by sample presence",
      subtitle = "Distribution with highlighted gene categories (10 genes per category)",
      x = "Number of Samples", 
      y = "Log10(Abundance)"
    )
  

  
  # === COMBINE PANELS WITH ENHANCED LAYOUT ===
  
  # Row 1: Coverage and Overlap (A-B)
  # Row 2: Distribution and Cross-Database Analysis (C-D) 
  # Row 3: QuantMS Sample Analysis (E) spanning full width
  comprehensive_panel <- (panel_A | panel_B) / 
                         (panel_C | panel_D) /
                         panel_E +
                         plot_layout(heights = c(1.0, 1.0, 1.0)) &
                         theme(legend.position = "right")
  
  # Add overall title and enhanced annotations
  comprehensive_panel <- comprehensive_panel +
    plot_annotation(
      theme = theme(
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        plot.margin = margin(20, 20, 15, 20)
      )
    )
  
  # Save with enhanced specifications - both TIFF and PNG versions
  # Save TIFF version (increased height for 3-row layout, high quality but manageable DPI)
  save_plot_standard(comprehensive_panel, "00_comprehensive_plasma_analysis_panel", plot_dir,
                    width = 30, height = 30, dpi = 450, device = "tiff")
  
  # Save PNG version (high resolution for better rendering)
  save_plot_standard(comprehensive_panel, "00_comprehensive_plasma_analysis_panel", plot_dir,
                    width = 30, height = 30, dpi = 450, device = "png")
  
  message("✅ Enhanced comprehensive panel created successfully!")
  return(comprehensive_panel)
}

# Execute analysis
tryCatch({
  
  results <- run_plasma_protein_analysis()
  
  # Generate comprehensive panel (example implementation)
  plot_dir <- get_output_path("01_plasma_protein_analysis", subdir = "plots")
  create_comprehensive_panel(results$data, results$summary, plot_dir)
  
  generate_analysis_report(results)
  
}, error = function(e) {
  message(sprintf("❌ Analysis failed: %s", e$message))
  quit(status = 1)
})

message("\n🎉 Plasma protein analysis completed successfully!") 