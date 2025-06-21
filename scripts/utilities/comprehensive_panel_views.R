#' Comprehensive Panel View Utilities for Blood Proteomics Analysis
#' 
#' This script defines comprehensive panel view functions for each analysis
#' that combine multiple visualizations into cohesive, publication-ready figures
#' 
#' @author Blood Proteomics Analysis Pipeline
#' 

library(patchwork)
library(ggplot2)
library(dplyr)

#' =============================================================================
#'                    01: PLASMA PROTEIN ANALYSIS PANEL
#' =============================================================================

#' Create comprehensive panel for Plasma Protein Analysis (Script 01)
#' 
#' Panel Layout:
#' +-------------------+-------------------+
#' | A: Data Coverage  | B: Technology     |
#' |    Bar Plot       |    Comparison     |
#' +-------------------+-------------------+
#' | C: Distribution Comparison (4 methods)|
#' |    Density plots: Log10 | Z-Score     |
#' |    Quantile Within | Quantile Across  |
#' +---------------------------------------+
#' | D: Consistency Metrics Summary        |
#' +---------------------------------------+
#'
create_plasma_comprehensive_panel <- function(normalized_data, summary_stats, plot_dir) {
  
  # A: Data Coverage by Source
  source_data <- summary_stats %>%
    filter(!source %in% c("Total Across Sources", "MS Technologies Combined")) %>%
    rename(count = unique_genes) %>%
    arrange(desc(count))
  
  panel_A <- create_protein_count_plot(
    source_data, 
    color_by = "technology",
    title = "(A) Plasma Proteins by Database",
    subtitle = "Unique genes detected per database"
  ) + theme(plot.title = element_text(size = 12, face = "bold"))
  
  # B: Technology Classification
  tech_summary <- summary_stats %>%
    filter(!source %in% c("Total Across Sources", "MS Technologies Combined")) %>%
    group_by(technology) %>%
    summarise(protein_count = n_distinct(normalized_data$gene[normalized_data$technology == first(technology)]), .groups = "drop")
  
  panel_B <- create_technology_comparison_plot(
    tech_summary,
    title = "(B) Technology Classification"
  ) + theme(plot.title = element_text(size = 12, face = "bold"))
  
  # C: Distribution Comparison (4 normalization methods)
  methods_data <- normalized_data %>%
    select(gene, source, log_abundance, z_score, quantile_normalized_within, quantile_normalized_across) %>%
    pivot_longer(cols = c(log_abundance, z_score, quantile_normalized_within, quantile_normalized_across),
                names_to = "method", values_to = "value") %>%
    mutate(
      method = factor(method,
                     levels = c("log_abundance", "z_score", "quantile_normalized_within", "quantile_normalized_across"),
                     labels = c("Log10 Original", "Z-Score", "Quantile Within", "Quantile Across"))
    )
  
  panel_C <- ggplot(methods_data, aes(x = value, fill = source)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~ method, scales = "free", ncol = 2) +
    scale_fill_manual(values = get_plot_colors("databases"), name = "Database") +
    theme_blood_proteomics() +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(C) Normalization Methods Comparison",
      subtitle = "Distribution density across databases for each normalization approach",
      x = "Normalized Value",
      y = "Density"
    )
  
  # D: Consistency Metrics Summary
  consistency_data <- data.frame(
    Method = c("Z-Score", "Quantile Within", "Quantile Across"),
    Consistency = c(
      calculate_consistency(normalized_data, "z_score"),
      calculate_consistency(normalized_data, "quantile_normalized_within"),
      calculate_consistency(normalized_data, "quantile_normalized_across")
    )
  ) %>%
    mutate(Method = reorder(Method, -Consistency))
  
  panel_D <- ggplot(consistency_data, aes(x = Method, y = Consistency, fill = Method)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sprintf("%.4f", Consistency)), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    theme_blood_proteomics() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(D) Cross-Database Consistency Scores",
      subtitle = "Lower scores indicate better harmonization (SD of means + SD of medians)",
      x = "Normalization Method",
      y = "Consistency Score"
    )
  
  # Combine panels
  comprehensive_panel <- (panel_A | panel_B) / 
                         panel_C / 
                         panel_D +
                         plot_layout(heights = c(1, 1.5, 1))
  
  comprehensive_panel <- comprehensive_panel +
    plot_annotation(
      title = "Comprehensive Plasma Protein Analysis",
      subtitle = paste0("Cross-database analysis of ", nrow(normalized_data), " protein measurements across ", 
                       length(unique(normalized_data$source)), " databases"),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  save_plot_standard(comprehensive_panel, "00_COMPREHENSIVE_plasma_analysis_panel", plot_dir,
                    width = 16, height = 14)
  
  return(comprehensive_panel)
}

#' =============================================================================
#'                    02: PEPTIDEATLAS QUANTIFICATION PANEL
#' =============================================================================

#' Create comprehensive panel for PeptideAtlas Quantification Analysis (Script 02)
#' 
#' Panel Layout:
#' +---------------------------+---------------------------+
#' | A: Distribution Comparison| B: Correlation Analysis  |
#' |    (Both methods)         |    (Scatter + stats)     |
#' +---------------------------+---------------------------+
#' | C: Method Statistics Summary & Recommendations       |
#' +-------------------------------------------------------+
#'
create_peptideatlas_comprehensive_panel <- function(peptideatlas_clean, plot_dir) {
  
  # A: Distribution Comparison
  dist_data <- peptideatlas_clean %>%
    select(n_observations, norm_PSMs_per_100K) %>%
    tidyr::pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
    mutate(
      metric = case_when(
        metric == "n_observations" ~ "Number of Observations",
        metric == "norm_PSMs_per_100K" ~ "Normalized PSMs per 100K"
      ),
      log_value = log10(value)
    )
  
  panel_A <- ggplot(dist_data, aes(x = metric, y = log_value, fill = metric)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.2, alpha = 0.7, outlier.alpha = 0.3) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
    theme_blood_proteomics() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(A) Quantification Methods Distribution",
      subtitle = "Comparison of abundance measurement approaches",
      x = "Quantification Method",
      y = "Log10(Value)"
    )
  
  # B: Correlation Analysis
  correlation_coef <- cor(log10(peptideatlas_clean$n_observations), 
                         log10(peptideatlas_clean$norm_PSMs_per_100K))
  
  panel_B <- ggplot(peptideatlas_clean, aes(x = log10(n_observations), y = log10(norm_PSMs_per_100K))) +
    geom_point(alpha = 0.6, color = "steelblue", size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
    theme_blood_proteomics() +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    labs(
      title = "(B) Methods Correlation Analysis",
      subtitle = sprintf("Pearson r = %.4f", correlation_coef),
      x = "Log10(Number of Observations)",
      y = "Log10(Normalized PSMs per 100K)"
    )
  
  # C: Summary Statistics Table (as plot)
  stats_summary <- data.frame(
    Metric = c("Sample Size", "Correlation", "Recommendation", "Rationale"),
    Value = c(
      paste(nrow(peptideatlas_clean), "proteins"),
      sprintf("r = %.4f", correlation_coef),
      "Use norm_PSMs_per_100K",
      "Normalized & comparable"
    )
  )
  
  panel_C <- ggplot(stats_summary, aes(x = 1, y = rev(seq_along(Metric)), label = paste(Metric, ":", Value))) +
    geom_text(hjust = 0, size = 4, fontface = c("bold", "plain", "bold", "plain")) +
    theme_void() +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    labs(title = "(C) Analysis Summary & Recommendation") +
    xlim(0, 2) +
    ylim(0, 5)
  
  # Combine panels
  comprehensive_panel <- (panel_A | panel_B) / panel_C +
                         plot_layout(heights = c(2, 1))
  
  comprehensive_panel <- comprehensive_panel +
    plot_annotation(
      title = "Comprehensive PeptideAtlas Quantification Analysis",
      subtitle = "Comparison of n_observations vs norm_PSMs_per_100K quantification methods",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  save_plot_standard(comprehensive_panel, "00_COMPREHENSIVE_peptideatlas_panel", plot_dir,
                    width = 14, height = 10)
  
  return(comprehensive_panel)
}

#' =============================================================================
#'                    03: BIOMARKER PLASMA ANALYSIS PANEL
#' =============================================================================

#' Create comprehensive panel for Biomarker Plasma Analysis (Script 03)
#' 
#' Panel Layout:
#' +-------------------+-------------------+-------------------+
#' | A: Biomarker      | B: Biomarker      | C: Cross-Database |
#' |    Coverage       |    Enrichment     |    Comparison     |
#' +-------------------+-------------------+-------------------+
#' | D: Distribution Comparison (3 normalization methods)     |
#' |    Log10 | Z-Score | Quantile Normalized                |
#' +-------------------------------------------------------+
#'
create_biomarker_comprehensive_panel <- function(all_databases_data, biomarker_genes, plot_dir) {
  
  # A: Biomarker Coverage by Database
  coverage_data <- all_databases_data %>%
    group_by(source) %>%
    summarise(
      total_proteins = n_distinct(gene),
      biomarkers_found = sum(gene %in% biomarker_genes),
      biomarker_percentage = (biomarkers_found / total_proteins) * 100,
      .groups = "drop"
    )
  
  panel_A <- ggplot(coverage_data, aes(x = reorder(source, biomarkers_found), y = biomarkers_found)) +
    geom_col(aes(fill = source), alpha = 0.8) +
    geom_text(aes(label = paste0(biomarkers_found, "\n(", round(biomarker_percentage, 1), "%)")), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    scale_fill_manual(values = get_plot_colors("databases")) +
    theme_blood_proteomics() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(A) Biomarker Coverage",
      subtitle = "Number and percentage of biomarkers per database",
      x = "Database",
      y = "Biomarkers Found"
    )
  
  # B: Biomarker vs Non-biomarker Abundance
  abundance_comparison <- all_databases_data %>%
    mutate(is_biomarker = gene %in% biomarker_genes) %>%
    group_by(source, is_biomarker) %>%
    summarise(
      median_abundance = median(log_abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(protein_type = ifelse(is_biomarker, "Biomarkers", "All Proteins"))
  
  panel_B <- ggplot(abundance_comparison, aes(x = source, y = median_abundance, fill = protein_type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("#69b3a2", "#E69F00")) +
    theme_blood_proteomics() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(B) Biomarker Abundance Enrichment",
      subtitle = "Median log10 abundance comparison",
      x = "Database",
      y = "Median Log10 Abundance",
      fill = "Protein Type"
    )
  
  # C: Cross-Database Biomarker Overlap
  biomarker_overlap <- all_databases_data %>%
    filter(gene %in% biomarker_genes) %>%
    group_by(gene) %>%
    summarise(databases_count = n_distinct(source), .groups = "drop") %>%
    count(databases_count, name = "biomarkers") %>%
    mutate(
      databases_count = factor(databases_count),
      percentage = (biomarkers / sum(biomarkers)) * 100
    )
  
  panel_C <- ggplot(biomarker_overlap, aes(x = databases_count, y = biomarkers)) +
    geom_col(fill = "#2B8CBE", alpha = 0.8) +
    geom_text(aes(label = paste0(biomarkers, "\n(", round(percentage, 1), "%)")), 
              vjust = -0.3, size = 3) +
    theme_blood_proteomics() +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    labs(
      title = "(C) Cross-Database Overlap",
      subtitle = "Biomarkers found in multiple databases",
      x = "Number of Databases",
      y = "Number of Biomarkers"
    )
  
  # D: Distribution Comparison for Biomarkers
  biomarker_distributions <- all_databases_data %>%
    filter(gene %in% biomarker_genes) %>%
    select(gene, source, log_abundance, z_score, quantile_normalized) %>%
    pivot_longer(cols = c(log_abundance, z_score, quantile_normalized),
                names_to = "method", values_to = "value") %>%
    mutate(
      method = factor(method,
                     levels = c("log_abundance", "z_score", "quantile_normalized"),
                     labels = c("Log10 Original", "Z-Score Normalized", "Quantile Normalized"))
    )
  
  panel_D <- ggplot(biomarker_distributions, aes(x = value, fill = source)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~ method, scales = "free", ncol = 3) +
    scale_fill_manual(values = get_plot_colors("databases"), name = "Database") +
    theme_blood_proteomics() +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(D) Biomarker Distribution Across Normalization Methods",
      subtitle = "Density distributions of biomarker abundances",
      x = "Normalized Value",
      y = "Density"
    )
  
  # Combine panels
  comprehensive_panel <- (panel_A | panel_B | panel_C) / panel_D +
                         plot_layout(heights = c(1, 1.2))
  
  comprehensive_panel <- comprehensive_panel +
    plot_annotation(
      title = "Comprehensive Biomarker Plasma Analysis",
      subtitle = paste0("Analysis of ", length(biomarker_genes), " biomarkers across ", 
                       length(unique(all_databases_data$source)), " plasma databases"),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  save_plot_standard(comprehensive_panel, "00_COMPREHENSIVE_biomarker_analysis_panel", plot_dir,
                    width = 16, height = 12)
  
  return(comprehensive_panel)
}

#' =============================================================================
#'                    04: SERUM PROTEIN ANALYSIS PANEL
#' =============================================================================

#' Create comprehensive panel for Serum Protein Analysis (Script 04)
#' 
#' Panel Layout:
#' +-------------------+-------------------+
#' | A: Database       | B: Technology     |
#' |    Coverage       |    Comparison     |
#' +-------------------+-------------------+
#' | C: Database Overlap (UpSet plot)     |
#' +---------------------------------------+
#'
create_serum_comprehensive_panel <- function(all_serum_data, stats_summary, plot_dir) {
  
  message("Creating comprehensive serum analysis panel...")
  
  # A: Database Coverage
  coverage_data <- data.frame(
    Database = c("GPMDB", "PAXDB", "HPA Immunoassay"),
    Count = c(stats_summary$gpmdb_serum, stats_summary$paxdb_serum, stats_summary$hpa_immunoassay_serum),
    Technology = c("MS", "MS", "Immunoassay")
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
  
  # B: Technology Summary
  tech_data <- data.frame(
    Category = c("MS Technologies", "Immunoassay", "All Combined"),
    Count = c(stats_summary$ms_technologies, stats_summary$hpa_immunoassay_serum, stats_summary$total_across_sources),
    Type = c("Combined", "Individual", "Total")
  )
  
  panel_B <- ggplot(tech_data, aes(x = reorder(Category, Count), y = Count, fill = Type)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 3.5) +
    coord_flip() +
    scale_fill_manual(values = c("Combined" = "#56B4E9", "Individual" = "#E15759", "Total" = "#009E73")) +
    theme_blood_proteomics() +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    labs(
      title = "(B) Technology Classification",
      subtitle = "Comparison of detection approaches",
      x = "Category",
      y = "Number of Proteins",
      fill = "Type"
    )
  
  # C: Database Overlap (UpSet plot)
  # Create gene lists for overlap analysis
  gene_lists <- list(
    GPMDB = all_serum_data$gene[all_serum_data$source == "GPMDB"],
    PAXDB = all_serum_data$gene[all_serum_data$source == "PAXDB"],
    HPA_Immunoassay = all_serum_data$gene[all_serum_data$source == "HPA Immunoassay"]
  )
  
  # Create a temporary file for the UpSet plot
  temp_upset_file <- tempfile(fileext = ".png")
  png(temp_upset_file, width = 12, height = 8, units = "in", res = 1200, bg = "white")
  print(UpSetR::upset(
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
  ))
  dev.off()
  
  # Read the UpSet plot as a raster image
  upset_img <- png::readPNG(temp_upset_file)
  upset_grob <- grid::rasterGrob(upset_img, interpolate = TRUE)
  
  # Clean up temporary file
  unlink(temp_upset_file)
  
  # Combine panels
  comprehensive_panel <- (panel_A | panel_B) / 
    wrap_elements(upset_grob) +
    plot_layout(heights = c(1, 1.2)) +
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

#' =============================================================================
#'                    05: CELL TYPE ANALYSIS PANEL
#' =============================================================================

#' Create comprehensive panel for Cell Type Analysis (Script 05)
#' 
#' Panel Layout:
#' +---------------------------+---------------------------+
#' | A: Cell Type Coverage     | B: Source Distribution    |
#' |    (Genes per cell type)  |    (Cell types per source)|
#' +---------------------------+---------------------------+
#' | C: Cell Type Abundance Profiles                      |
#' |    (Heatmap or grouped violin plots)                 |
#' +-------------------------------------------------------+
#' | D: Cross-Source Cell Type Comparison                 |
#' +-------------------------------------------------------+
#'
create_celltype_comprehensive_panel <- function(all_results, overall_summary, plot_dir) {
  
  # A: Cell Type Gene Coverage by Source (Stacked Bar)
  stacked_data <- all_results %>%
    group_by(celltype, source) %>%
    summarise(gene_count = n_distinct(gene), .groups = "drop") %>%
    mutate(celltype_display = format_celltype_names(celltype))

  # Order cell types by total gene count
  celltype_order <- stacked_data %>%
    group_by(celltype_display) %>%
    summarise(total_genes = sum(gene_count)) %>%
    arrange(total_genes) %>%
    pull(celltype_display)

  panel_A <- stacked_data %>%
    mutate(celltype_display = factor(celltype_display, levels = celltype_order)) %>%
    ggplot(aes(x = celltype_display, y = gene_count, fill = source)) +
    geom_col(alpha = 0.85) +
    geom_text(
      data = . %>% group_by(celltype_display) %>% summarise(total = sum(gene_count)),
      aes(x = celltype_display, y = total, label = scales::comma(total)),
      hjust = -0.1, size = 3, inherit.aes = FALSE
    ) +
    coord_flip() +
    scale_fill_viridis_d(option = "turbo", name = "Data Source") +
    theme_blood_proteomics() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(A) Gene Coverage by Cell Type and Source",
      subtitle = "Stacked bars: total unique genes per cell type, colored by contributing data source",
      x = "Cell Type",
      y = "Number of Genes"
    )

  # C: Cell Type Abundance Profiles (Top cell types)
  top_celltypes <- overall_summary %>%
    slice_max(total_genes, n = 8) %>%
    pull(celltype)

  abundance_profiles <- all_results %>%
    filter(celltype %in% top_celltypes) %>%
    mutate(
      log_intensity = log10(intensity + 1),
      celltype_display = format_celltype_names(celltype)
    )

  panel_C <- ggplot(abundance_profiles, aes(x = celltype_display, y = log_intensity, fill = celltype_display)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    scale_fill_viridis_d(option = "turbo", name = "Cell Type") +
    theme_blood_proteomics() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(B) Protein Abundance Profiles",
      subtitle = "Distribution of protein intensities by cell type (top 8)",
      x = "Cell Type",
      y = "Log10(Intensity + 1)"
    )

  # D: Cross-Source Comparison
  cross_source <- all_results %>%
    group_by(celltype, source) %>%
    summarise(
      median_intensity = median(log10(intensity + 1), na.rm = TRUE),
      gene_count = n_distinct(gene),
      .groups = "drop"
    ) %>%
    filter(celltype %in% top_celltypes)

  panel_D <- ggplot(cross_source, aes(x = format_celltype_names(celltype), y = median_intensity, fill = source)) +
    geom_col(position = "dodge", alpha = 0.8) +
    theme_blood_proteomics() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(C) Cross-Source Cell Type Comparison",
      subtitle = "Median protein intensities across data sources",
      x = "Cell Type",
      y = "Median Log10(Intensity + 1)",
      fill = "Data Source"
    )

  # Combine panels: Only new panel_A in top row, then panel_C and panel_D
  comprehensive_panel <- panel_A / panel_C / panel_D +
    plot_layout(heights = c(1.2, 1.2, 1.2))

  comprehensive_panel <- comprehensive_panel +
    plot_annotation(
      title = "Comprehensive Cell Type Analysis",
      subtitle = paste0("Analysis of ", length(unique(all_results$celltype)), " cell types across ", 
                       length(unique(all_results$source)), " data sources"),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  save_plot_standard(comprehensive_panel, "00_COMPREHENSIVE_celltype_analysis_panel", plot_dir,
                    width = 16, height = 16)

  return(comprehensive_panel)
}

#' =============================================================================
#'                    06: INTEGRATION ANALYSIS PANEL
#' =============================================================================

#' Create comprehensive panel for Integration Analysis (Script 06)
#' 
#' Panel Layout:
#' +---------------------------+---------------------------+
#' | A: Method Effectiveness   | B: Consistency Metrics   |
#' |    (Distribution overlay) |    (Bar chart)           |
#' +---------------------------+---------------------------+
#' | C: Cross-Database Harmonization Assessment           |
#' |    (Violin plots by method)                          |
#' +-------------------------------------------------------+
#' | D: Recommendation Summary & Statistics               |
#' +-------------------------------------------------------+
#'
create_integration_comprehensive_panel <- function(normalized_data, effectiveness_results, plot_dir) {
  
  # A: Method Effectiveness (Distribution Comparison)
  methods_data <- normalized_data %>%
    select(gene, source, z_score, quantile_normalized_within, quantile_normalized_across) %>%
    pivot_longer(cols = c(z_score, quantile_normalized_within, quantile_normalized_across),
                names_to = "method", values_to = "value") %>%
    mutate(
      method = factor(method,
                     levels = c("z_score", "quantile_normalized_within", "quantile_normalized_across"),
                     labels = c("Z-Score", "Quantile Within", "Quantile Across"))
    )
  
  panel_A <- ggplot(methods_data, aes(x = value, color = source)) +
    geom_density(alpha = 0.7, size = 0.8) +
    facet_wrap(~ method, scales = "free", ncol = 1) +
    scale_color_manual(values = get_plot_colors("databases"), name = "Database") +
    theme_blood_proteomics() +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(A) Integration Effectiveness",
      subtitle = "Distribution overlay comparison",
      x = "Normalized Value",
      y = "Density"
    )
  
  # B: Consistency Metrics
  consistency_data <- effectiveness_results$consistency %>%
    filter(str_detect(metric, "mean_sd")) %>%
    mutate(
      method = case_when(
        str_detect(metric, "zscore") ~ "Z-Score",
        str_detect(metric, "quantile_within") ~ "Quantile Within",
        str_detect(metric, "quantile_across") ~ "Quantile Across"
      )
    ) %>%
    group_by(method) %>%
    summarise(consistency_score = mean(value), .groups = "drop")
  
  panel_B <- ggplot(consistency_data, aes(x = reorder(method, -consistency_score), y = consistency_score, fill = method)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sprintf("%.4f", consistency_score)), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    theme_blood_proteomics() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(B) Consistency Metrics",
      subtitle = "Lower scores = better integration",
      x = "Normalization Method",
      y = "Consistency Score"
    )
  
  # C: Cross-Database Harmonization
  panel_C <- ggplot(methods_data, aes(x = source, y = value, fill = source)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, alpha = 0.8, outlier.alpha = 0.3) +
    facet_wrap(~ method, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = get_plot_colors("databases"), name = "Database") +
    theme_blood_proteomics() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "(C) Cross-Database Harmonization Assessment",
      subtitle = "Distribution shapes across databases for each method",
      x = "Database",
      y = "Normalized Value"
    )
  
  # D: Recommendation Summary
  best_method <- consistency_data %>%
    slice_min(consistency_score, n = 1) %>%
    pull(method)
  
  best_score <- consistency_data %>%
    slice_min(consistency_score, n = 1) %>%
    pull(consistency_score)
  
  summary_text <- data.frame(
    y = c(4, 3, 2, 1),
    label = c(
      paste("RECOMMENDED METHOD:", best_method),
      paste("Consistency Score:", sprintf("%.4f", best_score)),
      paste("Total Proteins:", nrow(normalized_data)),
      paste("Databases:", length(unique(normalized_data$source)))
    ),
    fontface = c("bold", "plain", "plain", "plain")
  )
  
  panel_D <- ggplot(summary_text, aes(x = 1, y = y, label = label)) +
    geom_text(aes(fontface = fontface), hjust = 0, size = 4) +
    theme_void() +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    labs(title = "(D) Integration Analysis Summary") +
    xlim(0, 3) +
    ylim(0, 5)
  
  # Combine panels
  comprehensive_panel <- (panel_A | panel_B) / panel_C / panel_D +
                         plot_layout(heights = c(1.5, 1.5, 0.8))
  
  comprehensive_panel <- comprehensive_panel +
    plot_annotation(
      title = "Comprehensive Integration Analysis",
      subtitle = "Comparison of normalization approaches for cross-database integration",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  save_plot_standard(comprehensive_panel, "00_COMPREHENSIVE_integration_analysis_panel", plot_dir,
                    width = 16, height = 16)
  
  return(comprehensive_panel)
}

#' =============================================================================
#'                           UTILITY FUNCTIONS
#' =============================================================================

#' Format cell type names for display
format_celltype_names <- function(celltype) {
  str_replace_all(celltype, "_", " ")
}

#' Helper function to calculate consistency (if not already defined)
calculate_consistency <- function(data, value_col) {
  if (!value_col %in% colnames(data)) return(NA)
  
  consistency_stats <- data %>%
    group_by(source) %>%
    summarise(
      mean_val = mean(.data[[value_col]], na.rm = TRUE),
      median_val = median(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  sd_of_means <- sd(consistency_stats$mean_val, na.rm = TRUE)
  sd_of_medians <- sd(consistency_stats$median_val, na.rm = TRUE)
  
  # Composite score. Lower is better.
  sd_of_means + sd_of_medians
}

message("✅ Comprehensive panel view utilities loaded successfully!") 