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
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "patchwork", "ggpubr", "UpSetR", "tibble", "ggupset", "purrr", "corrplot", "ggvenn")
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
  data_list <- load_multiple_sources(force_mapping = force_mapping)
  
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
  message("  - Applying log10 transformation...")
  data <- data %>%
    mutate(log_abundance = log10(abundance + 1))
  
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

#' Calculate consistency score for a normalization method
#' A lower score indicates better consistency across data sources.
#' @param data A data frame containing the normalized data.
#' @param value_col The name of the column with the normalization values (e.g., "z_score").
#' @return A numeric consistency score (sum of SD of means and SD of medians).
calculate_consistency <- function(data, value_col) {
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
  
  # Calculate and print consistency scores
  consistency_z <- calculate_consistency(normalized_data, "z_score")
  consistency_q_within <- calculate_consistency(normalized_data, "quantile_normalized_within")
  consistency_q_across <- calculate_consistency(normalized_data, "quantile_normalized_across")
  
  message("\n--- Normalization Consistency ---")
  message(sprintf("Z-Score Consistency (lower is better): %.4f", consistency_z))
  message(sprintf("Quantile Within Consistency (lower is better): %.4f", consistency_q_within))
  message(sprintf("Quantile Across Consistency (lower is better): %.4f", consistency_q_across))
  message("---------------------------------")
  
  return(summary_stats_extended)
}


#' Create all analysis plots
#' 
create_analysis_plots <- function(normalized_data, summary_stats) {
  
  plot_dir <- get_output_path("01_plasma_protein_analysis", subdir = "plots")
  
  # Plot 1: Individual source counts
  source_data <- summary_stats %>%
    filter(!source %in% c("Total Across Sources", "MS Technologies Combined")) %>%
    rename(count = unique_genes) %>%
    arrange(desc(count))
  
  p1 <- create_protein_count_plot(
    source_data, 
    color_by = "technology",
    title = "Plasma Proteins Quantified by Data Source",
    subtitle = "Number of unique genes detected in each database"
  )
  
  save_plot_standard(p1, "01_plasma_proteins_by_source", plot_dir)
  
  # Plot 2: Technology comparison
  tech_summary <- summary_stats %>%
    filter(!source %in% c("Total Across Sources", "MS Technologies Combined")) %>%
    group_by(technology) %>%
    summarise(protein_count = n_distinct(normalized_data$gene[normalized_data$technology == first(technology)]), .groups = "drop")

  p2 <- create_technology_comparison_plot(
    tech_summary,
    title = "Plasma Proteins by Technology Classification"
  )
  
  save_plot_standard(p2, "02_plasma_proteins_by_technology", plot_dir)
  
  # Plot 3: UpSet analysis for protein overlap
  p3_upset <- create_upset_analysis(normalized_data, plot_dir)

  # Plot 5: Plasma databases dot plot
  create_plasma_databases_dot_plot_analysis(normalized_data, plot_dir)
  
  # Plot 6: Cross-database correlation analysis
  create_cross_database_correlation_analysis(normalized_data, plot_dir)
  
  # Plot 7: Normalization distributions
  create_normalization_comparison_plots(normalized_data, plot_dir)
}

#' Create abundance distribution analysis plots for all normalization methods
#' 
create_normalization_comparison_plots <- function(data, plot_dir) {
  
  message("Creating normalization comparison plots...")
  
  # Violin plots
  p_v_log <- create_enhanced_violin_plot(data, "source", "log_abundance", "technology", "Abundance Distribution (Log10)", "Log10-transformed values")
  p_v_z <- create_enhanced_violin_plot(data, "source", "z_score", "technology", "Abundance Distribution (Z-Score)", "Within-database standardized")
  p_v_q_within <- create_enhanced_violin_plot(data, "source", "quantile_normalized_within", "technology", "Abundance Distribution (Quantile Within)", "Within-database quantile normalized")
  p_v_q_across <- create_enhanced_violin_plot(data, "source", "quantile_normalized_across", "technology", "Abundance Distribution (Quantile Across)", "Cross-database quantile normalized")
  
  save_plot_standard(p_v_log, "07_dist_violin_log10", plot_dir)
  save_plot_standard(p_v_z, "08_dist_violin_zscore", plot_dir)

  
  # Density plots
  p_d_log <- create_density_plot(data, "log_abundance", "source", "Quantification Density (Log10)")
  p_d_z <- create_density_plot(data, "z_score", "source", "Quantification Density (Z-Score)")
  p_d_q_within <- create_density_plot(data, "quantile_normalized_within", "source", "Quantification Density (Quantile Within)")
  p_d_q_across <- create_density_plot(data, "quantile_normalized_across", "source", "Quantification Density (Quantile Across)")
  
  combined_density <- (p_d_log + p_d_z) / (p_d_q_within + p_d_q_across) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save_plot_standard(combined_density, "11_dist_density_comparison", plot_dir, width = 14, height = 12)


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
  
  # Create the dot plot
  p_dot <- create_plasma_databases_dot_plot(plot_data)
  
  # Save the plot
  save_plot_standard(p_dot, "05_plasma_databases_comparison", plot_dir,
                    width = PLOT_CONFIG$dimensions$large_width,
                    height = PLOT_CONFIG$dimensions$default_height)
  
  message("✅ Plasma databases dot plot created successfully!")
  
  return(plot_data)
}

#' Create plasma databases dot plot
#' 
#' @param data Data frame with normalized protein data
#' @return ggplot object
#' 
create_plasma_databases_dot_plot <- function(data) {
  
  # Get database colors, with PeptideAtlas in black as reference
  db_colors <- get_plot_colors("databases")
  db_colors["PeptideAtlas"] <- "black"
  
  # Create the plot
  p <- ggplot(data, aes(x = order, y = z_score, color = source)) +
    geom_point(alpha = 0.7, size = 1) +
    scale_color_manual(values = db_colors, name = "Database") +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = "Plasma Protein Abundance Comparison Across Databases",
      subtitle = "Proteins ordered by PeptideAtlas z-score (reference database in black)",
      x = "Proteins (sorted by PeptideAtlas z-score)",
      y = "Z-score normalized abundance",
      caption = "Only proteins detected in PeptideAtlas are shown for comparison"
    ) +
    theme_blood_proteomics() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    )
  
  return(p)
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
  
  # Create correlation heatmap
  p_heatmap <- create_correlation_heatmap(correlation_matrix, data)
  
  # Create pairwise correlation scatter plots
  p_scatter <- create_pairwise_correlation_plots(correlation_data)
  
  # Create technology-based correlation summary
  p_tech_corr <- create_technology_correlation_summary(correlation_data)
  
  # Save all plots
  save_plot_standard(p_heatmap, "06_correlation_heatmap", plot_dir,
                    width = 10, height = 8)
  
  save_plot_standard(p_scatter, "06_pairwise_correlations", plot_dir,
                    width = 14, height = 10)
  
  save_plot_standard(p_tech_corr, "06_technology_correlations", plot_dir,
                    width = 12, height = 8)
  
  message("✅ Cross-database correlation analysis completed!")
  
  return(correlation_matrix)
}

#' Create correlation heatmap
#' 
create_correlation_heatmap <- function(correlation_matrix, data) {
  
  # Convert correlation matrix to long format for ggplot
  corr_data <- correlation_matrix %>%
    as.data.frame() %>%
    rownames_to_column("database1") %>%
    pivot_longer(cols = -database1, names_to = "database2", values_to = "correlation")
  
  # Add technology information
  tech_mapping <- data %>%
    select(source, technology) %>%
    distinct() %>%
    rename(database1 = source, tech1 = technology)
  
  tech_mapping2 <- tech_mapping %>%
    rename(database2 = database1, tech2 = tech1)
  
  corr_data <- corr_data %>%
    left_join(tech_mapping, by = "database1") %>%
    left_join(tech_mapping2, by = "database2") %>%
    mutate(
      same_technology = tech1 == tech2,
      correlation_type = case_when(
        database1 == database2 ~ "Self",
        same_technology ~ "Within Technology",
        !same_technology ~ "Cross Technology"
      )
    )
  
  # Create heatmap
  p <- ggplot(corr_data, aes(x = database1, y = database2, fill = correlation)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f", correlation)), 
              color = ifelse(abs(corr_data$correlation) > 0.6, "white", "black"),
              size = 3) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4",
                        midpoint = 0, limit = c(-1, 1), space = "Lab",
                        name = "Pearson\nCorrelation") +
    theme_blood_proteomics() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    ) +
    labs(
      title = "Cross-Database Z-Score Correlation Matrix",
      subtitle = "Correlation of protein abundances after z-score normalization",
      caption = "Only proteins present in multiple databases included"
    ) +
    coord_fixed()
  
  return(p)
}

#' Create correlation heatmap for comprehensive panel (simplified version)
#' 
create_correlation_heatmap_for_panel <- function(correlation_matrix, data) {
  
  # Convert correlation matrix to long format for ggplot
  corr_data <- correlation_matrix %>%
    as.data.frame() %>%
    rownames_to_column("database1") %>%
    pivot_longer(cols = -database1, names_to = "database2", values_to = "correlation")
  
  # Add technology information
  tech_mapping <- data %>%
    select(source, technology) %>%
    distinct() %>%
    rename(database1 = source, tech1 = technology)
  
  tech_mapping2 <- tech_mapping %>%
    rename(database2 = database1, tech2 = tech1)
  
  corr_data <- corr_data %>%
    left_join(tech_mapping, by = "database1") %>%
    left_join(tech_mapping2, by = "database2") %>%
    mutate(
      same_technology = tech1 == tech2,
      correlation_type = case_when(
        database1 == database2 ~ "Self",
        same_technology ~ "Within Technology",
        !same_technology ~ "Cross Technology"
      )
    )
  
  # Create enhanced heatmap for panel
  p <- ggplot(corr_data, aes(x = database1, y = database2, fill = correlation)) +
    geom_tile(color = "white", size = 0.4) +
    geom_text(aes(label = sprintf("%.2f", correlation)), 
              color = ifelse(abs(corr_data$correlation) > 0.6, "white", "black"),
              size = 8) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a73e8",
                        midpoint = 0, limit = c(-1, 1), space = "Lab",
                        guide = "none") +
    theme_blood_proteomics() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "none"
    ) +
    labs(
      title = "(G) CROSS-DATABASE CORRELATIONS",
      subtitle = "Pairwise correlation matrix (z-scores)",
      x = NULL,
      y = NULL
    ) +
    coord_fixed()
  
  return(p)
}

#' Create pairwise correlation scatter plots
#' 
create_pairwise_correlation_plots <- function(data) {
  
  # Create all pairwise combinations
  databases <- unique(data$source)
  combinations <- combn(databases, 2, simplify = FALSE)
  
  # Limit to most important comparisons to avoid too many plots
  important_combinations <- combinations[1:min(6, length(combinations))]
  
  plot_list <- map(important_combinations, function(pair) {
    
    pair_data <- data %>%
      filter(source %in% pair) %>%
      select(gene, source, z_score) %>%
      pivot_wider(names_from = source, values_from = z_score) %>%
      filter(!is.na(.data[[pair[1]]]) & !is.na(.data[[pair[2]]]))
    
    if (nrow(pair_data) > 10) {  # Only create plot if enough shared proteins
      correlation <- cor(pair_data[[pair[1]]], pair_data[[pair[2]]], use = "complete.obs")
      
      ggplot(pair_data, aes(x = .data[[pair[1]]], y = .data[[pair[2]]])) +
        geom_point(alpha = 0.6, size = 1) +
        geom_smooth(method = "lm", se = TRUE, color = "#2E86AB") +
        labs(
          title = paste0(pair[1], " vs ", pair[2]),
          subtitle = paste0("r = ", sprintf("%.3f", correlation), 
                           " (n = ", nrow(pair_data), " proteins)"),
          x = paste0(pair[1], " (Z-Score)"),
          y = paste0(pair[2], " (Z-Score)")
        ) +
        theme_blood_proteomics() +
        theme(plot.title = element_text(size = 10))
    } else {
      NULL
    }
  })
  
  # Remove NULL plots and combine
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  if (length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = 3) +
      plot_annotation(
        title = "Pairwise Database Correlations (Z-Score Normalized)",
        subtitle = "Scatter plots showing correlations between database pairs"
      )
    return(combined_plot)
  } else {
    return(NULL)
  }
}

#' Create technology-based correlation summary
#' 
create_technology_correlation_summary <- function(data) {
  
  # Calculate all pairwise correlations and categorize by technology
  databases <- unique(data$source)
  combinations <- combn(databases, 2, simplify = FALSE)
  
  correlation_summary <- map_dfr(combinations, function(pair) {
    
    pair_data <- data %>%
      filter(source %in% pair) %>%
      select(gene, source, z_score, technology) %>%
      pivot_wider(names_from = source, values_from = z_score) %>%
      filter(!is.na(.data[[pair[1]]]) & !is.na(.data[[pair[2]]]))
    
    if (nrow(pair_data) > 5) {
      correlation <- cor(pair_data[[pair[1]]], pair_data[[pair[2]]], use = "complete.obs")
      
      # Get technologies for each database
      tech1 <- data$technology[data$source == pair[1]][1]
      tech2 <- data$technology[data$source == pair[2]][1]
      
      tibble(
        database1 = pair[1],
        database2 = pair[2],
        tech1 = tech1,
        tech2 = tech2,
        correlation = correlation,
        n_proteins = nrow(pair_data),
        comparison_type = ifelse(tech1 == tech2, "Within Technology", "Cross Technology")
      )
    } else {
      NULL
    }
  }) %>%
  filter(!is.null(.))
  
  # Create summary plot
  p <- ggplot(correlation_summary, aes(x = comparison_type, y = correlation, fill = comparison_type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    scale_fill_manual(values = c("Within Technology" = "#4575b4", "Cross Technology" = "#d73027")) +
    theme_blood_proteomics() +
    theme(legend.position = "none") +
    labs(
      title = "Z-Score Correlation by Technology Type",
      subtitle = "Distribution of correlations within vs across technology types",
      x = "Comparison Type",
      y = "Pearson Correlation",
      caption = "Each point represents a pairwise database comparison"
    ) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
                fill = "white", color = "black")
  
  return(p)
}

#' Create Venn diagram for GPMDB vs PeptideAtlas comparison using ggvenn
#' 
create_ms_venn_diagram <- function(data) {
  
  message("Creating MS databases Venn diagram with ggvenn (GPMDB vs PeptideAtlas)...")
  
  # Get gene lists for GPMDB and PeptideAtlas
  gpmdb_genes <- data %>%
    filter(source == "GPMDB") %>%
    pull(gene) %>%
    unique()
  
  peptideatlas_genes <- data %>%
    filter(source == "PeptideAtlas") %>%
    pull(gene) %>%
    unique()
  
  # Create Venn diagram data as named list
  venn_data <- list(
    GPMDB = gpmdb_genes,
    PeptideAtlas = peptideatlas_genes
  )
  
  # Calculate statistics for annotations
  intersection_count <- length(intersect(gpmdb_genes, peptideatlas_genes))
  gpmdb_unique <- length(gpmdb_genes) - intersection_count
  peptideatlas_unique <- length(peptideatlas_genes) - intersection_count
  overlap_percent <- round(intersection_count / length(union(gpmdb_genes, peptideatlas_genes)) * 100, 1)
  
  # Create beautiful Venn diagram with ggvenn
  venn_plot <- ggvenn(
    venn_data,
    fill_color = c("#E31A1C", "#1F78B4"),
    stroke_color = "white",
    stroke_size = 0.8,
    stroke_alpha = 1,
    fill_alpha = 0.7,
    set_name_color = "black",
    set_name_size = 5,
    text_color = "black",
    text_size = 4.5
  ) +
  theme_void() +
  labs(
    title = "MS Databases Protein Overlap Analysis",
    subtitle = sprintf("GPMDB (%s proteins) vs PeptideAtlas (%s proteins)", 
                      scales::comma(length(gpmdb_genes)), 
                      scales::comma(length(peptideatlas_genes))),
    caption = sprintf("Shared: %s proteins (%.1f%%) | GPMDB unique: %s | PeptideAtlas unique: %s", 
                     scales::comma(intersection_count), 
                     overlap_percent,
                     scales::comma(gpmdb_unique),
                     scales::comma(peptideatlas_unique))
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
    plot.caption = element_text(size = 11, hjust = 0.5, color = "gray50", lineheight = 1.2, margin = margin(t = 15)),
    plot.margin = margin(20, 20, 20, 20),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
  
  return(venn_plot)
}

#' Create Venn diagram for panel using ggvenn
#' 
create_ms_venn_diagram_for_panel <- function(data) {
  
  message("Creating enhanced MS databases Venn diagram with ggvenn (GPMDB vs PeptideAtlas)...")
  
  # Get gene lists for GPMDB and PeptideAtlas
  gpmdb_genes <- data %>%
    filter(source == "GPMDB") %>%
    pull(gene) %>%
    unique()
  
  peptideatlas_genes <- data %>%
    filter(source == "PeptideAtlas") %>%
    pull(gene) %>%
    unique()
  
  # Create Venn diagram data as named list
  venn_data <- list(
    GPMDB = gpmdb_genes,
    PeptideAtlas = peptideatlas_genes
  )
  
  # Calculate statistics for annotations
  intersection_count <- length(intersect(gpmdb_genes, peptideatlas_genes))
  overlap_percent <- round(intersection_count / length(union(gpmdb_genes, peptideatlas_genes)) * 100, 1)
  
  # Create enhanced Venn diagram with ggvenn for panel
  venn_plot <- ggvenn(
    venn_data,
    fill_color = c("#E31A1C", "#1F78B4"),
    stroke_color = "white",
    stroke_size = 0.8,
    stroke_alpha = 1,
    fill_alpha = 0.75,
    set_name_color = "black",
    set_name_size = 4,
    text_color = "black",
    text_size = 3.8
  ) +
  theme_void() +
  labs(
    title = "(D) MS DATABASES COMPARISON",
    subtitle = sprintf("GPMDB (%s) vs PeptideAtlas (%s)", 
                      scales::comma(length(gpmdb_genes)), 
                      scales::comma(length(peptideatlas_genes))),
    caption = sprintf("Shared: %s proteins (%.1f%% overlap)", 
                     scales::comma(intersection_count), overlap_percent)
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "#2c3e50"),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#34495e"),
    plot.caption = element_text(size = 9, hjust = 0.5, color = "#7f8c8d", face = "bold"),
    plot.margin = margin(10, 10, 10, 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
  
  return(venn_plot)
}

#' Create abundance violin plot for comprehensive panel
#' 
create_abundance_boxplot_for_panel <- function(data) {
  
  message("Creating abundance distribution violin plot by data source...")
  
  # Create violin plot of z-score normalized abundance by source
  violin_panel <- ggplot(data, aes(x = reorder(source, z_score, median), y = z_score, fill = technology)) +
    geom_violin(alpha = 0.8, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, alpha = 0.9, outlier.size = 0.3, outlier.alpha = 0.5, 
                 show.legend = FALSE, color = "black", fill = "white") +
    scale_fill_manual(values = get_plot_colors("technology"), name = "Technology") +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    ) +
    labs(
      title = "(H) Abundance Distribution by Database",
      subtitle = "Z-score normalized protein abundances (violin + boxplot)",
      x = "Database (ordered by median abundance)",
      y = "Z-Score Normalized Abundance"
    )
  
  return(violin_panel)
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
  
  # Create technology mapping for colors
  tech_mapping <- data %>%
    select(source, technology) %>%
    distinct() %>%
    deframe()
  
  # Create color mapping based on technology
  tech_colors <- get_plot_colors("technology")
  set_colors <- tech_colors[tech_mapping[names(gene_lists)]]
  
  # Calculate intersection sizes
  upset_data <- UpSetR::fromList(gene_lists)
  intersection_sizes <- colSums(upset_data)
  
  # Filter to keep only intersections with 48 or more proteins
  keep_intersections <- names(intersection_sizes[intersection_sizes >= 48])
  
  # Create UpSet plot using fromList approach
  upset_plot <- UpSetR::upset(
    UpSetR::fromList(gene_lists),
    nsets = length(gene_lists),
    sets = names(gene_lists),
    keep.order = TRUE,
    order.by = "freq",
    decreasing = TRUE,
    text.scale = 1,
    point.size = 3,
    line.size = 1.2,
    mainbar.y.label = "Number of Proteins",
    sets.x.label = "Total Proteins per Database",
    sets.bar.color = "#4575b4",  # Use the same blue as panel version
    main.bar.color = "#4575b4",  # Use the same blue as panel version
    matrix.color = "#4575b4",    # Use the same blue as panel version
    nintersects = NA  # Show all intersections that meet our threshold
  )
  
  # Save as PNG with larger dimensions to accommodate all intersections
  png(file.path(plot_dir, "04_upset_protein_overlap.png"), 
      width = 15, height = 8, units = "in", res = 300, bg = "white")  # Increased width to accommodate all bars
  print(upset_plot)
  dev.off()
  
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
      plot.title = element_text(size = 20, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      axis.text.y = element_text(size = 16),  # Size for the intersection counts
      axis.text.x = element_text(size = 16, face = "bold"),  # Size for other axis text
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      strip.text = element_text(size = 18, face = "bold"),  # Increased size for database names
      strip.background = element_blank()
    ) +
    labs(
      title = "(C) Intersection between different databases",
      x = "Database Combinations",
      y = "Number of Proteins"
    )
  
  return(panel_plot)
}

#' Create plasma databases dot plot for comprehensive panel
#' 
#' @param data Data frame with normalized protein data
#' @return ggplot object
#' 
create_plasma_databases_dot_plot_for_panel <- function(data) {
  
  # Get database colors, with PeptideAtlas in black as reference
  db_colors <- get_plot_colors("databases")
  db_colors["PeptideAtlas"] <- "black"
  
  # Create simplified version for panel
  p <- ggplot(data, aes(x = order, y = z_score, color = source)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_manual(values = db_colors, name = "Database") +
    labs(
      title = "(F) Protein abundance correlation with PeptideAtlas",
      subtitle = "Z-score normalized values ordered by PeptideAtlas",
      x = "Proteins (ordered)",
      y = "Z-Score"
    ) +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  return(p)
}

#' Generic density plot function
create_density_plot <- function(data, value_col, group_col, title) {
  ggplot(data, aes(x = .data[[value_col]], fill = .data[[group_col]])) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = get_plot_colors("databases"), name = "Data Source") +
    labs(title = title, x = NULL, y = "Density") +
    theme_blood_proteomics() +
    theme(legend.position = "right")
}



#' Create detailed log vs z-score normalization comparison
#' 
create_log_vs_zscore_comparison <- function(data, plot_dir) {
  
  message("Creating log vs z-score normalization comparison plots...")
  
  # 1. Histogram comparison - Log10 values
  p_hist_log <- ggplot(data, aes(x = log_abundance, fill = source)) +
    geom_histogram(alpha = 0.7, bins = 50) +
    facet_wrap(~source, scales = "free", ncol = 2, 
               labeller = labeller(source = function(x) paste0("(", letters[seq_along(x)], ") ", x))) +
    scale_fill_manual(values = get_plot_colors("databases")) +
    theme_blood_proteomics() +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 10)
    ) +
    labs(
      title = "Distribution of Protein Quantification Values (Log10)",
      subtitle = "Log10-transformed values across plasma databases",
      x = "Log10(Quantification Value)",
      y = "Number of Proteins"
    )
  
  save_plot_standard(p_hist_log, "plasma_proteins_quantification_distributions_log10", 
                    plot_dir, width = 12, height = 10)
  
  # 2. Histogram comparison - Z-score normalized
  p_hist_zscore <- ggplot(data, aes(x = z_score, fill = source)) +
    geom_histogram(alpha = 0.7, bins = 50) +
    facet_wrap(~source, scales = "free", ncol = 2,
               labeller = labeller(source = function(x) paste0("(", letters[seq_along(x)], ") ", x))) +
    scale_fill_manual(values = get_plot_colors("databases")) +
    theme_blood_proteomics() +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 10)
    ) +
    labs(
      title = "Distribution of Z-Score Normalized Protein Quantification Values",
      subtitle = "Z-score normalized values for direct cross-database comparison",
      x = "Z-Score (standardized within each database)",
      y = "Number of Proteins",
      caption = "Z-scores calculated within each database: (log10(value) - mean) / sd"
    )
  
  save_plot_standard(p_hist_zscore, "plasma_proteins_quantification_distributions_zscore", 
                    plot_dir, width = 12, height = 10)
  
  # 3. Density plots for direct comparison (log10)
  p_density_log <- ggplot(data, aes(x = log_abundance, fill = source)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = get_plot_colors("databases")) +
    theme_blood_proteomics() +
    labs(
      title = "Protein Quantification Value Distributions Comparison (Log10)",
      subtitle = "Density plots of log10-transformed values across plasma databases",
      x = "Log10(Quantification Value)",
      y = "Density",
      fill = "Database"
    )
  
  save_plot_standard(p_density_log, "plasma_proteins_quantification_density_log10", plot_dir)
  
  # 4. Z-score normalized density plots for direct comparison
  p_density_zscore <- ggplot(data, aes(x = z_score, fill = source)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = get_plot_colors("databases")) +
    theme_blood_proteomics() +
    labs(
      title = "Z-Score Normalized Protein Quantification Value Distributions",
      subtitle = "Z-score normalized density plots for direct cross-database comparison",
      x = "Z-Score (standardized within each database)",
      y = "Density",
      fill = "Database",
      caption = "Z-scores calculated within each database: (log10(value) - mean) / sd"
    )
  
  save_plot_standard(p_density_zscore, "plasma_proteins_quantification_density_zscore", plot_dir)
  
  # 5. Side-by-side comparison plot
  p_comparison <- p_density_log + p_density_zscore + 
    plot_layout(ncol = 1) +
    plot_annotation(
      title = "Plasma Protein Quantification: Log10 vs Z-Score Normalized Comparison",
      subtitle = "Comparing original log10 values (top) vs z-score normalized values (bottom)",
      tag_levels = list(c('(a)', '(b)')),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  save_plot_standard(p_comparison, "plasma_proteins_log10_vs_zscore_comparison", 
                    plot_dir, width = 12, height = 12)
  
  # 6. Statistical summary comparison
  create_normalization_statistics_table(data, plot_dir)
  
  message("✅ Log vs z-score comparison plots created successfully!")
}

#' Create statistical summary table comparing normalization methods
#' 
create_normalization_statistics_table <- function(data, plot_dir) {
  
  # Calculate statistics for both log and z-score values
  stats_summary <- data %>%
    group_by(source, technology) %>%
    summarise(
      n_proteins = n(),
      # Log10 statistics
      log10_mean = round(mean(log_abundance, na.rm = TRUE), 3),
      log10_median = round(median(log_abundance, na.rm = TRUE), 3),
      log10_sd = round(sd(log_abundance, na.rm = TRUE), 3),
      log10_min = round(min(log_abundance, na.rm = TRUE), 3),
      log10_max = round(max(log_abundance, na.rm = TRUE), 3),
      log10_range = round(log10_max - log10_min, 3),
      # Z-score statistics  
      zscore_mean = round(mean(z_score, na.rm = TRUE), 3),
      zscore_median = round(median(z_score, na.rm = TRUE), 3),
      zscore_sd = round(sd(z_score, na.rm = TRUE), 3),
      zscore_min = round(min(z_score, na.rm = TRUE), 3),
      zscore_max = round(max(z_score, na.rm = TRUE), 3),
      zscore_range = round(zscore_max - zscore_min, 3),
      .groups = "drop"
    ) %>%
    arrange(desc(n_proteins))
  
  # Save the statistics table
  output_dir <- get_output_path("01_plasma_protein_analysis", subdir = "tables")
  write_csv(stats_summary, file.path(output_dir, "plasma_proteins_normalization_statistics.csv"))
  
  message(sprintf("📊 Normalization statistics saved: %d sources compared", nrow(stats_summary)))
  
  return(stats_summary)
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
    geom_text(aes(label = scales::comma(count)), hjust = -0.1, size = 8.0) +
    coord_flip() +
    scale_fill_manual(values = get_plot_colors("technology"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 20, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_blank(),
      legend.position = "none",
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      axis.text.y = element_text(size = 16, face = "bold"),
      panel.grid.major.x = element_line(color = "grey90", size = 0.3),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
    ) +
    labs(
      title = "(A) Number of proteins by database",
      x = "Database",
      y = "Number of Proteins"
    )
  
  # Panel B: Technology Classification (Enhanced)
  tech_summary <- summary_stats %>%
    filter(!source %in% c("Total Across Sources", "MS Technologies Combined")) %>%
    group_by(technology) %>%
    summarise(
      databases = n(),
      total_proteins = sum(unique_genes),
      .groups = "drop"
    )
  
  panel_B <- ggplot(tech_summary, aes(x = technology, y = total_proteins, fill = technology)) +
    geom_col(alpha = 0.85, width = 0.6) +
    geom_text(aes(label = paste0(scales::comma(total_proteins), "\n(", databases, " DBs)")), 
              vjust = -0.3, size = 8.0, lineheight = 0.9) +
    scale_fill_manual(values = get_plot_colors("technology"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)), labels = scales::comma) +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 20, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_blank(),
      legend.position = "none",
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(size = 18, face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
    ) +
    labs(
      title = "(B) Number of proteins by Technology",
      x = "Technology Type",
      y = "Total Proteins"
    )
  
  # === SECTION 2: OVERLAP ANALYSIS ===
  
  # Panel C: UpSet Plot for Protein Overlap (Enhanced)
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
  
  panel_C <- create_upset_plot_for_panel(gene_lists_panel, set_colors_panel) +
    theme(
      plot.title = element_text(size = 20, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16)
    ) +
    labs(
      title = "(C) Intersection between different databases",
      x = "Database Combinations",
      y = "Number of Proteins"
    )
  
  # Panel E: Abundance Distribution (Enhanced)
  panel_E <- ggplot(normalized_data, aes(x = reorder(source, z_score, median), y = z_score, fill = technology)) +
    geom_violin(alpha = 0.8, scale = "width", trim = TRUE, width = 0.6) +
    geom_boxplot(width = 0.1, alpha = 0.9, outlier.size = 0.4, outlier.alpha = 0.6, 
                 show.legend = FALSE, color = "black", fill = "white") +
    stat_summary(fun = median, geom = "point", shape = 20, size = 2, color = "red", alpha = 0.8) +
    scale_fill_manual(values = get_plot_colors("technology"), guide = "none") +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 18, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
      legend.position = "none",
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3)
    ) +
    labs(
      title = "(D) Distributions of abundances by databases",
      subtitle = "z-score normalized",
      x = NULL,
      y = "Z-Score"
    )
  
  # Panel F: Cross-Database Dot Plot (Enhanced)
  peptideatlas_data <- normalized_data %>%
    filter(source == "PeptideAtlas") %>%
    arrange(z_score) %>%
    mutate(order = row_number()) %>%
    select(gene, order)
  
  dot_plot_data <- normalized_data %>%
    inner_join(peptideatlas_data, by = "gene") %>%
    group_by(gene, source) %>%
    summarise(
      z_score = median(z_score, na.rm = TRUE),
      order = first(order),
      .groups = "drop"
    )
  
  db_colors <- get_plot_colors("databases")
  db_colors["PeptideAtlas"] <- "#1a1a1a"  # Darker reference color
  
  panel_F <- ggplot(dot_plot_data, aes(x = order, y = z_score, color = source)) +
    geom_point(alpha = 0.7, size = 1.0) +
    scale_color_manual(values = db_colors, name = "Database") +
    scale_x_continuous(labels = scales::comma) +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 20, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_text(size = 16),  # Added styling for subtitle
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3)
    ) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 4, alpha = 0.9))) +
    labs(
      title = "(E) Protein abundance correlation with PeptideAtlas",
      subtitle = "Proteins ordered by PeptideAtlas z-score",
      x = NULL,  # Removed x-axis title
      y = "Z-Score"
    )
  
  # Panel G: Correlation Heatmap (Enhanced)
  correlation_data_panel <- normalized_data %>%
    select(gene, source, z_score, technology) %>%
    group_by(gene) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    group_by(gene, source) %>%
    summarise(
      z_score = median(z_score, na.rm = TRUE),
      technology = first(technology),
      .groups = "drop"
    )
  
  correlation_matrix_data_panel <- correlation_data_panel %>%
    select(gene, source, z_score) %>%
    pivot_wider(names_from = source, values_from = z_score) %>%
    column_to_rownames("gene")
  
  correlation_matrix_panel <- cor(correlation_matrix_data_panel, use = "pairwise.complete.obs")
  
  panel_G <- create_correlation_heatmap_for_panel(correlation_matrix_panel, normalized_data) +
    theme(
      plot.title = element_text(size = 20, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 18, face = "bold"),
      legend.position = "none",
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
    ) +
    labs(
      title = "(F) Cross-database abundance correlation",
      subtitle = NULL
    )
  
  # === COMBINE PANELS WITH ENHANCED LAYOUT ===
  
  # Row 1: Coverage (A-B)
  # Row 2: Overlap and Distribution (C-E) 
  # Row 3: Quantification (F-G)
  comprehensive_panel <- (panel_A | panel_B) / 
                         (panel_C | panel_E) /
                         (panel_F | panel_G) +
                         plot_layout(heights = c(1.1, 1.0, 1.3),
                                   widths = c(1.5, 0.5)) &
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
  
  # Save with enhanced specifications
  save_plot_standard(comprehensive_panel, "00_comprehensive_plasma_analysis_panel", plot_dir,
                    width = 30, height = 22)  # Increased height to accommodate larger correlation matrix
  
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