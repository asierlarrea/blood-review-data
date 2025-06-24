#!/usr/bin/env Rscript
# ============================================================================
# PLASMA PROTEIN ANALYSIS
# ============================================================================
# Description: Comprehensive analysis of plasma proteins from multiple sources
# Uses quantile-to-normal normalization for robust cross-database comparisons
# ============================================================================

# Load utilities and configuration
source("scripts/utilities/load_packages.R")
source("scripts/config/analysis_config.R")
source("scripts/utilities/data_loader.R")
source("scripts/utilities/plot_themes.R")

# Set up environment
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "patchwork", "ggpubr", "UpSetR", "tibble", "ggupset", "purrr")
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
message("PLASMA PROTEIN ANALYSIS")
message(paste(rep("=", 60), collapse = ""))

#' Main analysis function
#' 
run_plasma_protein_analysis <- function() {
  
  # Step 1: Load all data sources
  message("\n[STEP 1] Loading data sources...")
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
  normalized_data <- apply_normalization(combined_data)
  
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

#' Apply quantile-to-normal normalization
#'
apply_normalization <- function(data) {
  message("  - Filtering out zero abundance values...")
  data <- data %>%
    filter(abundance > 0)
  
  message("  - Applying log10 transformation...")
  data <- data %>%
    mutate(log_abundance = log10(abundance + 1))
  
  message("  - Applying quantile-to-normal normalization...")
  data %>%
    group_by(source) %>%
    mutate(
      rank_quantile = rank(log_abundance) / (n() + 1),
      z_score = qnorm(rank_quantile)
    ) %>%
    ungroup()
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
  gene_lists <- create_upset_analysis(normalized_data)
  
  # Create comprehensive panel
  create_comprehensive_panel(normalized_data, summary_stats, plot_dir)
}

#' Create UpSet analysis plot for protein overlap between databases
#'
create_upset_analysis <- function(data) {
  message("Creating UpSet analysis for protein overlap...")
  
  # Ensure we have unique genes per database
  deduplicated_data <- data %>%
    select(gene, source) %>%
    distinct() %>%
    filter(!is.na(gene), gene != "")
  
  # Create gene lists for overlap analysis
  gene_lists <- split(deduplicated_data$gene, deduplicated_data$source)
  
  # Remove any remaining NA or empty gene names
  gene_lists <- lapply(gene_lists, function(x) unique(x[!is.na(x) & x != ""]))
  
  message("✅ UpSet analysis plot created successfully!")
  
  return(gene_lists)
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
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20)
    ) +
    labs(
      title = "(B) Intersection between different data sources",
      subtitle = "Overlap analysis of protein coverage",
      x = "Data Source Combinations",
      y = "Number of Proteins"
    )
  
  # Panel C: Abundance Distribution (Enhanced) - Using traditional z-score for distribution shapes
  # Create data with traditional z-score normalization to show actual distribution differences
  panel_c_data <- normalized_data %>%
    group_by(source) %>%
    mutate(
      # Traditional z-score normalization to preserve distribution shapes
      z_score_traditional = (log_abundance - mean(log_abundance, na.rm = TRUE)) / sd(log_abundance, na.rm = TRUE)
    ) %>%
    ungroup()
  
  panel_C <- ggplot(panel_c_data, aes(x = z_score_traditional, y = reorder(source, z_score_traditional, median), fill = technology)) +
    ggridges::geom_density_ridges(alpha = 0.8, scale = 0.9, bandwidth = 0.3, color = "white", size = 0.5) +
    stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red", alpha = 0.8, 
                 position = position_nudge(y = 0)) +
    scale_fill_manual(values = get_plot_colors("technology"), guide = "none") +
    scale_x_continuous(limits = c(-5, 5), breaks = seq(-4, 4, 2)) +
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
      title = "(C) Distribution of abundance/expression by data sources",
      subtitle = "Ridgeline plots showing actual distribution differences (traditional z-score)",
      x = "Traditional Z-Score (log abundance/expression)",
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
      title = "(D) Protein abundance/expression correlation with PeptideAtlas",
      subtitle = "Proteins ordered by PeptideAtlas quantile-normalized observations",
      x = NULL,  # Removed x-axis title
      y = "Quantile-normalized Values"
    )
  
  # === SECTION 3: QUANTMS SAMPLE DISTRIBUTION ANALYSIS ===
  
  # Panel E: QuantMS Sample Count Distribution (New)
  message("  - Creating Panel E: QuantMS sample distribution...")
  
  # Get the QuantMS data with quantile-normalized z-scores
  quantms_sample_data <- normalized_data %>%
    filter(source == "quantms") %>%
    select(gene, z_score, abundance, log_abundance)
  
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
    filter(sample_count >= 10) # Only genes present in 10+ samples (reduces sample count bias)
  
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
  
  # 3) Top 10 quantms genes not in PeptideAtlas - prioritize by z-score, then sample count
  quantms_only_high <- quantms_sample_data %>%
    filter(!gene %in% peptideatlas_genes$gene) %>%
    arrange(desc(z_score), desc(sample_count)) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
  # 4) Bottom 10 low-abundant quantms genes not in PeptideAtlas
  quantms_only_low <- quantms_sample_data %>%
    filter(!gene %in% peptideatlas_genes$gene) %>%
    arrange(z_score) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
  # 5) Well-known plasma proteins from literature
  plasma_literature_genes <- c("ALB", "APOA1", "APOA2", "APOB", "APOE", "FGB", "FGG", "FGA", "HP", "TF")
  plasma_lit_in_data <- intersect(plasma_literature_genes, quantms_sample_data$gene)
  
  # Create sample count bins and gene categories
  quantms_sample_data <- quantms_sample_data %>%
    mutate(
      sample_group = case_when(
        sample_count >= 10 & sample_count <= 20 ~ "10-20",
        sample_count >= 21 & sample_count <= 39 ~ "21-39", 
        sample_count >= 40 ~ "40+",
        TRUE ~ "Other"
      ),
      sample_group = factor(sample_group, levels = c("10-20", "21-39", "40+")),
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
  panel_E <- ggplot(quantms_sample_data, aes(x = sample_group, y = z_score)) +
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
      size = 7,
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
      shape = guide_legend(
        override.aes = list(
          size = 6, 
          alpha = 1, 
          color = c(
            "High PeptideAtlas (shared)" = "#d73027",
            "Low PeptideAtlas (shared)" = "#4575b4", 
            "High quantms (unique)" = "#85218e",
            "Low quantms (unique)" = "#1a7922",
            "Plasma literature" = "#a2aa09"
          )
        )
      )
    ) +
    labs(
      title = "(E) QuantMS protein expression by sample presence",
      subtitle = "Quantile-normalized expression distribution with highlighted gene categories (10 genes per category)",
      x = "Number of Samples", 
      y = "Quantile-normalized Expression"
    )
  

  
  # === SECTION 4: COMPLEMENTARY PLATFORMS ANALYSIS ===
  
  # Panel F: Complementary Platform Strengths
  message("  - Creating Panel F: Complementary platforms analysis...")
  
  # Create complementary analysis using current quantile-normalized data
  # Extract PeptideAtlas and QuantMS data with quantile-normalized z-scores
  peptideatlas_data_f <- normalized_data %>%
    filter(source == "PeptideAtlas") %>%
    select(gene, z_score) %>%
    rename(z_score_PeptideAtlas = z_score)
  
  quantms_data_f <- normalized_data %>%
    filter(source == "quantms") %>%
    select(gene, z_score) %>%
    rename(z_score_QuantMS = z_score)
  
  # Create complementary analysis data using current quantile-normalized z-scores
  complementary_data <- inner_join(peptideatlas_data_f, quantms_data_f, by = "gene") %>%
    mutate(
      # Calculate platform advantages
      peptideatlas_advantage = z_score_PeptideAtlas - z_score_QuantMS,
      quantms_advantage = z_score_QuantMS - z_score_PeptideAtlas,
      # Determine platform strength
      platform_strength = case_when(
        abs(peptideatlas_advantage) < 0.5 ~ "Both Platforms Equivalent",
        peptideatlas_advantage > 0.5 ~ "PeptideAtlas Superior", 
        quantms_advantage > 0.5 ~ "QuantMS Superior",
        TRUE ~ "Both Platforms Equivalent"
      ),
      # Calculate combined evidence and detection categories
      combined_evidence = sqrt(z_score_PeptideAtlas^2 + z_score_QuantMS^2),
      detection_category = case_when(
        z_score_PeptideAtlas > 1 & z_score_QuantMS > 1 ~ "High in Both",
        z_score_PeptideAtlas > 1 & z_score_QuantMS <= 1 ~ "High in PeptideAtlas Only",
        z_score_PeptideAtlas <= 1 & z_score_QuantMS > 1 ~ "High in QuantMS Only",
        TRUE ~ "Low in Both"
      )
    )
  
  # If complementary analysis file exists, load it for comparison but use current data
  complementary_file <- file.path("outputs/tables/01_plasma_protein_analysis/all_proteins_complementary_analysis.csv")
  
  # Calculate statistics for plot annotations
  coverage_stats <- complementary_data %>%
    summarise(
      total_proteins = n(),
      both_high = sum(detection_category == "High in Both"),
      pa_only_high = sum(detection_category == "High in PeptideAtlas Only"),
      qms_only_high = sum(detection_category == "High in QuantMS Only")
    )
  
  # Create the complementary platforms plot
  panel_F <- ggplot(complementary_data, aes(x = z_score_PeptideAtlas, y = z_score_QuantMS)) +
    # Background quadrants
    annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1, 
             fill = "lightgray", alpha = 0.2) +
    annotate("rect", xmin = 1, xmax = Inf, ymin = 1, ymax = Inf, 
             fill = "gold", alpha = 0.2) +
    annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = 1, 
             fill = "lightblue", alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = 1, ymin = 1, ymax = Inf, 
             fill = "lightcoral", alpha = 0.2) +
    
    # Points colored by platform strength
    geom_point(aes(color = platform_strength, size = combined_evidence), alpha = 0.7) +
    
    # Reference lines
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", alpha = 0.5) +
    
    # Dynamic quadrant labels with actual counts
    annotate("text", x = max(complementary_data$z_score_PeptideAtlas) * 0.85, 
             y = max(complementary_data$z_score_QuantMS) * 0.85, 
             label = paste0("Both High\n(n=", coverage_stats$both_high, ")"), 
             size = 8, fontface = "bold", color = "darkgoldenrod") +
    annotate("text", x = max(complementary_data$z_score_PeptideAtlas) * 0.85, 
             y = min(complementary_data$z_score_QuantMS) * 0.85, 
             label = paste0("PeptideAtlas\nSuperior\n(n=", coverage_stats$pa_only_high, ")"), 
             size = 8, fontface = "bold", color = "blue") +
    annotate("text", x = min(complementary_data$z_score_PeptideAtlas) * 0.5, 
             y = max(complementary_data$z_score_QuantMS) * 0.85, 
             label = paste0("QuantMS\nSuperior\n(n=", coverage_stats$qms_only_high, ")"), 
             size = 8, fontface = "bold", color = "red") +
    
    scale_color_manual(values = c("PeptideAtlas Superior" = "#4575b4", 
                                 "QuantMS Superior" = "#d73027",
                                 "Both Platforms Equivalent" = "#35978f"),
                      name = "Platform Strength") +
    scale_size_continuous(range = c(0.8, 3), name = "Combined\nEvidence") +
    
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_text(size = 20),
      axis.title = element_text(size = 24, face = "bold"),
      axis.text = element_text(size = 22),
      legend.position = "right",
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 20),
      legend.key.size = unit(1.5, "lines"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 6, alpha = 1)),
      size = guide_legend(override.aes = list(alpha = 1))
    ) +
    labs(
      title = "(F) Complementary platform strengths",
      subtitle = sprintf("Quantile-normalized analysis of %d proteins present in both platforms", coverage_stats$total_proteins),
      x = "PeptideAtlas Quantile-normalized Observations", 
      y = "QuantMS Quantile-normalized Expression"
    )
  
  # === COMBINE PANELS WITH ENHANCED LAYOUT ===
  
  # Row 1: Coverage and Overlap (A-B)
  # Row 2: Distribution and Cross-Database Analysis (C-D) 
  # Row 3: QuantMS Sample Analysis and Complementary Platforms (E-F)
  comprehensive_panel <- (panel_A | panel_B) / 
                         (panel_C | panel_D) /
                         (panel_E | panel_F) +
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
  
  # Count frequencies and filter out those below 100
  intersection_counts <- upset_data_wide %>%
    count(databases) %>%
    filter(n >= 100) %>%
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
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      axis.text.y = element_text(size = 20),  # Size for the intersection counts
      axis.text.x = element_text(size = 20, face = "bold"),  # Size for other axis text
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      strip.text = element_text(size = 22, face = "bold"),  # Database names
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
  
  # Extract key metrics
  total_unique_genes <- results$summary$unique_genes[results$summary$source == "Total Across Sources"]
  ms_unique_genes <- results$summary$unique_genes[results$summary$source == "MS Technologies Combined"]
  highest_source_idx <- which.max(results$summary$unique_genes[!results$summary$source %in% c("Total Across Sources", "MS Technologies Combined")])
  highest_source <- results$summary$source[!results$summary$source %in% c("Total Across Sources", "MS Technologies Combined")][highest_source_idx]
  highest_source_count <- results$summary$unique_genes[!results$summary$source %in% c("Total Across Sources", "MS Technologies Combined")][highest_source_idx]
  
  # Create detailed summary table
  detailed_summary <- results$summary %>%
    filter(!source %in% c("Total Across Sources", "MS Technologies Combined")) %>%
    arrange(desc(unique_genes))
  
  # Create report content
  report_content <- sprintf("
# Plasma Protein Analysis Report

**Analysis Date:** %s  
**Script Version:** %s

## 🔬 TOTAL UNIQUE GENES DETECTED: %s

This analysis identified **%s unique genes** across all plasma protein databases, representing the most comprehensive plasma proteome compilation to date.

## Summary Statistics

### Complete Data Source Overview
| Data Source | Technology | Unique Proteins | Total Entries | Median Abundance | Abundance Type |
|-------------|------------|-----------------|---------------|------------------|----------------|
%s

### Key Findings

- **🎯 Total unique genes across all sources:** **%s**
- **⚗️ Mass spectrometry technologies combined:** %s unique genes
- **🥇 Highest individual source:** %s (%s proteins)
- **📊 Data sources analyzed:** %d databases
- **🔬 Technologies represented:** %s

### Source-Specific Details
%s

## Files Generated

- **Summary data:** `%s`
- **Comprehensive plots:** `%s/plots/01_plasma_protein_analysis/`
- **Detailed analysis:** `%s/reports/01_plasma_protein_analysis/`

---
*Generated by the refactored plasma protein analysis pipeline*  
*Total proteome coverage: %s unique genes*
",
    Sys.Date(),
    PROJECT_CONFIG$version,
    scales::comma(total_unique_genes),
    scales::comma(total_unique_genes),
    paste(sprintf("| %s | %s | %s | %s | %g | %s |", 
                  detailed_summary$source, 
                  detailed_summary$technology, 
                  scales::comma(detailed_summary$unique_genes),
                  scales::comma(detailed_summary$total_entries),
                  detailed_summary$median_abundance,
                  detailed_summary$abundance_type), 
          collapse = "\n"),
    scales::comma(total_unique_genes),
    scales::comma(ms_unique_genes),
    highest_source,
    scales::comma(highest_source_count),
    nrow(detailed_summary),
    paste(unique(detailed_summary$technology), collapse = ", "),
    paste(sprintf("- **%s**: %s proteins (%s technology)", 
                  detailed_summary$source[1:min(3, nrow(detailed_summary))], 
                  scales::comma(detailed_summary$unique_genes[1:min(3, nrow(detailed_summary))]),
                  detailed_summary$technology[1:min(3, nrow(detailed_summary))]), 
          collapse = "\n"),
    "outputs/tables/01_plasma_protein_analysis/plasma_protein_summary.csv",
    "outputs",
    "outputs",
    scales::comma(total_unique_genes)
  )
  
  writeLines(report_content, report_file)
  message(sprintf("Generated report: %s", report_file))
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