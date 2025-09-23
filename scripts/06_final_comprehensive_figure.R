#!/usr/bin/env Rscript
# ============================================================================
# FINAL COMPREHENSIVE FIGURE
# ============================================================================
# Description: Creates a final comprehensive figure combining specific panels
# from plasma analysis (B, D, E) and biomarker analysis (A) into a new layout:
# Left column: A, B, C (from plasma analysis B, D, E)
# Right side: D (biomarker heatmap from biomarker analysis A)
# ============================================================================

# Load utilities and configuration
source("scripts/utilities/load_packages.R")
source("scripts/config/analysis_config.R")
source("scripts/utilities/data_loader.R")
source("scripts/utilities/plot_themes.R")

# Set up environment
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "patchwork", "ggpubr", "UpSetR", "tibble", "ggupset", "purrr", "ggrepel", "grid", "png")
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
message("FINAL COMPREHENSIVE FIGURE GENERATION")
message(paste(rep("=", 60), collapse = ""))

#' Create final comprehensive figure
#' 
create_final_comprehensive_figure <- function() {
  
  # Step 1: Load all data sources (same as plasma analysis)
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
  
  # Step 2: Combine and normalize data (same as plasma analysis)
  message("\n[STEP 2] Combining and normalizing data...")
  combined_data <- combine_data_sources(data_list)
  normalized_data <- apply_normalization(combined_data)
  
  # Step 3: Load biomarker data
  message("\n[STEP 3] Loading biomarker data...")
  biomarker_file <- "data/metadata/biomarkers_list.csv"
  biomarkers <- read_csv(biomarker_file, show_col_types = FALSE)
  biomarker_genes <- unique(biomarkers$gene_name)
  
  # Step 4: Create individual panels
  message("\n[STEP 4] Creating individual panels...")
  
  # Panel A: UpSet Plot (from original Panel B)
  panel_A <- create_upset_panel(normalized_data)
  
  # Panel B: Cross-Database Dot Plot (from original Panel D)
  panel_B <- create_dot_plot_panel(normalized_data)
  
  # Panel C: QuantMS Sample Distribution (from original Panel E)
  panel_C <- create_quantms_panel(normalized_data)
  
  # Panel D: Biomarker Heatmap (from biomarker analysis Panel A)
  panel_D <- create_biomarker_heatmap(normalized_data, biomarker_genes)
  
  # Step 5: Combine panels in the new layout
  message("\n[STEP 5] Combining panels in final layout...")
  final_figure <- combine_final_panels(panel_A, panel_B, panel_C, panel_D)
  
  # Step 6: Save the figure
  message("\n[STEP 6] Saving final figure...")
  save_final_figure(final_figure)
  
  message("\n✅ Final comprehensive figure created successfully!")
  
  return(final_figure)
}

#' Apply quantile-to-normal normalization (same as plasma analysis)
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

#' Create UpSet panel (Panel A - from original Panel B)
#'
create_upset_panel <- function(normalized_data) {
  message("  - Creating Panel A: UpSet plot...")
  
  # Create UpSet data
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
  
  panel_A <- create_upset_plot_for_panel(gene_lists_panel, set_colors_panel) +
    theme(
      plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_text(size = 18),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      plot.margin = margin(20, 20, 30, 20)  # Add bottom margin for spacing
    ) +
    labs(
      title = "(A) Intersection between different data sources",
      x = "Data Source Combinations",
      y = "Number of Proteins"
    )
  
  return(panel_A)
}

#' Create dot plot panel (Panel B - from original Panel D)
#'
create_dot_plot_panel <- function(normalized_data) {
  message("  - Creating Panel B: Cross-database dot plot...")
  
  # Create dot plot data
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
  
  panel_B <- ggplot(dot_plot_data, aes(x = order, y = z_score, color = source)) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_color_manual(values = db_colors, name = "Data Source") +
    scale_x_continuous(labels = scales::comma) +
    theme_blood_proteomics() +
    theme(
      plot.title = element_text(size = 28, face = "bold", color = "#2c3e50"),
      plot.subtitle = element_text(size = 20),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 20),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      plot.margin = margin(30, 20, 30, 20)  # Add top and bottom margins for spacing
    ) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 6, alpha = 0.9))) +
    labs(
      title = "(B) Protein abundance/detection frequency correlation with PeptideAtlas",
      x = NULL,
      y = "Quantile-normalized Values"
    )
  
  return(panel_B)
}

#' Create QuantMS panel (Panel C - from original Panel E)
#'
create_quantms_panel <- function(normalized_data) {
  message("  - Creating Panel C: QuantMS sample distribution...")
  
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
    filter(sample_count >= 10) # Only genes present in 10+ samples
  
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
  top_peptideatlas_in_quantms <- peptideatlas_genes %>%
    filter(gene %in% quantms_sample_data$gene) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
  low_peptideatlas_in_quantms <- peptideatlas_genes %>%
    filter(gene %in% quantms_sample_data$gene) %>%
    slice_tail(n = 10) %>%
    pull(gene)
  
  quantms_only_high <- quantms_sample_data %>%
    filter(!gene %in% peptideatlas_genes$gene) %>%
    arrange(desc(z_score), desc(sample_count)) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
  quantms_only_low <- quantms_sample_data %>%
    filter(!gene %in% peptideatlas_genes$gene) %>%
    arrange(z_score) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
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
  
  # Define colors and shapes for gene categories
  category_colors <- c(
    "High PeptideAtlas (shared)" = "#d73027",
    "Low PeptideAtlas (shared)" = "#4575b4",  
    "High quantms (unique)" = "#85218e",
    "Low quantms (unique)" = "#1a7922",
    "Plasma literature" = "#a2aa09",
    "Other" = "grey40"
  )
  
  category_shapes <- c(
    "High PeptideAtlas (shared)" = 16,
    "Low PeptideAtlas (shared)" = 17,
    "High quantms (unique)" = 15,
    "Low quantms (unique)" = 18,
    "Plasma literature" = 8,
    "Other" = 20
  )
  
  # Create violin plot with highlighted gene categories
  panel_C <- ggplot(quantms_sample_data, aes(x = sample_group, y = z_score)) +
    geom_violin(alpha = 0.6, scale = "width", trim = TRUE, width = 0.7, fill = "lightgrey", color = "grey60") +
    geom_jitter(data = filter(quantms_sample_data, gene_category == "Other"),
                alpha = 0.3, size = 2.5, width = 0.15, color = "grey50") +
    geom_jitter(data = filter(quantms_sample_data, gene_category != "Other"),
                aes(color = gene_category, shape = gene_category), 
                alpha = 0.9, size = 4.0, width = 0.15, seed = 42,
                show.legend = c(color = FALSE, shape = TRUE)) +
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
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      legend.position = "bottom",
      legend.title = element_text(size = 24, face = "bold"),
      legend.text = element_text(size = 22),
      legend.key.size = unit(2.0, "lines"),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      plot.margin = margin(30, 20, 20, 20)  # Add top margin for spacing
    ) +
    guides(
      shape = guide_legend(
        nrow = 2,
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
      title = "(C) quantms protein abundance by sample presence",
      x = "Number of Samples", 
      y = "Quantile-normalized Abundance"
    )
  
  return(panel_C)
}

#' Create biomarker heatmap panel (Panel D - from biomarker analysis Panel A)
#'
create_biomarker_heatmap <- function(normalized_data, biomarker_genes) {
  message("  - Creating Panel D: Biomarker heatmap...")
  
  # Get unique sources
  unique_sources <- unique(normalized_data$source)
  
  # Collect z-scores from all databases
  zscore_data <- list()
  for (source_name in unique_sources) {
    source_data <- normalized_data %>% filter(source == source_name)
    if (nrow(source_data) > 0) {
      zscore_data[[source_name]] <- source_data %>%
        select(gene, z_score)
    }
  }
  
  # Create a matrix of z-scores for biomarkers
  biomarker_matrix <- matrix(NA, 
                           nrow = length(biomarker_genes), 
                           ncol = length(zscore_data),
                           dimnames = list(biomarker_genes, names(zscore_data)))
  
  # Fill the matrix with z-scores
  for (db in names(zscore_data)) {
    db_data <- zscore_data[[db]]
    biomarker_matrix[, db] <- db_data$z_score[match(biomarker_genes, db_data$gene)]
  }
  
  # Calculate detection frequency for each biomarker
  detection_freq <- rowSums(!is.na(biomarker_matrix)) / ncol(biomarker_matrix) * 100
  
  # Calculate mean z-score (when detected)
  mean_zscore <- rowMeans(biomarker_matrix, na.rm = TRUE)
  
  # Create data frame for plotting
  plot_data <- as.data.frame(biomarker_matrix) %>%
    tibble::rownames_to_column("Biomarker") %>%
    tidyr::pivot_longer(-Biomarker, names_to = "Database", values_to = "Z_score") %>%
    mutate(
      Database = factor(Database, levels = names(zscore_data)),
      Biomarker = factor(Biomarker, levels = biomarker_genes[order(mean_zscore, decreasing = TRUE)]),
      Detection = ifelse(is.na(Z_score), "Not detected", "Detected"),
      Z_score = ifelse(is.na(Z_score), 0, Z_score),
      Z_score_capped = pmin(pmax(Z_score, -3), 3)
    )
  
  # Create the profile matrix plot
  panel_D <- ggplot(plot_data, aes(x = Database, y = Biomarker)) +
    geom_tile(aes(fill = Z_score_capped), color = "white", linewidth = 0.25, width = 0.65, height = 0.85, alpha = 1.0) +
    geom_point(data = filter(plot_data, Detection == "Not detected"),
               color = "grey80", size = 0.75) +
    scale_fill_gradient2(
      low = "#2166AC", 
      mid = "#E8E8E8",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-3, 3),
      name = "Z-score",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    geom_text(
      data = unique(plot_data[c("Biomarker")]) %>%
        mutate(
          Database = "Detection\nFrequency",
          label = sprintf("%.0f%%", detection_freq[as.character(Biomarker)])
        ),
      aes(label = label),
      hjust = 0.5,
      size = 9
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
      axis.text.y = element_text(size = 22, face = "bold"),  # Increased from 18 to 22
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 20, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 18),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.key.size = unit(2.0, "cm"),
      panel.spacing = unit(1, "lines"),
      plot.margin = margin(20, 40, 20, 20)  # Add more right margin to separate from edge
    ) +
    labs(
      title = "(D) Biomarker Detection Profile Across Databases"
    )
  
  return(panel_D)
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
      axis.text.y = element_text(size = 20),
      axis.text.x = element_text(size = 20, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      strip.text = element_text(size = 22, face = "bold"),
      strip.background = element_blank()
    ) +
    labs(
      title = "(A) Intersection between different data sources",
      x = "Data Source Combinations",
      y = "Number of Proteins"
    )
  
  return(panel_plot)
}

#' Combine panels in the final layout
#'
combine_final_panels <- function(panel_A, panel_B, panel_C, panel_D) {
  message("  - Combining panels in final layout...")
  
  # Left column: A, B, C (stacked vertically with increased spacing)
  left_column <- panel_A / panel_B / panel_C
  
  # Create an empty spacer
  spacer <- ggplot() + theme_void()
  
  # Right side: D (heatmap) - give it more space
  right_side <- panel_D
  
  # Combine left, spacer, and right with adjusted widths and heights
  final_figure <- left_column | spacer | right_side
  
  # Add overall styling with increased spacing and height adjustments
  final_figure <- final_figure +
    plot_annotation(
      theme = theme(
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        plot.margin = margin(20, 30, 15, 20)  # Increased right margin for better separation
      )
    ) +
    plot_layout(
      widths = c(1, 0.1, 1.2),  # Reduce Panel D width from 1.8 to 1.4
      heights = c(1, 1, 1, 1.2)  # Give Panel D 20% more height
    )
  
  return(final_figure)
}

#' Save the final figure
#'
save_final_figure <- function(final_figure) {
  output_dir <- "outputs/plots/06_final_plots"
  
  # Save PNG version
  ggsave(
    file.path(output_dir, "final_figure_comprehensive_figure.png"),
    final_figure,
    width = 36,
    height = 24,
    dpi = 600,
    bg = "white",
    device = "png"
  )
  
  # Save TIFF version
  ggsave(
    file.path(output_dir, "final_figure_comprehensive_figure.tiff"),
    final_figure,
    width = 36,
    height = 24,
    dpi = 600,
    bg = "white",
    device = "tiff",
    compression = "lzw"
  )
  
  message(sprintf("✅ Final figure saved to: %s", output_dir))
}

# Execute the analysis
tryCatch({
  
  final_figure <- create_final_comprehensive_figure()
  
}, error = function(e) {
  message(sprintf("❌ Analysis failed: %s", e$message))
  quit(status = 1)
})

message("\n🎉 Final comprehensive figure generation completed successfully!")
