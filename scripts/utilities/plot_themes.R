# ============================================================================
# STANDARDIZED PLOT THEMES AND FUNCTIONS
# ============================================================================
# Description: Centralized plotting utilities for consistent visualizations
# Eliminates repetitive ggplot code and ensures consistent styling
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(patchwork)
})

# Source configuration
source("scripts/config/analysis_config.R")

#' Create standardized theme for blood proteomics plots
#' 
#' @param base_size Base font size
#' @param grid_lines Include major grid lines (default: TRUE)
#' @return ggplot2 theme object
#' 
theme_blood_proteomics <- function(base_size = NULL, grid_lines = TRUE) {
  
  if (is.null(base_size)) {
    base_size <- PLOT_CONFIG$theme$base_size
  }
  
  base_theme <- theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(
        size = PLOT_CONFIG$theme$title_size, 
        face = "bold", 
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size = PLOT_CONFIG$theme$subtitle_size, 
        hjust = 0.5,
        margin = margin(b = 15)
      ),
      axis.title = element_text(size = PLOT_CONFIG$theme$axis_title_size),
      axis.text = element_text(size = PLOT_CONFIG$theme$axis_text_size),
      legend.title = element_text(size = PLOT_CONFIG$theme$axis_title_size),
      legend.text = element_text(size = PLOT_CONFIG$theme$legend_text_size),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(size = PLOT_CONFIG$theme$axis_title_size, face = "bold")
    )
  
  if (!grid_lines) {
    base_theme <- base_theme +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  } else {
    base_theme <- base_theme +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()
      )
  }
  
  return(base_theme)
}

#' Create protein count bar plot
#' 
#' @param data Data frame with source and count columns
#' @param color_by Column to use for coloring (default: "technology")
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @return ggplot object
#' 
create_protein_count_plot <- function(data, color_by = "technology", title = NULL, subtitle = NULL) {
  
  colors <- get_plot_colors(color_by)
  
  p <- ggplot(data, aes(x = reorder(source, count), y = count, fill = .data[[color_by]])) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = comma(count)), hjust = -0.1, size = 3.5) +
    coord_flip() +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
    theme_blood_proteomics(grid_lines = FALSE) +
    labs(
      title = title %||% "Protein Counts by Data Source",
      subtitle = subtitle,
      x = "Data Source",
      y = "Number of Proteins",
      fill = str_to_title(color_by)
    )
  
  return(p)
}

#' Create violin plot with improved aesthetics
#' 
#' @param data Data frame with gene, abundance, and grouping columns
#' @param x_var Column for x-axis grouping
#' @param y_var Column for y-axis values
#' @param fill_var Column for fill colors
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param colors Named vector of colors (optional)
#' @return ggplot object
#' 
create_enhanced_violin_plot <- function(data, x_var, y_var, fill_var, 
                                       title = NULL, subtitle = NULL, colors = NULL) {
  
  if (is.null(colors)) {
    colors <- get_plot_colors("technology")
  }
  
  p <- ggplot(data, aes(x = reorder(.data[[x_var]], .data[[y_var]], FUN = median), 
                        y = .data[[y_var]], fill = .data[[fill_var]])) +
    # Enhanced violin plots
    geom_violin(alpha = 0.8, 
                scale = "area",
                trim = TRUE,
                adjust = 1.2,
                draw_quantiles = c(0.25, 0.5, 0.75),
                color = "gray30",
                size = 0.3) +
    # Refined boxplot overlay
    geom_boxplot(width = 0.15, 
                 alpha = 0.6,
                 outlier.size = 1.2,
                 outlier.alpha = 0.7,
                 outlier.color = "gray20",
                 color = "gray30",
                 size = 0.5) +
    # Median points for emphasis
    stat_summary(fun = median, 
                 geom = "point", 
                 shape = 21, 
                 size = 2.5, 
                 fill = "white", 
                 color = "gray20", 
                 stroke = 1) +
    coord_flip() +
    scale_fill_manual(values = colors) +
    theme_blood_proteomics() +
    labs(
      title = title,
      subtitle = subtitle,
      x = str_to_title(str_replace_all(x_var, "_", " ")),
      y = str_to_title(str_replace_all(y_var, "_", " ")),
      fill = str_to_title(str_replace_all(fill_var, "_", " "))
    )
  
  return(p)
}

#' Create biomarker distribution plot
#' 
#' @param data Data frame with protein abundance data
#' @param biomarker_genes Vector of biomarker gene symbols
#' @param abundance_col Column name for abundance values
#' @param source_col Column name for data source
#' @param title Plot title
#' @return ggplot object
#' 
create_biomarker_distribution_plot <- function(data, biomarker_genes, 
                                              abundance_col = "abundance",
                                              source_col = "source",
                                              title = NULL) {
  
  plot_data <- data %>%
    mutate(
      log_abundance = log10(.data[[abundance_col]] + PROJECT_CONFIG$analysis$log_transform_offset),
      is_biomarker = gene %in% biomarker_genes
    )
  
  biomarker_count <- sum(plot_data$is_biomarker)
  total_count <- nrow(plot_data)
  percentage <- round(biomarker_count / length(biomarker_genes) * 100, 1)
  
  colors <- get_plot_colors("biomarker")
  
  p <- ggplot(plot_data) +
    geom_density(aes(x = log_abundance, fill = is_biomarker), alpha = 0.6) +
    scale_fill_manual(
      values = colors,
      labels = c("Other proteins", "Biomarkers"),
      name = "Protein Type"
    ) +
    theme_blood_proteomics() +
    labs(
      title = title %||% "Protein Abundance Distribution",
      subtitle = sprintf("Biomarkers detected: %d (%.1f%%)", biomarker_count, percentage),
      x = expression(log[10]("Abundance + 1")),
      y = "Density"
    )
  
  return(p)
}

#' Create technology comparison plot
#' 
#' @param summary_data Data frame with technology summary statistics
#' @param title Plot title
#' @return ggplot object
#' 
create_technology_comparison_plot <- function(summary_data, title = NULL) {
  
  colors <- get_plot_colors("technology")
  
  p <- ggplot(summary_data, aes(x = reorder(technology, protein_count), 
                                y = protein_count, fill = technology)) +
    geom_col(alpha = 0.8, width = 0.7) +
    geom_text(aes(label = comma(protein_count)), hjust = -0.1, size = 4, fontface = "bold") +
    coord_flip() +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
    theme_blood_proteomics(grid_lines = FALSE) +
    labs(
      title = title %||% "Protein Detection by Technology",
      x = "Technology",
      y = "Number of Proteins",
      fill = "Technology"
    )
  
  return(p)
}

#' Apply z-score normalization to abundance data
#' 
#' @param data Data frame with abundance data
#' @param abundance_col Column name for abundance values
#' @param group_col Column name for grouping (e.g., "source")
#' @return Data frame with z_score column added
#' 
apply_zscore_normalization <- function(data, abundance_col = "abundance", group_col = "source") {
  
  data %>%
    mutate(log_abundance = log10(.data[[abundance_col]] + PROJECT_CONFIG$analysis$log_transform_offset)) %>%
    group_by(.data[[group_col]]) %>%
    mutate(
      z_score = (log_abundance - mean(log_abundance, na.rm = TRUE)) / sd(log_abundance, na.rm = TRUE)
    ) %>%
    ungroup()
}

#' Save plot with standardized settings
#' 
#' @param plot ggplot object
#' @param filename Output filename (without extension)
#' @param output_dir Output directory
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution in dots per inch (default: 600 for high-quality output)
#' @param device Output device (e.g., "tiff", "png", "pdf")
#' @note For publication-quality figures, use DPI 600 or higher
#' 
save_plot_standard <- function(plot, filename, output_dir, 
                              width = NULL, height = NULL, 
                              dpi = 600, device = "tiff") {
  
  if (is.null(width)) width <- PLOT_CONFIG$dimensions$default_width
  if (is.null(height)) height <- PLOT_CONFIG$dimensions$default_height
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set file extension based on device
  ext <- switch(device,
                "tiff" = "tiff",
                "png" = "png",
                "pdf" = "pdf",
                "tiff")  # default to tiff
  
  output_file <- file.path(output_dir, paste0(filename, ".", ext))
  
  # Save with specific device settings
  if (device == "tiff") {
    ggsave(
      filename = output_file,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      device = "tiff",
      compression = "lzw",  # Use LZW compression for smaller file size
      bg = "white"
    )
  } else {
    ggsave(
      filename = output_file,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      device = device,
      bg = "white"
    )
  }
  
  message(sprintf("Saved plot: %s", output_file))
}

#' Create ranked abundance plot (waterfall style)
#' 
#' @param data Data frame with gene and abundance columns
#' @param biomarker_genes Vector of biomarker gene symbols
#' @param source_name Name of the data source
#' @param abundance_col Column name for abundance values
#' @return ggplot object
#' 
create_ranked_abundance_plot <- function(data, biomarker_genes, source_name, 
                                        abundance_col = "abundance") {
  
  plot_data <- data %>%
    arrange(desc(.data[[abundance_col]])) %>%
    mutate(
      rank = row_number(),
      is_biomarker = gene %in% biomarker_genes,
      log_abundance = log10(.data[[abundance_col]] + PROJECT_CONFIG$analysis$log_transform_offset)
    )
  
  biomarker_subset <- plot_data[plot_data$is_biomarker, ]
  
  p <- ggplot(plot_data, aes(x = rank, y = log_abundance)) +
    # Background area
    geom_area(fill = "lightgrey", alpha = 0.7) +
    # Biomarker highlights
    geom_segment(data = biomarker_subset, 
                 aes(x = rank, xend = rank, y = 0, yend = log_abundance),
                 color = BIOMARKER_CONFIG$highlight_color, alpha = 0.8, size = 0.8) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = function(x) parse(text = paste0("10^", x))) +
    theme_blood_proteomics(grid_lines = FALSE) +
    labs(
      title = paste(source_name, "Protein Abundance Ranking"),
      subtitle = paste0("n = ", nrow(plot_data), " proteins | ", 
                       sum(plot_data$is_biomarker), " biomarkers highlighted"),
      x = "Protein Rank (highest to lowest)",
      y = expression(log[10]("Abundance"))
    )
  
  return(p)
} 