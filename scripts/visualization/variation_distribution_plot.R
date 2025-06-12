#!/usr/bin/env Rscript

# Improved Variation Distribution Plot
# Comparing PEA and MS/MS variability measures with enhanced aesthetics

# Load required packages
source("scripts/utilities/load_packages.R")
required_packages <- c("ggplot2", "dplyr", "tidyr", "scales", "viridis", "RColorBrewer")
load_packages(required_packages)

# Create output directory
output_dir <- "outputs/plots/01_Database_Analysis/PEA_MSMS_Variability"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load the merged data from previous analysis
merged_data_file <- file.path(output_dir, "merged_pea_msms_data.csv")

if (file.exists(merged_data_file)) {
  cat("Loading merged data from previous analysis...\n")
  merged_data <- read.csv(merged_data_file, stringsAsFactors = FALSE)
  cat("Loaded data for", nrow(merged_data), "proteins\n")
  
  # Prepare data for distribution plot
  variation_long <- merged_data %>%
    select(Gene_Clean, Variation.between.individuals, Variation.within.individuals, CV_MSMS_mean) %>%
    pivot_longer(cols = c(Variation.between.individuals, Variation.within.individuals, CV_MSMS_mean),
                 names_to = "Variation_Type", values_to = "Variation_Value") %>%
    mutate(
      Variation_Type = case_when(
        Variation_Type == "Variation.between.individuals" ~ "PEA\nBetween Individuals",
        Variation_Type == "Variation.within.individuals" ~ "PEA\nWithin Individuals", 
        Variation_Type == "CV_MSMS_mean" ~ "MS/MS\nCoefficient of Variation"
      )
    ) %>%
    filter(!is.na(Variation_Value) & Variation_Value > 0)
  
  # Define beautiful color palette
  colors <- c("#1f77b4", "#ff7f0e", "#2ca02c")  # Professional blue, orange, green
  
  # Create enhanced distribution plot
  p <- ggplot(variation_long, aes(x = Variation_Type, y = Variation_Value, fill = Variation_Type)) +
    # Add violin plots with better transparency
    geom_violin(alpha = 0.8, trim = FALSE, scale = "width", color = "white", linewidth = 0.5) +
    
    # Add box plots with outliers
    geom_boxplot(width = 0.15, alpha = 0.9, outlier.size = 0.8, outlier.alpha = 0.6,
                 color = "white", linewidth = 0.4) +
    
    # Add median points
    stat_summary(fun = median, geom = "point", size = 3, color = "white", shape = 21, stroke = 1.5) +
    
    # Use manual colors
    scale_fill_manual(values = colors) +
    
    # Use log scale for better visualization
    scale_y_log10(
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      breaks = scales::trans_breaks("log10", function(x) 10^x, n = 6)
    ) +
    
    # Enhanced labels
    labs(
      title = "Distribution of Protein Variability Measures",
      subtitle = paste("Comparison across", length(unique(variation_long$Gene_Clean)), "proteins with both PEA and MS/MS data"),
      x = "Measurement Method",
      y = expression("Variation Value (log"[10]*" scale)"),
      caption = "Violin plots show data distribution; box plots show quartiles; white dots indicate medians"
    ) +
    
    # Clean theme with white background
    theme_minimal() +
    theme(
      # Background and panel
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      
      # Text elements
      text = element_text(color = "grey20"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15), color = "grey40"),
      plot.caption = element_text(size = 9, hjust = 0.5, margin = margin(t = 10), color = "grey50"),
      
      # Axis elements
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
      axis.text.x = element_text(size = 11, face = "bold", color = "grey30"),
      axis.text.y = element_text(size = 10, color = "grey40"),
      axis.line = element_line(color = "grey80", linewidth = 0.5),
      axis.ticks = element_line(color = "grey80", linewidth = 0.3),
      
      # Legend
      legend.position = "none",
      
      # Margins
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Add summary statistics as text annotations
  summary_stats <- variation_long %>%
    group_by(Variation_Type) %>%
    summarise(
      median_val = median(Variation_Value, na.rm = TRUE),
      mean_val = mean(Variation_Value, na.rm = TRUE),
      q25 = quantile(Variation_Value, 0.25, na.rm = TRUE),
      q75 = quantile(Variation_Value, 0.75, na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    )
  
  # Save the improved plot
  ggsave(file.path(output_dir, "variation_distributions_improved.png"), p, 
         width = 10, height = 8, dpi = 300, bg = "white")
  
  # Also save as PDF for publication quality
  ggsave(file.path(output_dir, "variation_distributions_improved.pdf"), p, 
         width = 10, height = 8, bg = "white")
  
  # Create a version with summary statistics table
  stats_table <- summary_stats %>%
    mutate(
      Type = gsub("\n", " ", Variation_Type),
      Median = round(median_val, 3),
      Mean = round(mean_val, 3),
      IQR = paste0("[", round(q25, 3), " - ", round(q75, 3), "]"),
      N = n
    ) %>%
    select(Type, N, Median, Mean, IQR)
  
  # Save summary statistics
  write.csv(stats_table, file.path(output_dir, "variation_distribution_stats.csv"), row.names = FALSE)
  
  # Print summary to console
  cat("\nSummary Statistics:\n")
  cat("==================\n")
  for (i in 1:nrow(stats_table)) {
    cat(sprintf("%-25s: N=%3d, Median=%.3f, Mean=%.3f, IQR=[%.3f-%.3f]\n",
                stats_table$Type[i], stats_table$N[i], 
                summary_stats$median_val[i], summary_stats$mean_val[i],
                summary_stats$q25[i], summary_stats$q75[i]))
  }
  
  cat("\nFiles generated:\n")
  cat("- variation_distributions_improved.png (high-resolution PNG)\n")
  cat("- variation_distributions_improved.pdf (publication-quality PDF)\n")
  cat("- variation_distribution_stats.csv (summary statistics)\n")
  cat("\nPlot saved with enhanced aesthetics and clean white background!\n")
  
} else {
  cat("Error: Merged data file not found. Please run the correlation analysis script first.\n")
  cat("Expected file:", merged_data_file, "\n")
} 
