#!/usr/bin/env Rscript

# Improved Ranked Abundance Plots for Plasma Protein Concentrations
# Creates waterfall-style plots showing protein concentrations ranked from highest to lowest
# with biomarker genes highlighted in red
# 
# IMPROVEMENTS:
# - Fixed biomarker counting bug (unique genes per source)
# - Updated deprecated ggplot2 aesthetics
# - Better error handling and data validation
# - Improved log transformation logic
# - Enhanced reporting and statistics
# - Better dynamic range calculation

# Load required packages (install via install_dependencies.R if needed)
source("load_packages.R")
required_packages <- c("ggplot2", "dplyr", "scales", "gridExtra", "viridis")
load_packages(required_packages)

# Create output directory
output_dir <- "plots/05_Ranked_Abundance_Improved"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to load biomarker genes with validation
load_biomarker_genes <- function() {
  if(!file.exists("biomarkers_list.tsv")) {
    stop("biomarkers_list.tsv not found!")
  }
  
  biomarkers <- read.delim("biomarkers_list.tsv", stringsAsFactors = FALSE)
  
  # Validate biomarkers file
  if(!"To" %in% colnames(biomarkers)) {
    stop("biomarkers_list.tsv must have a 'To' column with gene symbols")
  }
  
  biomarker_genes <- unique(toupper(biomarkers$To))
  biomarker_genes <- biomarker_genes[!is.na(biomarker_genes) & biomarker_genes != ""]
  
  message(paste("Loaded", length(biomarker_genes), "unique biomarker genes"))
  return(biomarker_genes)
}

# Function to calculate safe log transformation with appropriate offset
safe_log_transform <- function(abundance, abundance_type) {
  # Use appropriate offset based on abundance type and data range
  min_val <- min(abundance[abundance > 0], na.rm = TRUE)
  
  # Handle vector abundance_type by taking the first value
  abundance_type_single <- abundance_type[1]
  
  if(abundance_type_single == "mg/L") {
    offset <- min(0.001, min_val / 10)
  } else if(abundance_type_single == "ppm") {
    offset <- min(0.0001, min_val / 10)
  } else {
    offset <- 1
  }
  
  return(log10(abundance + offset))
}

# Function to validate and clean source data
validate_source_data <- function(data, source_name) {
  original_rows <- nrow(data)
  
  # Remove rows with missing or invalid abundance
  data <- data[!is.na(data$abundance) & data$abundance > 0, ]
  
  # Remove rows with missing gene symbols
  data <- data[!is.na(data$gene_symbol) & data$gene_symbol != "", ]
  
  filtered_rows <- nrow(data)
  
  if(filtered_rows < original_rows) {
    message(paste("  ", source_name, ": Filtered", 
                  original_rows - filtered_rows, "invalid rows"))
  }
  
  if(filtered_rows == 0) {
    warning(paste("No valid data remaining for", source_name))
    return(NULL)
  }
  
  return(data)
}

# Function to create ranked abundance plot for a single source
create_ranked_plot <- function(data, source_name, biomarker_genes, abundance_type) {
  # Validate and clean data
  data <- validate_source_data(data, source_name)
  if(is.null(data)) return(NULL)
  
  # Remove duplicates and rank by abundance
  source_data <- data %>%
    arrange(desc(abundance)) %>%
    mutate(
      rank = row_number(),
      is_biomarker = gene_symbol %in% biomarker_genes,
      log_abundance = safe_log_transform(abundance, abundance_type)
    )
  
  # Get biomarker subset for highlighting
  biomarker_subset <- source_data[source_data$is_biomarker, ]
  
  # Calculate statistics for subtitle
  n_unique_biomarkers <- length(unique(biomarker_subset$gene_symbol))
  total_biomarkers <- length(biomarker_genes)
  
  # Create the plot
  p <- ggplot(source_data, aes(x = rank, y = log_abundance)) +
    # Background proteins as grey area
    geom_area(fill = "lightgrey", alpha = 0.7) +
    # Add horizontal grid lines
    geom_hline(yintercept = seq(floor(min(source_data$log_abundance)), 
                               ceiling(max(source_data$log_abundance)), 1), 
               color = "white", alpha = 0.8, linewidth = 0.3) +
    # Highlight biomarkers as red vertical lines (FIXED: use linewidth instead of size)
    geom_segment(data = biomarker_subset, 
                aes(x = rank, xend = rank, y = 0, yend = log_abundance),
                color = "red", alpha = 0.8, linewidth = 0.8) +
    # Formatting
    scale_y_continuous(
      name = paste0("Concentration (", abundance_type, ", log10 scale)"),
      labels = function(x) {
        # Simplified label function
        parse(text = paste0("10^", round(x, 1)))
      }
    ) +
    scale_x_continuous(
      name = "Protein Rank (highest to lowest)",
      labels = scales::comma
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = paste0(source_name, " Plasma Protein Concentrations"),
      subtitle = paste0("n = ", nrow(source_data), " proteins | ", 
                       n_unique_biomarkers, "/", total_biomarkers, 
                       " biomarkers (", round(100*n_unique_biomarkers/total_biomarkers, 1), 
                       "%) highlighted in red")
    )
  
  return(p)
}

# Function to calculate safe dynamic range
calculate_dynamic_range <- function(abundance) {
  valid_abundance <- abundance[!is.na(abundance) & abundance > 0]
  
  if(length(valid_abundance) < 2) return(NA)
  
  min_val <- min(valid_abundance)
  max_val <- max(valid_abundance)
  
  # Avoid log(0) issues
  if(min_val <= 0) {
    min_val <- min(valid_abundance[valid_abundance > 0])
  }
  
  return(round(log10(max_val / min_val), 2))
}

# Function to create improved summary statistics
create_improved_summary <- function(plasma_data, biomarker_genes) {
  summary_stats <- plasma_data %>%
    group_by(source) %>%
    summarise(
      n_proteins = n_distinct(gene_symbol),
      # FIX: Count unique biomarkers per source, not total occurrences
      unique_biomarkers = n_distinct(gene_symbol[gene_symbol %in% biomarker_genes]),
      total_biomarkers = length(biomarker_genes),
      biomarker_coverage_pct = round(100 * unique_biomarkers / length(biomarker_genes), 1),
      min_abundance = min(abundance, na.rm = TRUE),
      max_abundance = max(abundance, na.rm = TRUE),
      median_abundance = median(abundance, na.rm = TRUE),
      mean_abundance = mean(abundance, na.rm = TRUE),
      sd_abundance = sd(abundance, na.rm = TRUE),
      dynamic_range = calculate_dynamic_range(abundance),
      abundance_type = first(abundance_type),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_proteins))
  
  return(summary_stats)
}

# Function to create combined plots with adaptive layout
create_all_ranked_plots <- function() {
  message("Loading plasma expression data...")
  
  # Load plasma data from previous analysis
  if(!file.exists("plots/04_Biomarker_Analysis/plasma_expression_data.csv")) {
    stop("Plasma expression data not found. Please run biomarker_plasma_analysis.R first.")
  }
  
  plasma_data <- read.csv("plots/04_Biomarker_Analysis/plasma_expression_data.csv", 
                         stringsAsFactors = FALSE)
  
  # Validate plasma data
  if(nrow(plasma_data) == 0) {
    stop("Plasma data is empty!")
  }
  
  required_cols <- c("gene_symbol", "abundance", "source", "abundance_type")
  missing_cols <- setdiff(required_cols, colnames(plasma_data))
  if(length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Filter out HPA Immunoassay data
  original_sources <- unique(plasma_data$source)
  plasma_data <- plasma_data[plasma_data$source != "HPA_Immunoassay", ]
  filtered_sources <- unique(plasma_data$source)
  
  message(paste("Excluded HPA_Immunoassay from analysis"))
  message(paste("Original sources:", paste(original_sources, collapse = ", ")))
  message(paste("Filtered sources:", paste(filtered_sources, collapse = ", ")))
  
  # Load biomarker genes
  biomarker_genes <- load_biomarker_genes()
  
  message(paste("Loaded", nrow(plasma_data), "plasma measurements"))
  message(paste("Creating ranked plots for", length(unique(plasma_data$source)), "sources"))
  
  # Create individual plots for each source
  sources <- unique(plasma_data$source)
  plot_list <- list()
  
  for(source in sources) {
    message(paste("Creating plot for", source, "..."))
    
    source_data <- plasma_data[plasma_data$source == source, ]
    
    # Remove duplicates by taking mean abundance for duplicate genes
    source_data <- source_data %>%
      group_by(gene_symbol) %>%
      summarise(
        abundance = mean(abundance, na.rm = TRUE),
        abundance_type = first(abundance_type),
        .groups = 'drop'
      )
    
    # Create individual plot
    plot_obj <- create_ranked_plot(source_data, source, biomarker_genes, 
                                 source_data$abundance_type[1])
    
    if(!is.null(plot_obj)) {
      plot_list[[source]] <- plot_obj
      
      # Save individual plot (using PNG to save disk space)
      safe_name <- gsub("[^A-Za-z0-9_]", "_", source)
      png(file.path(output_dir, paste0("ranked_abundance_", safe_name, ".png")), 
          width = 12, height = 8, units = "in", res = 300)
      print(plot_list[[source]])
      dev.off()
    } else {
      warning(paste("Skipping", source, "due to data issues"))
    }
  }
  
  # Create a combined plot with adaptive layout
  if(length(plot_list) > 0) {
    message("Creating combined plot...")
    
    # Determine optimal layout
    n_plots <- length(plot_list)
    if(n_plots <= 2) {
      ncol <- 1; nrow <- n_plots
    } else if(n_plots <= 4) {
      ncol <- 2; nrow <- ceiling(n_plots/2)
    } else {
      ncol <- 2; nrow <- ceiling(n_plots/2)
    }
    
    # Adjust figure size based on layout
    fig_height <- max(10, nrow * 5)
    fig_width <- max(12, ncol * 8)
    
    png(file.path(output_dir, "ranked_abundance_all_sources.png"), 
        width = fig_width, height = fig_height, units = "in", res = 300)
    
    grid.arrange(grobs = plot_list, ncol = ncol)
    dev.off()
  }
  
  # Create improved summary statistics
  summary_stats <- create_improved_summary(plasma_data, biomarker_genes)
  
  write.csv(summary_stats, file.path(output_dir, "ranked_abundance_summary_improved.csv"), 
           row.names = FALSE)
  
  # Create an enhanced comparison plot
  message("Creating enhanced comparison plots...")
  
  # Dynamic range comparison
  png(file.path(output_dir, "dynamic_range_comparison.png"), 
      width = 12, height = 8, units = "in", res = 300)
  
  p_comparison <- ggplot(summary_stats, aes(x = reorder(source, dynamic_range), 
                                           y = dynamic_range, fill = source)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_viridis_d(name = "Source") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = "Protein Concentration Dynamic Range by Source",
      subtitle = "Log10 orders of magnitude between highest and lowest concentrations",
      x = "Data Source",
      y = "Dynamic Range (log10 orders of magnitude)"
    ) +
    geom_text(aes(label = paste0(dynamic_range, " orders")), 
              hjust = -0.1, size = 3)
  
  print(p_comparison)
  dev.off()
  
  # Biomarker coverage comparison
  png(file.path(output_dir, "biomarker_coverage_comparison.png"), 
      width = 12, height = 8, units = "in", res = 300)
  
  p_coverage <- ggplot(summary_stats, aes(x = reorder(source, biomarker_coverage_pct), 
                                         y = biomarker_coverage_pct, fill = source)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_viridis_d(name = "Source") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = "Biomarker Coverage by Source",
      subtitle = paste("Percentage of", length(biomarker_genes), "biomarker genes detected"),
      x = "Data Source",
      y = "Biomarker Coverage (%)"
    ) +
    geom_text(aes(label = paste0(biomarker_coverage_pct, "%")), 
              hjust = -0.1, size = 3) +
    ylim(0, 100)
  
  print(p_coverage)
  dev.off()
  
  # Print improved summary
  message("\n" , rep("=", 60))
  message("IMPROVED Ranked Abundance Analysis Summary")
  message(rep("=", 60))
  message("\nData Quality Summary:")
  message(paste("✓ Total plasma measurements:", nrow(plasma_data)))
  message(paste("✓ Unique genes across all sources:", length(unique(plasma_data$gene_symbol))))
  message(paste("✓ Data sources analyzed:", length(unique(plasma_data$source))))
  message(paste("✓ Biomarker genes in reference set:", length(biomarker_genes)))
  
  message("\nDetailed Source Statistics:")
  for(i in 1:nrow(summary_stats)) {
    row <- summary_stats[i, ]
    message(sprintf("\n%s:", row$source))
    message(sprintf("  • Proteins detected: %s", scales::comma(row$n_proteins)))
    message(sprintf("  • Biomarkers found: %d/%d (%.1f%%)", 
                   row$unique_biomarkers, row$total_biomarkers, row$biomarker_coverage_pct))
    message(sprintf("  • Abundance range: %.3g - %.3g %s", 
                   row$min_abundance, row$max_abundance, row$abundance_type))
    message(sprintf("  • Dynamic range: %.1f orders of magnitude", row$dynamic_range))
  }
  
  return(list(plots = plot_list, summary = summary_stats))
}

# Function to generate detailed report
generate_detailed_report <- function(results) {
  summary_stats <- results$summary
  
  report_file <- file.path(output_dir, "analysis_report.txt")
  
  sink(report_file)
  cat("IMPROVED RANKED ABUNDANCE ANALYSIS - DETAILED REPORT\n")
  cat(rep("=", 60), "\n")
  cat("Generated:", Sys.time(), "\n\n")
  
  cat("ANALYSIS OVERVIEW:\n")
  cat("-----------------\n")
  cat("This analysis creates ranked abundance plots (waterfall plots) showing\n")
  cat("protein concentrations from highest to lowest, with biomarker genes\n")
  cat("highlighted in red to show their position in the overall distribution.\n\n")
  
  cat("IMPROVEMENTS IN THIS VERSION:\n")
  cat("----------------------------\n")
  cat("• Fixed biomarker counting bug (now counts unique genes per source)\n")
  cat("• Updated deprecated ggplot2 aesthetics (size → linewidth)\n")
  cat("• Enhanced data validation and error handling\n")
  cat("• Improved log transformation with appropriate offsets\n")
  cat("• Better dynamic range calculation (handles edge cases)\n")
  cat("• Adaptive plot layout for different numbers of sources\n")
  cat("• More detailed reporting and statistics\n\n")
  
  cat("SOURCE COMPARISON:\n")
  cat("------------------\n")
  print(summary_stats)
  
  sink()
  
  message(paste("\nDetailed report saved to:", report_file))
}

# Main execution
message("Starting IMPROVED ranked abundance plot generation...")
results <- create_all_ranked_plots()
generate_detailed_report(results)
message("Improved ranked abundance plots completed!")
message(paste("Results saved to:", output_dir)) 