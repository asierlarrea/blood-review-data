#!/usr/bin/env Rscript

# Ranked Abundance Plots for Plasma Protein Concentrations
# Creates waterfall-style plots showing protein concentrations ranked from highest to lowest
# with biomarker genes highlighted in red

# Set CRAN mirror
if(is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

# Function to install packages if not available
install_if_missing <- function(packages) {
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      install.packages(pkg, dependencies = TRUE)
      if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
        stop(paste("Failed to install package:", pkg))
      }
    }
  }
}

# Install required packages
required_packages <- c("ggplot2", "dplyr", "scales", "gridExtra", "viridis")
install_if_missing(required_packages)

# Create output directory
output_dir <- "plots/05_Ranked_Abundance"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to load biomarker genes
load_biomarker_genes <- function() {
  if(!file.exists("biomarkers_list.tsv")) {
    stop("biomarkers_list.tsv not found!")
  }
  
  biomarkers <- read.delim("biomarkers_list.tsv", stringsAsFactors = FALSE)
  biomarker_genes <- unique(toupper(biomarkers$To))
  biomarker_genes <- biomarker_genes[!is.na(biomarker_genes) & biomarker_genes != ""]
  
  return(biomarker_genes)
}

# Function to create ranked abundance plot for a single source
create_ranked_plot <- function(data, source_name, biomarker_genes, abundance_type) {
  # Remove duplicates and rank by abundance
  source_data <- data %>%
    arrange(desc(abundance)) %>%
    mutate(
      rank = row_number(),
      is_biomarker = gene_symbol %in% biomarker_genes,
      log_abundance = log10(abundance + ifelse(abundance_type == "mg/L", 0.001, 1))
    )
  
  # Get biomarker subset for highlighting
  biomarker_subset <- source_data[source_data$is_biomarker, ]
  
  # Create the plot
  p <- ggplot(source_data, aes(x = rank, y = log_abundance)) +
    # Background proteins as grey area
    geom_area(fill = "lightgrey", alpha = 0.7) +
    # Add grid
    geom_hline(yintercept = seq(floor(min(source_data$log_abundance)), 
                               ceiling(max(source_data$log_abundance)), 1), 
               color = "white", alpha = 0.8, size = 0.3) +
    # Highlight biomarkers as red vertical lines
    geom_segment(data = biomarker_subset, 
                aes(x = rank, xend = rank, y = 0, yend = log_abundance),
                color = "red", alpha = 0.8, size = 0.8) +
    # Formatting
    scale_y_continuous(
      name = paste0("Concentration (", abundance_type, ", log scale)"),
      labels = function(x) {
        if(abundance_type == "mg/L") {
          parse(text = paste0("10^", x))
        } else if(abundance_type == "ppm") {
          parse(text = paste0("10^", x))
        } else if(abundance_type == "Normalized PSMs per 100K") {
          parse(text = paste0("10^", x))
        } else {
          parse(text = paste0("10^", x))
        }
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
                       sum(source_data$is_biomarker), " biomarkers highlighted in red")
    )
  
  return(p)
}

# Function to create combined plots
create_all_ranked_plots <- function() {
  message("Loading plasma expression data...")
  
  # Load plasma data from previous analysis
  if(!file.exists("plots/04_Biomarker_Analysis/plasma_expression_data.csv")) {
    stop("Plasma expression data not found. Please run biomarker_plasma_analysis.R first.")
  }
  
  plasma_data <- read.csv("plots/04_Biomarker_Analysis/plasma_expression_data.csv", 
                         stringsAsFactors = FALSE)
  
  # Load biomarker genes
  biomarker_genes <- load_biomarker_genes()
  
  message(paste("Loaded", nrow(plasma_data), "plasma measurements"))
  message(paste("Creating ranked plots for", length(unique(plasma_data$source)), "sources"))
  
  # Create individual plots for each source
  sources <- unique(plasma_data$source)
  plot_list <- list()
  
  for(source in sources) {
    message(paste("Creating plot for", source))
    
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
    plot_list[[source]] <- create_ranked_plot(source_data, source, biomarker_genes, 
                                             source_data$abundance_type[1])
    
    # Save individual plot
    safe_name <- gsub("[^A-Za-z0-9_]", "_", source)
    tiff(file.path(output_dir, paste0("ranked_abundance_", safe_name, ".tiff")), 
         width = 12, height = 8, units = "in", res = 600, compression = "lzw")
    print(plot_list[[source]])
    dev.off()
  }
  
  # Create a combined plot with all sources
  message("Creating combined plot...")
  
  tiff(file.path(output_dir, "ranked_abundance_all_sources.tiff"), 
       width = 16, height = 20, units = "in", res = 600, compression = "lzw")
  
  # Arrange plots in a grid
  if(length(plot_list) <= 4) {
    grid.arrange(grobs = plot_list, ncol = 2)
  } else {
    grid.arrange(grobs = plot_list, ncol = 2)
  }
  
  dev.off()
  
  # Create summary statistics
  summary_stats <- plasma_data %>%
    group_by(source) %>%
    summarise(
      n_proteins = n_distinct(gene_symbol),
      n_biomarkers = sum(gene_symbol %in% biomarker_genes),
      biomarker_pct = round(100 * n_biomarkers / length(biomarker_genes), 1),
      min_abundance = min(abundance, na.rm = TRUE),
      max_abundance = max(abundance, na.rm = TRUE),
      median_abundance = median(abundance, na.rm = TRUE),
      dynamic_range = round(log10(max_abundance / min_abundance), 2),
      abundance_type = first(abundance_type),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_proteins))
  
  write.csv(summary_stats, file.path(output_dir, "ranked_abundance_summary.csv"), 
           row.names = FALSE)
  
  # Create a comparison plot showing dynamic ranges
  message("Creating dynamic range comparison...")
  
  tiff(file.path(output_dir, "dynamic_range_comparison.tiff"), 
       width = 12, height = 8, units = "in", res = 600, compression = "lzw")
  
  p_comparison <- ggplot(summary_stats, aes(x = reorder(source, dynamic_range), 
                                           y = dynamic_range, fill = source)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_viridis_d(name = "Source") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = "Protein Concentration Dynamic Range by Source",
      subtitle = "Orders of magnitude between highest and lowest concentrations",
      x = "Data Source",
      y = "Dynamic Range (log10 orders of magnitude)"
    ) +
    geom_text(aes(label = paste0(dynamic_range, " orders")), 
              hjust = -0.1, size = 3)
  
  print(p_comparison)
  dev.off()
  
  # Print summary
  message("\nRanked Abundance Analysis Summary:")
  message("===================================")
  print(summary_stats)
  
  return(list(plots = plot_list, summary = summary_stats))
}

# Main execution
message("Starting ranked abundance plot generation...")
results <- create_all_ranked_plots()
message("Ranked abundance plots completed!")
message(paste("Results saved to:", output_dir)) 