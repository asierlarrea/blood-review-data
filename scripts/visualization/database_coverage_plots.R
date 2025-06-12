#!/usr/bin/env Rscript

# Database Coverage Analysis - Descriptive Plots for Each Database
# Focus on plasma and serum data as described in the manuscript

# Load required packages
source("scripts/utilities/load_packages.R")
required_packages <- c("ggplot2", "dplyr", "tidyr", "viridis", "RColorBrewer", 
                      "scales", "gridExtra", "stringr", "cowplot", "ggrepel")
load_packages(required_packages)

# Create output directory
output_dir <- "outputs/plots/01_Database_Analysis/Individual_Database_Coverage"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define consistent colors and themes
database_colors <- c(
  "PeptideAtlas" = "#E31A1C",
  "HPA (MS)" = "#1F78B4", 
  "HPA (Immunoassay)" = "#33A02C",
  "HPA (PEA)" = "#FF7F00",
  "PaxDB" = "#6A3D9A",
  "GPMDB" = "#B15928"
)

# Custom theme for plots
custom_theme <- theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.3)
  )

# Function to read and process database files
read_database_files <- function() {
  message("Reading database files...")
  
  databases <- list()
  
  # 1. PeptideAtlas
  if(file.exists(get_data_path("PeptideAtlas.csv"))) {
    peptideatlas <- read.csv(get_data_path("PeptideAtlas.csv"), stringsAsFactors = FALSE)
    peptideatlas_clean <- peptideatlas %>%
      filter(!is.na(norm_PSMs_per_100K) & norm_PSMs_per_100K > 0) %>%
      mutate(
        Database = "PeptideAtlas",
        Protein_Count = n(),
        Abundance = norm_PSMs_per_100K,
        Source = "MS"
      )
    databases$PeptideAtlas <- peptideatlas_clean
    message(paste("PeptideAtlas:", nrow(peptideatlas_clean), "proteins"))
  }
  
  # 2. HPA MS data
  if(file.exists(get_data_path("HPA_MS.csv"))) {
    hpa_ms <- read.csv(get_data_path("HPA_MS.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    # Extract concentration values and convert to numeric (handle special characters)
    hpa_ms$Concentration_Clean <- gsub("[^0-9.]", "", hpa_ms$Concentration)
    hpa_ms$Concentration_Value <- suppressWarnings(as.numeric(hpa_ms$Concentration_Clean))
    hpa_ms_clean <- hpa_ms %>%
      filter(!is.na(Concentration_Value) & Concentration_Value > 0) %>%
      mutate(
        Database = "HPA (MS)",
        Abundance = Concentration_Value,
        Source = "MS"
      )
    databases$HPA_MS <- hpa_ms_clean
    message(paste("HPA MS:", nrow(hpa_ms_clean), "proteins"))
  }
  
  # 3. HPA Immunoassay data
  if(file.exists(get_data_path("HPA_Immunoassay.csv"))) {
    hpa_immuno <- read.csv(get_data_path("HPA_Immunoassay.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    # Create abundance values for immunoassay data
    hpa_immuno_clean <- hpa_immuno %>%
      filter(!is.na(Gene) & Gene != "") %>%
      mutate(
        Database = "HPA (Immunoassay)",
        Abundance = row_number(), # Rank-based abundance
        Source = "Immunoassay"
      )
    databases$HPA_Immunoassay <- hpa_immuno_clean
    message(paste("HPA Immunoassay:", nrow(hpa_immuno_clean), "proteins"))
  }
  
  # 4. HPA PEA data
  if(file.exists(get_data_path("HPA_PEA.csv"))) {
    hpa_pea <- read.csv(get_data_path("HPA_PEA.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    hpa_pea_clean <- hpa_pea %>%
      filter(!is.na(Gene) & Gene != "") %>%
      mutate(
        Database = "HPA (PEA)",
        Abundance = row_number(), # Rank-based abundance
        Source = "PEA"
      )
    databases$HPA_PEA <- hpa_pea_clean
    message(paste("HPA PEA:", nrow(hpa_pea_clean), "proteins"))
  }
  
  # 5. PaxDB Plasma
  if(file.exists(get_data_path("PaxDb_Plasma.csv"))) {
    paxdb <- read.csv(get_data_path("PaxDb_Plasma.csv"), stringsAsFactors = FALSE)
    paxdb_clean <- paxdb %>%
      filter(!is.na(abundance) & abundance > 0) %>%
      mutate(
        Database = "PaxDB (Plasma)",
        Abundance = abundance,
        Source = "MS"
      )
    databases$PaxDB <- paxdb_clean
    message(paste("PaxDB:", nrow(paxdb_clean), "proteins"))
  }
  
  # 6. GPMDB Plasma
  if(file.exists(get_data_path("GPMDB_Plasma.csv"))) {
    tryCatch({
      # Read GPMDB file properly (skip header lines)
      gpmdb_lines <- readLines(get_data_path("GPMDB_Plasma.csv"), encoding = "UTF-8")
      header_line <- which(grepl("^#,accession,total", gpmdb_lines))
      
      if(length(header_line) > 0) {
        gpmdb <- read.csv(get_data_path("GPMDB_Plasma.csv"), skip = header_line, header = FALSE, 
                         stringsAsFactors = FALSE, fileEncoding = "UTF-8")
        colnames(gpmdb) <- c("rank", "accession", "total", "log_e", "EC", "description")
        
        gpmdb_clean <- gpmdb %>%
          filter(!is.na(total) & is.numeric(total) & total > 0) %>%
          mutate(
            Database = "GPMDB (Plasma)",
            Abundance = as.numeric(total),
            Source = "MS"
          )
        databases$GPMDB <- gpmdb_clean
        message(paste("GPMDB:", nrow(gpmdb_clean), "proteins"))
      }
    }, error = function(e) {
      message(paste("Warning: Could not process GPMDB file:", e$message))
    })
  }
  
  return(databases)
}

# Function to create database summary statistics
create_database_summary <- function(databases) {
  summary_stats <- data.frame()
  
  for(db_name in names(databases)) {
    db_data <- databases[[db_name]]
    
    if(nrow(db_data) > 0) {
      stats <- data.frame(
        Database = db_name,
        Total_Proteins = nrow(db_data),
        Source = unique(db_data$Source)[1],
        Min_Abundance = min(db_data$Abundance, na.rm = TRUE),
        Max_Abundance = max(db_data$Abundance, na.rm = TRUE),
        Median_Abundance = median(db_data$Abundance, na.rm = TRUE),
        Mean_Abundance = mean(db_data$Abundance, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      summary_stats <- rbind(summary_stats, stats)
    }
  }
  
  return(summary_stats)
}

# Plot 1: Database Coverage Comparison
plot_database_coverage <- function(summary_stats) {
  # Reorder databases by protein count
  summary_stats$Database <- factor(summary_stats$Database, 
                                  levels = summary_stats$Database[order(summary_stats$Total_Proteins)])
  
  # Create bar plot
  p1 <- ggplot(summary_stats, aes(x = Database, y = Total_Proteins, fill = Database)) +
    geom_col(alpha = 0.8, width = 0.7) +
    scale_fill_manual(values = database_colors) +
    scale_y_continuous(labels = comma_format()) +
    labs(
      title = "Protein Coverage Across Databases",
      subtitle = "Total number of proteins quantified in plasma/serum",
      x = "Database",
      y = "Number of Proteins",
      caption = "Data sources: PeptideAtlas (MS), HPA (MS/Immunoassay/PEA), PaxDB (MS), GPMDB (MS)"
    ) +
    custom_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    geom_text(aes(label = comma(Total_Proteins)), 
              vjust = -0.5, size = 3.5, fontface = "bold")
  
  return(p1)
}

# Plot 2: Abundance Distribution by Database
plot_abundance_distributions <- function(databases) {
  # Combine data for plotting
  combined_data <- data.frame()
  
  for(db_name in names(databases)) {
    db_data <- databases[[db_name]]
    if(nrow(db_data) > 0 && "Abundance" %in% colnames(db_data)) {
      temp_data <- data.frame(
        Database = db_name,
        Abundance = log10(db_data$Abundance + 1),
        Source = db_data$Source[1],
        stringsAsFactors = FALSE
      )
      combined_data <- rbind(combined_data, temp_data)
    }
  }
  
  # Create violin plots
  p2 <- ggplot(combined_data, aes(x = Database, y = Abundance, fill = Database)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.size = 0.5) +
    scale_fill_manual(values = database_colors) +
    labs(
      title = "Protein Abundance Distributions",
      subtitle = "Log10-transformed abundance values across databases",
      x = "Database",
      y = "Log10(Abundance + 1)",
      caption = "Note: Different databases use different abundance units and scales"
    ) +
    custom_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  return(p2)
}

# Plot 3: Technology Platform Comparison
plot_technology_platforms <- function(summary_stats) {
  # Group by technology
  tech_summary <- summary_stats %>%
    group_by(Source) %>%
    summarise(
      Databases = n(),
      Total_Proteins = sum(Total_Proteins),
      Avg_Proteins = round(mean(Total_Proteins)),
      .groups = 'drop'
    )
  
  p3 <- ggplot(tech_summary, aes(x = Source, y = Total_Proteins, fill = Source)) +
    geom_col(alpha = 0.8, width = 0.6) +
    scale_fill_manual(values = c("MS" = "#2E86AB", "Immunoassay" = "#A23B72", "PEA" = "#F18F01")) +
    scale_y_continuous(labels = comma_format()) +
    labs(
      title = "Protein Coverage by Technology Platform",
      subtitle = "Cumulative proteins across databases by detection method",
      x = "Technology Platform",
      y = "Total Proteins (Cumulative)",
      caption = "MS: Mass Spectrometry, PEA: Proximity Extension Assay"
    ) +
    custom_theme +
    theme(legend.position = "none") +
    geom_text(aes(label = paste0(comma(Total_Proteins), "\n(", Databases, " DBs)")), 
              vjust = -0.5, size = 3.5, fontface = "bold")
  
  return(p3)
}

# Plot 4: Individual Database Details
plot_individual_database_details <- function(databases, summary_stats) {
  plots_list <- list()
  
  for(db_name in names(databases)) {
    db_data <- databases[[db_name]]
    db_stats <- summary_stats[summary_stats$Database == db_name, ]
    
    if(nrow(db_data) > 0 && nrow(db_stats) > 0) {
      # Create histogram of abundance distribution
      if("Abundance" %in% colnames(db_data)) {
        p <- ggplot(db_data, aes(x = log10(Abundance + 1))) +
          geom_histogram(bins = 30, fill = database_colors[db_name], alpha = 0.7, color = "white") +
          geom_vline(aes(xintercept = log10(db_stats$Median_Abundance + 1)), 
                    color = "red", linetype = "dashed", size = 1) +
          labs(
            title = paste(db_name, "- Protein Abundance Distribution"),
            subtitle = paste("Total proteins:", comma(db_stats$Total_Proteins), 
                           "| Median abundance:", scientific(db_stats$Median_Abundance, digits = 2)),
            x = "Log10(Abundance + 1)",
            y = "Number of Proteins",
            caption = paste("Source:", db_stats$Source, "| Range:", 
                          scientific(db_stats$Min_Abundance, digits = 2), "to", 
                          scientific(db_stats$Max_Abundance, digits = 2))
          ) +
          custom_theme
      } else {
        # Simple count plot for databases without abundance data
        p <- ggplot(data.frame(x = 1, y = db_stats$Total_Proteins), aes(x = x, y = y)) +
          geom_col(fill = database_colors[db_name], alpha = 0.7, width = 0.5) +
          labs(
            title = paste(db_name, "- Protein Count"),
            subtitle = paste("Total proteins:", comma(db_stats$Total_Proteins)),
            x = "",
            y = "Number of Proteins",
            caption = paste("Source:", db_stats$Source)
          ) +
          custom_theme +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      
      plots_list[[db_name]] <- p
    }
  }
  
  return(plots_list)
}

# Plot 5: Database Characteristics Radar Plot
plot_database_characteristics <- function(summary_stats) {
  # Normalize metrics for radar plot
  radar_data <- summary_stats %>%
    mutate(
      Norm_Proteins = (Total_Proteins - min(Total_Proteins)) / (max(Total_Proteins) - min(Total_Proteins)),
      Norm_Range = log10(Max_Abundance/Min_Abundance),
      Norm_Range = (Norm_Range - min(Norm_Range, na.rm = TRUE)) / (max(Norm_Range, na.rm = TRUE) - min(Norm_Range, na.rm = TRUE))
    ) %>%
    select(Database, Norm_Proteins, Norm_Range) %>%
    pivot_longer(cols = c(Norm_Proteins, Norm_Range), names_to = "Metric", values_to = "Value")
  
  # Create grouped bar plot instead of radar (simpler and clearer)
  p5 <- ggplot(radar_data, aes(x = Database, y = Value, fill = Metric)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Norm_Proteins" = "#1B9E77", "Norm_Range" = "#D95F02"),
                     labels = c("Protein Coverage", "Dynamic Range")) +
    labs(
      title = "Database Characteristics Comparison",
      subtitle = "Normalized protein coverage and dynamic range",
      x = "Database",
      y = "Normalized Score (0-1)",
      fill = "Characteristic"
    ) +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p5)
}

# Main execution function
main <- function() {
  message("Starting database coverage analysis...")
  
  # Read database files
  databases <- read_database_files()
  
  if(length(databases) == 0) {
    stop("No database files found or processed successfully.")
  }
  
  # Create summary statistics
  summary_stats <- create_database_summary(databases)
  
  # Generate plots
  message("Generating plots...")
  
  # Plot 1: Coverage comparison
  p1 <- plot_database_coverage(summary_stats)
  ggsave(file.path(output_dir, "01_database_coverage_comparison.png"), 
         p1, width = 10, height = 6, dpi = 300, bg = "white")
  
  # Plot 2: Abundance distributions
  if(any(sapply(databases, function(x) "Abundance" %in% colnames(x)))) {
    p2 <- plot_abundance_distributions(databases)
    ggsave(file.path(output_dir, "02_abundance_distributions.png"), 
           p2, width = 12, height = 6, dpi = 300, bg = "white")
  }
  
  # Plot 3: Technology platforms
  p3 <- plot_technology_platforms(summary_stats)
  ggsave(file.path(output_dir, "03_technology_platforms.png"), 
         p3, width = 8, height = 6, dpi = 300, bg = "white")
  
  # Plot 4: Individual database details
  individual_plots <- plot_individual_database_details(databases, summary_stats)
  for(db_name in names(individual_plots)) {
    ggsave(file.path(output_dir, paste0("04_", db_name, "_details.png")), 
           individual_plots[[db_name]], width = 8, height = 6, dpi = 300, bg = "white")
  }
  
  # Plot 5: Database characteristics
  p5 <- plot_database_characteristics(summary_stats)
  ggsave(file.path(output_dir, "05_database_characteristics.png"), 
         p5, width = 10, height = 6, dpi = 300, bg = "white")
  
  # Create combined overview plot
  if(length(individual_plots) >= 4) {
    combined_plot <- plot_grid(
      p1, p3,
      individual_plots[[1]], individual_plots[[2]],
      ncol = 2, nrow = 2,
      labels = c("A", "B", "C", "D"),
      label_size = 14
    )
    
    ggsave(file.path(output_dir, "00_database_overview_combined.png"), 
           combined_plot, width = 16, height = 12, dpi = 300, bg = "white")
  }
  
  # Save summary statistics
  write.csv(summary_stats, file.path(output_dir, "database_summary_statistics.csv"), row.names = FALSE)
  
  # Create detailed report
  cat("DATABASE COVERAGE ANALYSIS REPORT\n", 
      "=================================\n\n",
      file = file.path(output_dir, "database_analysis_report.txt"))
  
  for(i in 1:nrow(summary_stats)) {
    stats <- summary_stats[i, ]
    cat(sprintf("%s:\n", stats$Database),
        sprintf("  - Total proteins: %s\n", comma(stats$Total_Proteins)),
        sprintf("  - Technology: %s\n", stats$Source),
        sprintf("  - Abundance range: %s to %s\n", 
                scientific(stats$Min_Abundance, digits = 2),
                scientific(stats$Max_Abundance, digits = 2)),
        sprintf("  - Median abundance: %s\n\n", scientific(stats$Median_Abundance, digits = 2)),
        file = file.path(output_dir, "database_analysis_report.txt"), append = TRUE)
  }
  
  message(paste("Analysis complete! Results saved to:", output_dir))
  message("Generated files:")
  message("  - 01_database_coverage_comparison.png")
  message("  - 02_abundance_distributions.png") 
  message("  - 03_technology_platforms.png")
  message("  - 04_[Database]_details.png (for each database)")
  message("  - 05_database_characteristics.png")
  message("  - 00_database_overview_combined.png")
  message("  - database_summary_statistics.csv")
  message("  - database_analysis_report.txt")
  
  return(list(databases = databases, summary_stats = summary_stats))
}

# Run the analysis
if(!interactive()) {
  main()
} 
