#!/usr/bin/env Rscript

# HPA PEA Variation vs MS/MS Variability Correlation Analysis
# Correlating PEA variation metrics with MS/MS coefficient of variation

# Load required packages
source("scripts/utilities/load_packages.R")
required_packages <- c("ggplot2", "dplyr", "tidyr", "scales", "corrplot", 
                      "ggpubr", "viridis", "RColorBrewer", "gridExtra", "cowplot")
load_packages(required_packages)

# Create output directory
output_dir <- "outputs/plots/01_Database_Analysis/PEA_MSMS_Variability"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to calculate coefficient of variation
calculate_cv <- function(values) {
  values_numeric <- as.numeric(values)
  values_clean <- values_numeric[!is.na(values_numeric) & values_numeric > 0]
  if (length(values_clean) < 2) return(NA)
  cv <- sd(values_clean) / mean(values_clean)
  return(cv)
}

# Function to safely extract numeric columns
extract_intensity_columns <- function(df, exclude_cols) {
  numeric_cols <- sapply(df, function(x) {
    if (is.character(x) || is.factor(x)) {
      # Try to convert to numeric
      numeric_test <- suppressWarnings(as.numeric(as.character(x)))
      return(sum(!is.na(numeric_test)) > 0)
    }
    return(is.numeric(x))
  })
  
  # Get column names that are numeric and not in exclude list
  intensity_cols <- names(df)[numeric_cols & !names(df) %in% exclude_cols]
  return(intensity_cols)
}

cat("Loading and processing datasets...\n")

# Load HPA PEA data
hpa_pea <- read.csv(get_data_path("HPA_PEA.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
cat("HPA PEA loaded:", nrow(hpa_pea), "proteins\n")

# Clean up gene names in HPA PEA
hpa_pea$Gene_Clean <- toupper(trimws(hpa_pea$Gene))

# Initialize list to store MS/MS variability data
msms_variability <- list()

# Process datasets with multiple intensity columns
datasets_to_process <- list(
  "PXD040957_CD8" = list(file = "PXD040957_CD8.csv", exclude = c("Protein.IDs", "Gene.Names", "Gene.names", "Protein.names")),
  "PXD040957_Macrophages" = list(file = "PXD040957_Macrophages.csv", exclude = c("Protein.IDs", "Gene.Names", "Gene.names", "Protein.names")),
  "PXD025174" = list(file = "PXD025174.csv", exclude = c("Protein.IDs", "Gene.Names", "Gene.names", "Protein.names"))
)

for (dataset_name in names(datasets_to_process)) {
  file_path <- datasets_to_process[[dataset_name]]$file
  exclude_cols <- datasets_to_process[[dataset_name]]$exclude
  
  if (file.exists(file_path)) {
    cat(paste("Processing", dataset_name, "...\n"))
    
    tryCatch({
      # Read dataset
      data <- read.csv(file_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
      
      # Find intensity columns
      intensity_cols <- extract_intensity_columns(data, exclude_cols)
      
      if (length(intensity_cols) >= 2) {
        cat(paste("  Found", length(intensity_cols), "intensity columns\n"))
        
        # Extract gene names - try different possible column names
        gene_col <- NULL
        possible_gene_cols <- c("Gene.Names", "Gene.names", "Gene", "GENE", "gene")
        for (col in possible_gene_cols) {
          if (col %in% names(data)) {
            gene_col <- col
            break
          }
        }
        
        if (!is.null(gene_col)) {
          # Calculate CV for each protein
          cv_data <- data %>%
            rowwise() %>%
            mutate(
              Gene_Clean = toupper(trimws(strsplit(as.character(.data[[gene_col]]), ";")[[1]][1])),
              CV_MSMS = calculate_cv(c_across(all_of(intensity_cols))),
              Dataset = dataset_name,
              N_Measurements = sum(!is.na(c_across(all_of(intensity_cols))))
            ) %>%
            filter(!is.na(CV_MSMS) & CV_MSMS > 0 & Gene_Clean != "" & Gene_Clean != "NA") %>%
            select(Gene_Clean, CV_MSMS, Dataset, N_Measurements)
          
          msms_variability[[dataset_name]] <- as.data.frame(cv_data)
          cat(paste("  Calculated CV for", nrow(cv_data), "proteins\n"))
        } else {
          cat(paste("  No gene column found in", dataset_name, "\n"))
        }
      } else {
        cat(paste("  Insufficient intensity columns in", dataset_name, "\n"))
      }
    }, error = function(e) {
      cat(paste("  Error processing", dataset_name, ":", e$message, "\n"))
    })
  }
}

# Combine all MS/MS variability data
if (length(msms_variability) > 0) {
  msms_combined <- do.call(rbind, msms_variability)
  cat("Combined MS/MS data:", nrow(msms_combined), "protein measurements\n")
  
  # Average CV across datasets for proteins measured in multiple datasets
  msms_avg <- msms_combined %>%
    group_by(Gene_Clean) %>%
    summarise(
      CV_MSMS_mean = mean(CV_MSMS, na.rm = TRUE),
      CV_MSMS_median = median(CV_MSMS, na.rm = TRUE),
      N_Datasets = n(),
      Total_Measurements = sum(N_Measurements),
      .groups = 'drop'
    )
  
  cat("Unique proteins with MS/MS CV:", nrow(msms_avg), "\n")
  
  # Merge with HPA PEA data
  merged_data <- hpa_pea %>%
    inner_join(msms_avg, by = "Gene_Clean") %>%
    filter(!is.na(Variation.between.individuals) & 
           !is.na(Variation.within.individuals) & 
           !is.na(CV_MSMS_mean))
  
  cat("Proteins with both PEA and MS/MS variability data:", nrow(merged_data), "\n")
  
  if (nrow(merged_data) > 0) {
    # Custom theme for plots
    custom_theme <- theme_minimal() +
      theme(
        text = element_text(size = 12, family = "Arial"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
      )
    
    # 1. Correlation plot: PEA between-individual variation vs MS/MS CV
    p1 <- ggplot(merged_data, aes(x = Variation.between.individuals, y = CV_MSMS_mean)) +
      geom_point(alpha = 0.6, size = 2, color = "#2E86AB") +
      geom_smooth(method = "lm", color = "#A23B72", fill = "#F18F01", alpha = 0.3) +
      stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.9, size = 4) +
      labs(
        title = "PEA Between-Individual Variation vs MS/MS Coefficient of Variation",
        subtitle = paste("Correlation analysis for", nrow(merged_data), "proteins"),
        x = "PEA Variation Between Individuals",
        y = "MS/MS Coefficient of Variation (mean)"
      ) +
      custom_theme
    
    # 2. Correlation plot: PEA within-individual variation vs MS/MS CV
    p2 <- ggplot(merged_data, aes(x = Variation.within.individuals, y = CV_MSMS_mean)) +
      geom_point(alpha = 0.6, size = 2, color = "#2E86AB") +
      geom_smooth(method = "lm", color = "#A23B72", fill = "#F18F01", alpha = 0.3) +
      stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.9, size = 4) +
      labs(
        title = "PEA Within-Individual Variation vs MS/MS Coefficient of Variation",
        subtitle = paste("Correlation analysis for", nrow(merged_data), "proteins"),
        x = "PEA Variation Within Individuals",
        y = "MS/MS Coefficient of Variation (mean)"
      ) +
      custom_theme
    
    # 3. PEA variation comparison (between vs within individuals)
    p3 <- ggplot(merged_data, aes(x = Variation.between.individuals, y = Variation.within.individuals)) +
      geom_point(alpha = 0.6, size = 2, color = "#2E86AB") +
      geom_smooth(method = "lm", color = "#A23B72", fill = "#F18F01", alpha = 0.3) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.9, size = 4) +
      labs(
        title = "PEA Variation: Between vs Within Individuals",
        subtitle = "Dashed line indicates equal variation",
        x = "Variation Between Individuals",
        y = "Variation Within Individuals"
      ) +
      custom_theme
    
    # 4. Distribution comparison plot
    variation_long <- merged_data %>%
      select(Gene_Clean, Variation.between.individuals, Variation.within.individuals, CV_MSMS_mean) %>%
      pivot_longer(cols = c(Variation.between.individuals, Variation.within.individuals, CV_MSMS_mean),
                   names_to = "Variation_Type", values_to = "Variation_Value") %>%
      mutate(
        Variation_Type = case_when(
          Variation_Type == "Variation.between.individuals" ~ "PEA Between Individuals",
          Variation_Type == "Variation.within.individuals" ~ "PEA Within Individuals",
          Variation_Type == "CV_MSMS_mean" ~ "MS/MS CV"
        )
      )
    
    p4 <- ggplot(variation_long, aes(x = Variation_Type, y = Variation_Value, fill = Variation_Type)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
      stat_summary(fun = median, geom = "point", size = 2, color = "white") +
      scale_fill_manual(values = c("#2E86AB", "#A23B72", "#F18F01")) +
      scale_y_log10(labels = scales::scientific_format()) +
      labs(
        title = "Distribution of Variation Metrics",
        subtitle = "Comparison of PEA and MS/MS variability measures",
        x = "Variation Type",
        y = "Variation Value (log scale)"
      ) +
      custom_theme +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # 5. Hexagonal density plot for better visualization of overlap
    p5 <- ggplot(merged_data, aes(x = Variation.between.individuals, y = CV_MSMS_mean)) +
      geom_hex(bins = 20, alpha = 0.8) +
      geom_smooth(method = "lm", color = "#A23B72", se = TRUE, alpha = 0.3) +
      scale_fill_viridis_c(name = "Count", trans = "log10") +
      stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.9, size = 4, color = "white") +
      labs(
        title = "Density Plot: PEA Between-Individual Variation vs MS/MS CV",
        subtitle = "Hexagonal binning shows data density",
        x = "PEA Variation Between Individuals",
        y = "MS/MS Coefficient of Variation (mean)"
      ) +
      custom_theme
    
    # Save individual plots
    ggsave(file.path(output_dir, "pea_between_vs_msms_cv.png"), p1, width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "pea_within_vs_msms_cv.png"), p2, width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "pea_between_vs_within.png"), p3, width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "variation_distributions.png"), p4, width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "density_between_vs_msms.png"), p5, width = 10, height = 8, dpi = 300)
    
    # Create summary panel
    summary_panel <- plot_grid(
      p1, p2, p3, p4,
      ncol = 2, nrow = 2,
      labels = c("A", "B", "C", "D"),
      label_size = 16
    )
    
    ggsave(file.path(output_dir, "pea_msms_variability_summary.png"), summary_panel, 
           width = 16, height = 12, dpi = 300)
    
    # Calculate correlation statistics
    cor_between_msms <- cor.test(merged_data$Variation.between.individuals, merged_data$CV_MSMS_mean, method = "pearson")
    cor_within_msms <- cor.test(merged_data$Variation.within.individuals, merged_data$CV_MSMS_mean, method = "pearson")
    cor_between_within <- cor.test(merged_data$Variation.between.individuals, merged_data$Variation.within.individuals, method = "pearson")
    
    # Generate statistical report
    report <- paste0(
      "HPA PEA vs MS/MS Variability Correlation Analysis Report\n",
      "=" %+% paste(rep("=", 55), collapse = "") %+% "\n\n",
      "Dataset Summary:\n",
      "- HPA PEA proteins: ", nrow(hpa_pea), "\n",
      "- MS/MS datasets processed: ", length(msms_variability), "\n",
      "- Total proteins with MS/MS CV: ", nrow(msms_avg), "\n",
      "- Proteins with both PEA and MS/MS data: ", nrow(merged_data), "\n\n",
      
      "Correlation Results:\n",
      "1. PEA Between-Individual Variation vs MS/MS CV:\n",
      "   - Pearson r = ", round(cor_between_msms$estimate, 3), "\n",
      "   - p-value = ", format(cor_between_msms$p.value, scientific = TRUE, digits = 3), "\n",
      "   - 95% CI: [", round(cor_between_msms$conf.int[1], 3), ", ", round(cor_between_msms$conf.int[2], 3), "]\n\n",
      
      "2. PEA Within-Individual Variation vs MS/MS CV:\n",
      "   - Pearson r = ", round(cor_within_msms$estimate, 3), "\n",
      "   - p-value = ", format(cor_within_msms$p.value, scientific = TRUE, digits = 3), "\n",
      "   - 95% CI: [", round(cor_within_msms$conf.int[1], 3), ", ", round(cor_within_msms$conf.int[2], 3), "]\n\n",
      
      "3. PEA Between vs Within Individual Variation:\n",
      "   - Pearson r = ", round(cor_between_within$estimate, 3), "\n",
      "   - p-value = ", format(cor_between_within$p.value, scientific = TRUE, digits = 3), "\n",
      "   - 95% CI: [", round(cor_between_within$conf.int[1], 3), ", ", round(cor_between_within$conf.int[2], 3), "]\n\n",
      
      "Summary Statistics:\n",
      "PEA Between-Individual Variation:\n",
      "   - Mean: ", round(mean(merged_data$Variation.between.individuals), 3), "\n",
      "   - Median: ", round(median(merged_data$Variation.between.individuals), 3), "\n",
      "   - Range: [", round(min(merged_data$Variation.between.individuals), 3), ", ", 
                    round(max(merged_data$Variation.between.individuals), 3), "]\n\n",
      
      "PEA Within-Individual Variation:\n",
      "   - Mean: ", round(mean(merged_data$Variation.within.individuals), 3), "\n",
      "   - Median: ", round(median(merged_data$Variation.within.individuals), 3), "\n",
      "   - Range: [", round(min(merged_data$Variation.within.individuals), 3), ", ", 
                    round(max(merged_data$Variation.within.individuals), 3), "]\n\n",
      
      "MS/MS Coefficient of Variation:\n",
      "   - Mean: ", round(mean(merged_data$CV_MSMS_mean), 3), "\n",
      "   - Median: ", round(median(merged_data$CV_MSMS_mean), 3), "\n",
      "   - Range: [", round(min(merged_data$CV_MSMS_mean), 3), ", ", 
                    round(max(merged_data$CV_MSMS_mean), 3), "]\n\n",
      
      "Interpretation:\n",
      "- PEA between-individual variation shows ", 
      if(cor_between_msms$estimate > 0) "positive" else "negative", 
      " correlation with MS/MS CV (r=", round(cor_between_msms$estimate, 3), ")\n",
      "- PEA within-individual variation shows ", 
      if(cor_within_msms$estimate > 0) "positive" else "negative", 
      " correlation with MS/MS CV (r=", round(cor_within_msms$estimate, 3), ")\n",
      "- Between-individual variation is generally ", 
      if(mean(merged_data$Variation.between.individuals) > mean(merged_data$Variation.within.individuals)) "higher" else "lower",
      " than within-individual variation\n\n",
      
      "Files Generated:\n",
      "- pea_between_vs_msms_cv.png: Correlation plot (between-individual vs MS/MS)\n",
      "- pea_within_vs_msms_cv.png: Correlation plot (within-individual vs MS/MS)\n",
      "- pea_between_vs_within.png: PEA variation comparison\n",
      "- variation_distributions.png: Distribution comparison\n",
      "- density_between_vs_msms.png: Hexagonal density plot\n",
      "- pea_msms_variability_summary.png: Four-panel summary figure\n"
    )
    
    # Define %+% operator for string concatenation
    `%+%` <- function(a, b) paste0(a, b)
    
    writeLines(report, file.path(output_dir, "correlation_analysis_report.txt"))
    
    # Save merged data for further analysis
    write.csv(merged_data, file.path(output_dir, "merged_pea_msms_data.csv"), row.names = FALSE)
    
    cat("Analysis complete! Files saved to:", output_dir, "\n")
    cat("Found", nrow(merged_data), "proteins with both PEA and MS/MS variability data\n")
    cat("Correlation between PEA (between-individual) and MS/MS CV: r =", round(cor_between_msms$estimate, 3), "\n")
    cat("Correlation between PEA (within-individual) and MS/MS CV: r =", round(cor_within_msms$estimate, 3), "\n")
    
  } else {
    cat("No overlapping proteins found between PEA and MS/MS datasets\n")
  }
  
} else {
  cat("No MS/MS datasets with sufficient variability data found\n")
}

cat("Script completed.\n") 
