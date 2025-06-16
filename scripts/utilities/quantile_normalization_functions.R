# ============================================================================
# QUANTILE NORMALIZATION FUNCTIONS FOR PROTEOMICS DATA
# ============================================================================
# Description: Implementation of quantile normalization methods for 
#              cross-database comparison of proteomics data
# Author: Generated for blood-review-data analysis
# ============================================================================

library(dplyr)
library(tidyr)
library(tibble)

# Load required packages
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  message("Warning: preprocessCore not available. Some functions will use base R alternatives.")
}

#' Manual quantile normalization implementation
#' 
#' @param mat Matrix with rows = genes, columns = samples/databases
#' @return Quantile normalized matrix
#' 
normalize_quantiles_manual <- function(mat) {
  
  # Validate input
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("Input must be a matrix or data.frame")
  }
  
  # Convert to matrix if data.frame
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }
  
  # Check if matrix contains atomic values
  if (!is.numeric(mat)) {
    stop("Matrix must contain numeric values only")
  }
  
  # Handle missing values
  na_mask <- is.na(mat)
  
  # For each column, compute the sorted values and their ranks
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  
  if (n_cols == 0 || n_rows == 0) {
    warning("Empty matrix provided - returning as-is")
    return(mat)
  }
  
  # Calculate the reference distribution (average of all sorted columns)
  # Add error handling for the sort operation
  sorted_cols <- list()
  for (i in 1:n_cols) {
    col_data <- mat[, i]
    non_na_data <- col_data[!is.na(col_data)]
    
    # Check if column data is atomic
    if (!is.atomic(non_na_data)) {
      stop(paste("Column", i, "contains non-atomic data. Check for list-columns in input."))
    }
    
    if (length(non_na_data) > 0) {
      sorted_cols[[i]] <- sort(non_na_data)
    } else {
      sorted_cols[[i]] <- numeric(0)
    }
  }
  
  # Handle columns with different numbers of non-NA values
  col_lengths <- sapply(sorted_cols, length)
  if (all(col_lengths == 0)) {
    warning("All columns contain only NA values - returning original matrix")
    return(mat)
  }
  
  max_length <- max(col_lengths)
  
  # Pad shorter columns with NA at the end for alignment
  padded_cols <- lapply(sorted_cols, function(x) {
    if (length(x) < max_length) {
      c(x, rep(NA, max_length - length(x)))
    } else {
      x
    }
  })
  
  # Calculate reference distribution (mean across columns at each quantile)
  padded_matrix <- do.call(cbind, padded_cols)
  reference_dist <- rowMeans(padded_matrix, na.rm = TRUE)
  
  # Apply quantile normalization to each column
  normalized_mat <- mat
  for (j in 1:n_cols) {
    col_data <- mat[, j]
    non_na <- !is.na(col_data)
    
    if (sum(non_na) > 0) {
      # Get ranks of non-NA values
      ranks <- rank(col_data[non_na], ties.method = "average")
      
      # Map ranks to reference distribution
      # Scale ranks to range [1, length(reference_dist)]
      if (length(ranks) > 1) {
        scaled_ranks <- (ranks - 1) * (length(reference_dist) - 1) / (length(ranks) - 1) + 1
      } else {
        scaled_ranks <- length(reference_dist) / 2  # Single value goes to median
      }
      
      # Interpolate to get normalized values
      if (length(reference_dist) >= 2) {
        normalized_values <- approx(x = 1:length(reference_dist), 
                                   y = reference_dist, 
                                   xout = scaled_ranks, 
                                   method = "linear")$y
      } else {
        # Fallback for single value case
        normalized_values <- rep(reference_dist[1], length(scaled_ranks))
      }
      
      normalized_mat[non_na, j] <- normalized_values
    }
  }
  
  return(normalized_mat)
}

#' Quantile Normalization using preprocessCore (Bioconductor)
#' 
#' @param data Data frame with abundance data
#' @param abundance_col Column name for abundance values
#' @param group_col Column name for grouping (e.g., "source", "database")
#' @param log_transform Whether to log-transform data before normalization
#' @return Data frame with quantile_normalized column added
#' 
apply_quantile_normalization_bioc <- function(data, abundance_col = "abundance", 
                                             group_col = "source", log_transform = TRUE) {
  
  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("preprocessCore package is required for this function. Install with BiocManager::install('preprocessCore')")
  }
  
  # Prepare data
  if (log_transform) {
    data <- data %>%
      mutate(log_abundance = log10(.data[[abundance_col]] + 1))
    value_col <- "log_abundance"
  } else {
    value_col <- abundance_col
  }
  
  # Create matrix for quantile normalization
  # Rows = genes, Columns = databases/sources
  wide_data <- data %>%
    select(gene, !!sym(group_col), !!sym(value_col)) %>%
    tidyr::pivot_wider(names_from = !!sym(group_col), values_from = !!sym(value_col), 
                       values_fill = NA) %>%
    column_to_rownames("gene")
  
  # Apply quantile normalization
  normalized_matrix <- preprocessCore::normalize.quantiles(as.matrix(wide_data), copy = TRUE)
  rownames(normalized_matrix) <- rownames(wide_data)
  colnames(normalized_matrix) <- colnames(wide_data)
  
  # Convert back to long format
  normalized_long <- normalized_matrix %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    tidyr::pivot_longer(cols = -gene, names_to = group_col, values_to = "quantile_normalized") %>%
    filter(!is.na(quantile_normalized))
  
  # Merge back with original data
  result <- data %>%
    left_join(normalized_long, by = c("gene", group_col))
  
  return(result)
}

#' Quantile Normalization using base R (no Bioconductor dependency)
#' 
#' @param data Data frame with abundance data
#' @param abundance_col Column name for abundance values
#' @param group_col Column name for grouping (e.g., "source", "database")
#' @param log_transform Whether to log-transform data before normalization
#' @return Data frame with quantile_normalized column added
#' 
apply_quantile_normalization_base <- function(data, abundance_col = "abundance", 
                                             group_col = "source", log_transform = TRUE) {
  
  # Prepare data
  if (log_transform) {
    data <- data %>%
      mutate(log_abundance = log10(.data[[abundance_col]] + 1))
    value_col <- "log_abundance"
  } else {
    value_col <- abundance_col
  }
  
  # Create matrix for quantile normalization
  wide_data <- data %>%
    select(gene, !!sym(group_col), !!sym(value_col)) %>%
    tidyr::pivot_wider(names_from = !!sym(group_col), values_from = !!sym(value_col), 
                       values_fill = NA) %>%
    column_to_rownames("gene")
  
  # Apply quantile normalization (base R implementation)
  normalized_matrix <- normalize_quantiles_manual(as.matrix(wide_data))
  rownames(normalized_matrix) <- rownames(wide_data)
  colnames(normalized_matrix) <- colnames(wide_data)
  
  # Convert back to long format
  normalized_long <- normalized_matrix %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    tidyr::pivot_longer(cols = -gene, names_to = group_col, values_to = "quantile_normalized") %>%
    filter(!is.na(quantile_normalized))
  
  # Merge back with original data
  result <- data %>%
    left_join(normalized_long, by = c("gene", group_col))
  
  return(result)
}

#' Apply quantile normalization within each database/source individually
#' Similar to z-score normalization but using quantile normalization within each group
#' 
#' @param data Data frame with abundance data
#' @param abundance_col Column name for abundance values
#' @param group_col Column name for grouping (e.g., "source", "database")
#' @param log_transform Whether to log-transform data before normalization
#' @return Data frame with quantile_normalized_within column added
#' 
apply_quantile_normalization_within <- function(data, abundance_col = "abundance", 
                                               group_col = "source", log_transform = TRUE) {
  
  # Prepare data
  if (log_transform) {
    data <- data %>%
      mutate(log_abundance = log10(.data[[abundance_col]] + 1))
    value_col <- "log_abundance"
  } else {
    value_col <- abundance_col
  }
  
  # Apply quantile normalization within each group separately
  result <- data %>%
    group_by(.data[[group_col]]) %>%
    mutate(
      quantile_normalized_within = {
        values <- .data[[value_col]]
        non_na_values <- values[!is.na(values)]
        
        if (length(non_na_values) <= 1) {
          # If only one value or all NA, return as-is
          values
        } else {
          # Create a reference distribution from this group's values
          sorted_values <- sort(non_na_values)
          
          # For each value, find its rank and map to the reference distribution
          ranks <- rank(values, ties.method = "average", na.last = "keep")
          
          # Map ranks to evenly spaced quantiles (0 to 1)
          quantiles <- (ranks - 1) / (length(non_na_values) - 1)
          
          # Map quantiles back to the reference distribution
          # Use the empirical distribution of this group as reference
          reference_quantiles <- seq(0, 1, length.out = length(sorted_values))
          
          # Interpolate to get normalized values
          approx(x = reference_quantiles, y = sorted_values, 
                xout = quantiles, method = "linear", rule = 2)$y
        }
      }
    ) %>%
    ungroup()
  
  return(result)
}

#' Compare normalization methods
#' 
#' @param data Data frame with abundance data
#' @param abundance_col Column name for abundance values
#' @param group_col Column name for grouping
#' @return Data frame with all normalization methods applied
#' 
compare_normalization_methods <- function(data, abundance_col = "abundance", group_col = "source") {
  
  message("Applying different normalization methods...")
  
  # 1. Z-score normalization (existing method - within databases)
  data_zscore <- apply_zscore_normalization(data, abundance_col, group_col)
  
  # 2. Quantile normalization within databases (new approach)
  data_quantile_within <- apply_quantile_normalization_within(data, abundance_col, group_col)
  
  # 3. Quantile normalization across databases (integration approach)
  data_quantile_across <- apply_quantile_normalization_base(data, abundance_col, group_col)
  
  # 4. Combine all results
  result <- data_zscore %>%
    left_join(
      data_quantile_within %>% select(gene, !!sym(group_col), quantile_normalized_within),
      by = c("gene", group_col)
    ) %>%
    left_join(
      data_quantile_across %>% select(gene, !!sym(group_col), quantile_normalized),
      by = c("gene", group_col)
    ) %>%
    # Rename for clarity
    rename(
      quantile_normalized_across = quantile_normalized
    )
  
  return(result)
}

#' Apply quantile normalization to simple data format
#' 
#' @param data Data frame with gene and value columns
#' @param value_col Column name for values to normalize
#' @param group_col Column name for grouping
#' @return Data frame with quantile_normalized column added
#' 
apply_quantile_normalization_simple <- function(data, value_col = "log_value", group_col = "Database") {
  
  # Handle cases where we only have one group (no cross-group normalization possible)
  unique_groups <- unique(data[[group_col]])
  if (length(unique_groups) <= 1) {
    message("Only one group found - returning original data with quantile_normalized = value")
    data$quantile_normalized <- data[[value_col]]
    return(data)
  }
  
  # Check if data has gene column - if not, use row-based approach
  if (!"gene" %in% colnames(data)) {
    return(apply_quantile_normalization_no_genes(data, value_col, group_col))
  }
  
  # Handle duplicate gene-group combinations by aggregating (taking mean)
  # This prevents pivot_wider from creating list-columns
  aggregated_data <- data %>%
    select(gene, !!sym(group_col), !!sym(value_col)) %>%
    filter(!is.na(!!sym(value_col))) %>%  # Remove NA values before aggregation
    group_by(gene, !!sym(group_col)) %>%
    summarise(
      aggregated_value = mean(!!sym(value_col), na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create matrix for quantile normalization
  wide_data <- aggregated_data %>%
    pivot_wider(names_from = !!sym(group_col), values_from = aggregated_value, 
                values_fill = NA)
  
  # Extract gene names and create matrix
  gene_names <- wide_data$gene
  value_matrix <- as.matrix(wide_data[, -1])
  rownames(value_matrix) <- gene_names
  
  # Apply quantile normalization
  normalized_matrix <- normalize_quantiles_manual(value_matrix)
  
  # Convert back to long format
  normalized_long <- normalized_matrix %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = group_col, values_to = "quantile_normalized") %>%
    filter(!is.na(quantile_normalized))
  
  # Merge back with original data - handling duplicates properly
  # First, add the aggregated values to the original data for matching
  data_with_agg <- data %>%
    left_join(aggregated_data %>% rename(!!paste0(value_col, "_agg") := aggregated_value), 
              by = c("gene", group_col))
  
  # Then add the quantile normalized values
  result <- data_with_agg %>%
    left_join(normalized_long, by = c("gene", group_col)) %>%
    select(-ends_with("_agg"))  # Remove the temporary aggregated column
  
  return(result)
}

#' Apply quantile normalization when no gene identifiers are available
#' 
#' @param data Data frame with value and group columns (no gene column)
#' @param value_col Column name for values to normalize
#' @param group_col Column name for grouping
#' @return Data frame with quantile_normalized column added
#' 
apply_quantile_normalization_no_genes <- function(data, value_col = "log_value", group_col = "Database") {
  
  # Handle cases where we only have one group
  unique_groups <- unique(data[[group_col]])
  if (length(unique_groups) <= 1) {
    message("Only one group found - returning original data with quantile_normalized = value")
    data$quantile_normalized <- data[[value_col]]
    return(data)
  }
  
  # For data without gene identifiers, we'll apply a simpler approach
  # that normalizes the distribution shapes while preserving within-group rankings
  
  result_data <- data
  result_data$quantile_normalized <- NA
  
  # Get all non-NA values across all groups to create reference distribution
  all_values <- data[[value_col]][!is.na(data[[value_col]])]
  
  if (length(all_values) == 0) {
    message("No valid values found - returning original data")
    result_data$quantile_normalized <- data[[value_col]]
    return(result_data)
  }
  
  # Create reference distribution from all values
  reference_quantiles <- quantile(all_values, probs = seq(0, 1, length.out = 1000), na.rm = TRUE)
  
  # Apply quantile normalization to each group
  for (group_name in unique_groups) {
    group_indices <- which(data[[group_col]] == group_name & !is.na(data[[value_col]]))
    group_values <- data[[value_col]][group_indices]
    
    if (length(group_values) > 0) {
      # Get ranks of values within this group
      group_ranks <- rank(group_values, ties.method = "average")
      
      # Convert ranks to quantiles (0 to 1)
      group_quantiles <- (group_ranks - 1) / (length(group_ranks) - 1)
      
      # Handle edge case where all values are the same
      if (length(group_ranks) == 1) {
        group_quantiles <- 0.5  # Use median
      }
      
      # Map quantiles to reference distribution
      normalized_values <- approx(x = seq(0, 1, length.out = length(reference_quantiles)), 
                                 y = reference_quantiles, 
                                 xout = group_quantiles, 
                                 method = "linear")$y
      
      result_data$quantile_normalized[group_indices] <- normalized_values
    }
  }
  
  return(result_data)
}

#' Evaluate normalization effectiveness
#' 
#' @param data Data frame with normalized data
#' @param group_col Column name for grouping
#' @return Data frame with normalization assessment metrics
#' 
evaluate_normalization_effectiveness <- function(data, group_col = "source") {
  
  # Calculate metrics for each normalization method
  metrics <- data %>%
    group_by(!!sym(group_col)) %>%
    summarise(
      n_proteins = n(),
      
      # Original log values
      log_mean = mean(log_abundance, na.rm = TRUE),
      log_median = median(log_abundance, na.rm = TRUE),
      log_sd = sd(log_abundance, na.rm = TRUE),
      log_cv = log_sd / abs(log_mean), # Coefficient of variation
      
      # Z-score normalized
      zscore_mean = mean(z_score, na.rm = TRUE),
      zscore_median = median(z_score, na.rm = TRUE),
      zscore_sd = sd(z_score, na.rm = TRUE),
      zscore_range = max(z_score, na.rm = TRUE) - min(z_score, na.rm = TRUE),
      
      # Quantile normalized within databases
      quantile_within_mean = mean(quantile_normalized_within, na.rm = TRUE),
      quantile_within_median = median(quantile_normalized_within, na.rm = TRUE),
      quantile_within_sd = sd(quantile_normalized_within, na.rm = TRUE),
      quantile_within_range = max(quantile_normalized_within, na.rm = TRUE) - min(quantile_normalized_within, na.rm = TRUE),
      
      # Quantile normalized across databases
      quantile_across_mean = mean(quantile_normalized_across, na.rm = TRUE),
      quantile_across_median = median(quantile_normalized_across, na.rm = TRUE),
      quantile_across_sd = sd(quantile_normalized_across, na.rm = TRUE),
      quantile_across_range = max(quantile_normalized_across, na.rm = TRUE) - min(quantile_normalized_across, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  # Calculate cross-database consistency metrics
  consistency_metrics <- data.frame(
    metric = c("mean_sd_zscore", "mean_sd_quantile_within", "mean_sd_quantile_across", 
               "median_sd_zscore", "median_sd_quantile_within", "median_sd_quantile_across"),
    value = c(
      sd(metrics$zscore_mean),              # Variability in means across databases (Z-score)
      sd(metrics$quantile_within_mean),     # Variability in means across databases (Quantile within)
      sd(metrics$quantile_across_mean),     # Variability in means across databases (Quantile across)
      sd(metrics$zscore_median),            # Variability in medians across databases (Z-score)
      sd(metrics$quantile_within_median),   # Variability in medians across databases (Quantile within)
      sd(metrics$quantile_across_median)    # Variability in medians across databases (Quantile across)
    ),
    interpretation = c(
      "Lower is better - indicates more consistent means across databases",
      "Lower is better - indicates more consistent means across databases",
      "Lower is better - indicates more consistent means across databases",
      "Lower is better - indicates more consistent medians across databases",
      "Lower is better - indicates more consistent medians across databases",
      "Lower is better - indicates more consistent medians across databases"
    )
  )
  
  return(list(
    by_database = metrics,
    consistency = consistency_metrics
  ))
}

#' Create normalization comparison plots
#' 
#' @param data Data frame with all normalization methods applied
#' @param group_col Column name for grouping
#' @param output_dir Output directory for plots
#' 
create_normalization_comparison_plots <- function(data, group_col = "source", output_dir) {
  
  # 1. Distribution comparison plot
  plot_data_long <- data %>%
    select(gene, !!sym(group_col), log_abundance, z_score, quantile_normalized_within, quantile_normalized_across) %>%
    pivot_longer(cols = c(log_abundance, z_score, quantile_normalized_within, quantile_normalized_across),
                names_to = "normalization_method",
                values_to = "value") %>%
    mutate(
      normalization_method = case_when(
        normalization_method == "log_abundance" ~ "Log10 Original",
        normalization_method == "z_score" ~ "Z-Score Normalized",
        normalization_method == "quantile_normalized_within" ~ "Quantile Within Database",
        normalization_method == "quantile_normalized_across" ~ "Quantile Across Databases"
      )
    )
  
  # Distribution comparison
  p_dist_compare <- ggplot(plot_data_long, aes(x = value, fill = !!sym(group_col))) +
    geom_density(alpha = 0.6) +
    facet_grid(normalization_method ~ ., scales = "free") +
    scale_fill_manual(values = get_plot_colors("databases")) +
    theme_blood_proteomics() +
    labs(
      title = "Comparison of Normalization Methods",
      subtitle = "Distribution of protein abundances across databases",
      x = "Normalized Value",
      y = "Density",
      fill = "Database"
    )
  
  save_plot_standard(p_dist_compare, "normalization_methods_comparison", output_dir,
                    width = 12, height = 10)
  
  # 2. Box plot comparison
  p_box_compare <- ggplot(plot_data_long, aes(x = !!sym(group_col), y = value, fill = !!sym(group_col))) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
    facet_wrap(~ normalization_method, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = get_plot_colors("databases")) +
    theme_blood_proteomics() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Database Distribution Comparison Across Normalization Methods",
      subtitle = "Box plots showing distribution characteristics",
      x = "Database",
      y = "Normalized Value",
      fill = "Database"
    )
  
  save_plot_standard(p_box_compare, "normalization_boxplot_comparison", output_dir,
                    width = 10, height = 12)
  

}

#' Generate comprehensive normalization analysis report
#' 
#' @param data Data frame with abundance data
#' @param abundance_col Column name for abundance values
#' @param group_col Column name for grouping
#' @param output_dir Output directory
#' @return List with analysis results
#' 
run_quantile_normalization_analysis <- function(data, abundance_col = "abundance", 
                                               group_col = "source", output_dir) {
  
  message(paste(rep("=", 60), collapse=""))
  message("QUANTILE NORMALIZATION COMPREHENSIVE ANALYSIS")
  message(paste(rep("=", 60), collapse=""))
  
  # 1. Apply all normalization methods
  message("\n[STEP 1] Applying normalization methods...")
  normalized_data <- compare_normalization_methods(data, abundance_col, group_col)
  
  # 2. Evaluate effectiveness
  message("\n[STEP 2] Evaluating normalization effectiveness...")
  effectiveness <- evaluate_normalization_effectiveness(normalized_data, group_col)
  
  # 3. Create visualizations
  message("\n[STEP 3] Creating comparison visualizations...")
  create_normalization_comparison_plots(normalized_data, group_col, output_dir)
  
  # 4. Save results
  message("\n[STEP 4] Saving analysis results...")
  
  # Save effectiveness metrics
  write_csv(effectiveness$by_database, 
           file.path(output_dir, "normalization_effectiveness_by_database.csv"))
  write_csv(effectiveness$consistency, 
           file.path(output_dir, "normalization_consistency_metrics.csv"))
  
  # Save normalized data
  write_csv(normalized_data, 
           file.path(output_dir, "all_normalization_methods_data.csv"))
  
  # 5. Generate summary report
  generate_normalization_summary_report(effectiveness, output_dir)
  
  message("\n✅ Quantile normalization analysis completed successfully!")
  
  return(list(
    data = normalized_data,
    effectiveness = effectiveness
  ))
}

#' Generate normalization summary report
#' 
#' @param effectiveness List with effectiveness metrics
#' @param output_dir Output directory
#' 
generate_normalization_summary_report <- function(effectiveness, output_dir) {
  
  report_file <- file.path(output_dir, "quantile_normalization_analysis_report.txt")
  
  sink(report_file)
  
  cat("QUANTILE NORMALIZATION ANALYSIS REPORT\n")
  cat("=====================================\n\n")
  
  cat("ANALYSIS OVERVIEW:\n")
  cat("This analysis compares three normalization methods:\n")
  cat("1. Log10 transformation (baseline)\n")
  cat("2. Z-score normalization (current method)\n")
  cat("3. Quantile normalization (proposed method)\n\n")
  
  cat("NORMALIZATION EFFECTIVENESS BY DATABASE:\n")
  cat("=======================================\n")
  print(effectiveness$by_database)
  cat("\n")
  
  cat("CROSS-DATABASE CONSISTENCY METRICS:\n")
  cat("==================================\n")
  print(effectiveness$consistency)
  cat("\n")
  
  cat("INTERPRETATION:\n")
  cat("==============\n")
  cat("- Z-Score Normalization: Standardizes each database to mean=0, sd=1\n")
  cat("- Quantile Normalization: Forces all databases to have identical distributions\n")
  cat("- Lower consistency metric values indicate better cross-database normalization\n\n")
  
  # Determine best method
  zscore_consistency <- mean(effectiveness$consistency$value[c(1,3)])
  quantile_consistency <- mean(effectiveness$consistency$value[c(2,4)])
  
  cat("RECOMMENDATION:\n")
  cat("==============\n")
  if (quantile_consistency < zscore_consistency) {
    cat("✅ QUANTILE NORMALIZATION is recommended\n")
    cat(sprintf("   - Consistency score: %.4f (vs %.4f for Z-score)\n", 
               quantile_consistency, zscore_consistency))
    cat("   - Better cross-database harmonization\n")
    cat("   - More robust to outliers and distributional differences\n")
  } else {
    cat("✅ Z-SCORE NORMALIZATION remains adequate\n")
    cat(sprintf("   - Consistency score: %.4f (vs %.4f for Quantile)\n", 
               zscore_consistency, quantile_consistency))
    cat("   - Simpler implementation and interpretation\n")
    cat("   - Preserves within-database relative relationships\n")
  }
  
  cat("\nFILES GENERATED:\n")
  cat("===============\n")
  cat("• normalization_effectiveness_by_database.csv - Detailed metrics by database\n")
  cat("• normalization_consistency_metrics.csv - Cross-database consistency scores\n")
  cat("• all_normalization_methods_data.csv - Complete normalized dataset\n")
  cat("• normalization_methods_comparison.png - Distribution comparison plots\n")
  cat("• normalization_boxplot_comparison.png - Box plot comparisons\n")

  cat("• quantile_normalization_analysis_report.txt - This report\n")
  
  sink()
  
  message(sprintf("Summary report saved to: %s", report_file))
} 