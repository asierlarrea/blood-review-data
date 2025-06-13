#!/usr/bin/env Rscript

# Gene Deduplication Utility
# Purpose: Handle cases where multiple proteins map to the same gene
# Removes duplicates and uses median quantification values

suppressPackageStartupMessages({
  library(dplyr)
})

#' Deduplicate genes and aggregate quantification values
#' 
#' This function handles the common issue where multiple proteins map to the same gene.
#' It removes duplicate genes and uses the median quantification value when multiple
#' proteins map to the same gene.
#' 
#' @param data Data frame containing protein data
#' @param gene_col Character string specifying the gene column name
#' @param quant_col Character string specifying the quantification column name
#' @param additional_cols Character vector of additional columns to preserve (optional)
#' @param aggregation_method Character string specifying aggregation method ("median", "mean", "max", "min")
#' @return Data frame with deduplicated genes and aggregated quantification values
#' 
deduplicate_genes <- function(data, gene_col, quant_col, additional_cols = NULL, aggregation_method = "median") {
  
  # Validate inputs
  if (!gene_col %in% colnames(data)) {
    stop(sprintf("Gene column '%s' not found in data", gene_col))
  }
  if (!quant_col %in% colnames(data)) {
    stop(sprintf("Quantification column '%s' not found in data", quant_col))
  }
  if (!aggregation_method %in% c("median", "mean", "max", "min")) {
    stop("aggregation_method must be one of: 'median', 'mean', 'max', 'min'")
  }
  
  # Remove rows with missing gene or quantification values
  clean_data <- data %>%
    filter(!is.na(.data[[gene_col]]) & !is.na(.data[[quant_col]]) & 
           .data[[gene_col]] != "" & .data[[quant_col]] > 0)
  
  original_count <- nrow(clean_data)
  
  # Prepare aggregation formula
  agg_func <- switch(aggregation_method,
    "median" = function(x) median(x, na.rm = TRUE),
    "mean" = function(x) mean(x, na.rm = TRUE),
    "max" = function(x) max(x, na.rm = TRUE),
    "min" = function(x) min(x, na.rm = TRUE)
  )
  
  # Group by gene and aggregate
  if (is.null(additional_cols)) {
    # Simple case: only gene and quantification columns
    result <- clean_data %>%
      group_by(.data[[gene_col]]) %>%
      summarise(
        !!quant_col := agg_func(.data[[quant_col]]),
        protein_count = n(),
        .groups = "drop"
      )
  } else {
    # Complex case: preserve additional columns
    # For additional columns, we'll take the first value (or could use most common)
    group_cols <- c(gene_col, additional_cols)
    
    # First, get the main aggregation
    main_agg <- clean_data %>%
      group_by(.data[[gene_col]]) %>%
      summarise(
        !!quant_col := agg_func(.data[[quant_col]]),
        protein_count = n(),
        .groups = "drop"
      )
    
    # Then get the first occurrence of additional columns for each gene
    additional_data <- clean_data %>%
      group_by(.data[[gene_col]]) %>%
      slice(1) %>%
      select(all_of(c(gene_col, additional_cols))) %>%
      ungroup()
    
    # Merge the results
    result <- main_agg %>%
      left_join(additional_data, by = gene_col)
  }
  
  final_count <- nrow(result)
  duplicated_genes <- original_count - final_count
  
  # Report summary
  message(sprintf("Gene deduplication summary:"))
  message(sprintf("  - Original entries: %d", original_count))
  message(sprintf("  - Unique genes: %d", final_count))
  message(sprintf("  - Duplicated entries removed: %d", duplicated_genes))
  message(sprintf("  - Aggregation method: %s", aggregation_method))
  
  if (duplicated_genes > 0) {
    message(sprintf("  - %.1f%% of entries had gene duplicates", (duplicated_genes/original_count)*100))
  }
  
  return(result)
}

#' Quick function to get gene count statistics before and after deduplication
#' 
#' @param data Data frame containing protein data
#' @param gene_col Character string specifying the gene column name
#' @param quant_col Character string specifying the quantification column name (optional, for filtering)
#' @return List with before/after statistics
#' 
get_deduplication_stats <- function(data, gene_col, quant_col = NULL) {
  
  # Filter valid data
  if (!is.null(quant_col)) {
    clean_data <- data %>%
      filter(!is.na(.data[[gene_col]]) & !is.na(.data[[quant_col]]) & 
             .data[[gene_col]] != "" & .data[[quant_col]] > 0)
  } else {
    clean_data <- data %>%
      filter(!is.na(.data[[gene_col]]) & .data[[gene_col]] != "")
  }
  
  # Count duplicates
  gene_counts <- clean_data %>%
    group_by(.data[[gene_col]]) %>%
    summarise(protein_count = n(), .groups = "drop")
  
  stats <- list(
    total_entries = nrow(clean_data),
    unique_genes = nrow(gene_counts),
    genes_with_duplicates = sum(gene_counts$protein_count > 1),
    max_proteins_per_gene = max(gene_counts$protein_count),
    avg_proteins_per_gene = mean(gene_counts$protein_count)
  )
  
  return(stats)
}

#' Print formatted deduplication statistics
#' 
#' @param stats Statistics list from get_deduplication_stats()
#' @param title Optional title for the statistics display
#' 
print_deduplication_stats <- function(stats, title = "Gene Deduplication Statistics") {
  cat("\n", title, "\n")
  cat(rep("=", nchar(title)), "\n", sep = "")
  cat(sprintf("Total entries: %d\n", stats$total_entries))
  cat(sprintf("Unique genes: %d\n", stats$unique_genes))
  cat(sprintf("Genes with duplicates: %d\n", stats$genes_with_duplicates))
  cat(sprintf("Max proteins per gene: %d\n", stats$max_proteins_per_gene))
  cat(sprintf("Average proteins per gene: %.2f\n", stats$avg_proteins_per_gene))
  
  if (stats$total_entries > stats$unique_genes) {
    reduction_pct <- ((stats$total_entries - stats$unique_genes) / stats$total_entries) * 100
    cat(sprintf("Reduction in count: %.1f%%\n", reduction_pct))
  }
  cat("\n")
} 