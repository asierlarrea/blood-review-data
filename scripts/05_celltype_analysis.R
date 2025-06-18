#!/usr/bin/env Rscript
#
# Analysis of Cell Type Protein Expression Data
# Script: 05_celltype_analysis.R
# 
# Purpose: Analyze protein expression across different blood cell types from multiple sources:
# - PAXDB (cell type in filename: paxdb_b_cell.csv, paxdb_cd4.csv, etc.)
# - ProteomeXchange (cell type in columns: Intensity_B.memory_03_activated, etc.)
# 
# Key features:
# - Filter out proteins with multiple IDs (e.g., Q9C0C7;A0A075B6T1;Q9C0C7-4)
# - Map proteins to gene names using existing utilities
# - Calculate median intensities for multiple columns of the same cell type
# - Handle different data formats consistently

# Load utilities and set up output directories
source("scripts/utilities/load_packages.R")
source("scripts/config/analysis_config.R")
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", 
                      "ggthemes", "patchwork", "UpSetR", "ggupset", "RColorBrewer", "ggrepel", "tibble")
load_packages(required_packages)

# Ensure patchwork is explicitly loaded
suppressPackageStartupMessages(library(patchwork))

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

# Load gene mapping and deduplication utilities
source("scripts/data_processing/simple_id_mapping.R")
source("scripts/utilities/gene_deduplication.R")

# Load quantile normalization utility
source("scripts/utilities/quantile_normalization_functions.R")

# Load specialized processors
source("scripts/data_processing/proteomexchange_processor.R")
source("scripts/data_processing/simple_celltype_processor.R")

# Set output directory
output_dir <- get_output_path("", subdir = "celltype_analysis")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Helper function to check if protein ID is uniquely identified
is_unique_protein <- function(protein_ids) {
  # Check if protein ID contains semicolons (multiple IDs)
  return(!str_detect(protein_ids, ";"))
}

# Helper function to safely use %||% operator (null-coalescing)
`%||%` <- function(a, b) if (is.null(a)) b else a

message("=== Cell Type Analysis Starting ===")

# 1. Process PAXDB files (cell type in filename)
message("\n1. Processing PAXDB files...")

paxdb_files <- list.files("data/raw/paxdb", pattern = "^paxdb_.*\\.csv$", full.names = TRUE)
# Filter out plasma and serum files (not cell types)
paxdb_files <- paxdb_files[!str_detect(basename(paxdb_files), "(plasma|serum)")]
message(sprintf("  Found %d cell type files (excluding plasma/serum)", length(paxdb_files)))

paxdb_results <- list()

for (file_path in paxdb_files) {
  # Use the specialized processor
  result <- process_filename_celltype_file(
    file_path = file_path,
    celltype = NULL,  # Auto-detect from filename
    intensity_column = "abundance",
    id_column = "string_external_id",
    force_mapping = force_mapping
  )
  
  if (!is.null(result)) {
    celltype <- unique(result$celltype)[1]
    paxdb_results[[celltype]] <- result
  }
}

# 2. Process GPMDB files (cell type in filename)
message("\n2. Processing GPMDB files...")

gpmdb_files <- list.files("data/raw/gpmdb", pattern = "^gpmdb_.*\\.csv$", full.names = TRUE)
# Filter out plasma and serum files (not cell types)
gpmdb_files <- gpmdb_files[!str_detect(basename(gpmdb_files), "(plasma|serum)")]
message(sprintf("  Found %d cell type files (excluding plasma/serum)", length(gpmdb_files)))

gpmdb_results <- list()

for (file_path in gpmdb_files) {
  # Use the specialized processor
  result <- process_filename_celltype_file(
    file_path = file_path,
    celltype = NULL,  # Auto-detect from filename
    intensity_column = "total",  # GPMDB uses 'total' column
    id_column = "accession",  # GPMDB uses 'accession' column for protein IDs
    force_mapping = force_mapping
  )
  
  if (!is.null(result)) {
    celltype <- unique(result$celltype)[1]
    gpmdb_results[[celltype]] <- result
  }
}

# 3. Process ProteomeXchange files (cell type in columns)
message("\n3. Processing ProteomeXchange files...")

px_results <- list()

# Process pxd004352.csv (the large file with intensity columns)
px_file <- "data/raw/proteomexchange/pxd004352.csv"

if (file.exists(px_file)) {
  # Use the specialized processor for complex intensity files
  px_result <- process_proteomexchange_intensity_file(
    file_path = px_file,
    force_mapping = force_mapping
  )
  
  if (!is.null(px_result)) {
    # Split results by cell type
    for (celltype in unique(px_result$celltype)) {
      ct_data <- px_result %>% filter(celltype == !!celltype)
      px_results[[celltype]] <- ct_data
    }
  }
} else {
  message("  pxd004352.csv not found, skipping...")
}

# Process other ProteomeXchange files
other_px_files <- c(
  "data/raw/proteomexchange/pxd025174.csv",
  "data/raw/proteomexchange/pxd040957_cd8.csv", 
  "data/raw/proteomexchange/pxd040957_macrophages.csv"
)

for (file_path in other_px_files) {
  if (file.exists(file_path)) {
    filename <- basename(file_path)
    message(sprintf("  Processing %s", filename))
    
    # Use the simple processor for these files
    result <- process_simple_celltype_file(
      file_path = file_path,
      celltype_column_mapping = NULL,  # Auto-detect columns
      force_mapping = force_mapping
    )
    
    if (!is.null(result)) {
      # Split results by cell type and combine with existing results
      for (celltype in unique(result$celltype)) {
        ct_data <- result %>% filter(celltype == !!celltype)
        
        if (celltype %in% names(px_results)) {
          # Combine with existing data
          px_results[[celltype]] <- bind_rows(px_results[[celltype]], ct_data)
        } else {
          px_results[[celltype]] <- ct_data
        }
      }
    }
  } else {
    message(sprintf("  %s not found, skipping...", basename(file_path)))
  }
}

# 4. Combine and summarize results
message("\n4. Combining and summarizing results...")

# Ensure all results have the same column structure
standardize_results <- function(results_list) {
  if (length(results_list) == 0) return(NULL)
  
  # Combine all results first
  combined <- bind_rows(results_list)
  
  # Ensure all required columns exist
  required_cols <- c("gene", "intensity", "protein_count", "celltype", "source", "filename")
  
  for (col in required_cols) {
    if (!col %in% colnames(combined)) {
      if (col == "protein_count") {
        combined[[col]] <- 1  # Default to 1 if not specified
      } else {
        combined[[col]] <- NA  # Add missing columns as NA
      }
    }
  }
  
  # Select only the required columns in consistent order
  combined %>% select(all_of(required_cols))
}

# Standardize all result sets
paxdb_std <- standardize_results(paxdb_results)
gpmdb_std <- standardize_results(gpmdb_results)
px_std <- standardize_results(px_results)

# Combine all results
all_results <- bind_rows(paxdb_std, gpmdb_std, px_std)

# Group related cell subtypes into main categories
group_cell_types <- function(celltype) {
  case_when(
    # NK cells (all NK subtypes)
    str_detect(celltype, "NK|nk") ~ "NK_cells",
    
    # B cells (including plasma cells)
    str_detect(celltype, "B_cells|B\\.plasma") ~ "B_cells",
    
    # Monocytes/Macrophages (all monocyte subtypes)
    str_detect(celltype, "MO\\.|Monocyte|monocyte|Macrophage") ~ "Monocytes",
    
    # CD4 T cells (including subtypes)
    str_detect(celltype, "CD4|T4\\.|Th1|Th2|Th17") ~ "CD4_T_cells",
    
    # CD8 T cells (including subtypes)  
    str_detect(celltype, "CD8|T8\\.|CD8_T_cells") ~ "CD8_T_cells",
    
    # Dendritic cells
    str_detect(celltype, "DC|pDC|mDC") ~ "Dendritic_cells",
    
    # Other immune cells
    str_detect(celltype, "Neutrophil") ~ "Neutrophils",
    str_detect(celltype, "Eosinophil") ~ "Eosinophils", 
    str_detect(celltype, "Basophil") ~ "Basophils",
    
    # Blood cells
    str_detect(celltype, "Erythrocyte|erythrocyte") ~ "Erythrocytes",
    str_detect(celltype, "Thrombocyte|platelet") ~ "Thrombocytes",
    
    # Default: keep original name
    TRUE ~ celltype
  )
}

message("  Grouping related cell subtypes...")
# Show original cell types before grouping
original_celltypes <- unique(all_results$celltype)
message(sprintf("  Original cell types (%d): %s", 
               length(original_celltypes), 
               paste(sort(original_celltypes), collapse = ", ")))

# Apply grouping
all_results <- all_results %>%
  mutate(celltype = group_cell_types(celltype))

# Show grouped cell types  
grouped_celltypes <- unique(all_results$celltype)
message(sprintf("  Grouped cell types (%d): %s", 
               length(grouped_celltypes), 
               paste(sort(grouped_celltypes), collapse = ", ")))

# Aggregate intensities for grouped cell types by taking median
all_results <- all_results %>%
  group_by(gene, celltype, source, filename) %>%
  summarise(
    intensity = median(intensity, na.rm = TRUE),
    protein_count = sum(protein_count, na.rm = TRUE),
    .groups = "drop"
  )

# Create summary statistics
celltype_summary <- all_results %>%
  group_by(celltype, source) %>%
  summarise(
    gene_count = n(),
    filename = first(filename),
    .groups = "drop"
  ) %>%
  arrange(celltype, source)

# Overall summary by cell type
overall_summary <- all_results %>%
  group_by(celltype) %>%
  summarise(
    total_genes = length(unique(gene)),
    sources = n_distinct(source),
    source_list = paste(unique(source), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(total_genes))

# 5. Generate visualizations
message("\n5. Generating visualizations...")

# Set up plot output directory
plot_dir <- "outputs/plots/05_celltype_analysis"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Function to format cell type names for display
format_celltype_names <- function(celltype) {
  str_replace_all(celltype, "_", " ")
}

# Function to format source names for display
format_source_names <- function(source) {
  case_when(
    str_detect(source, "ProteomeXchange_pxd") ~ str_extract(source, "pxd\\d+") %>% str_to_upper(),
    str_detect(source, "gpmdb_") ~ "GPMDB",
    str_detect(source, "paxdb_") ~ "PAXDB",
    TRUE ~ source
  )
}

# Create CD8 T Cell correlation analysis
message("  Creating CD8 T cell correlation analysis...")
# Filter data for CD8 T cells and prepare for correlation analysis
cd8_data <- all_results %>%
  filter(str_detect(celltype, "CD8|T8")) %>%  # Include all CD8 T cell subtypes
  group_by(gene, source) %>%
  summarise(
    intensity = median(intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Calculate z-scores within each source
  group_by(source) %>%
  mutate(
    z_score = scale(log10(intensity))[,1]
  ) %>%
  ungroup()

# Get unique sources and create empty matrices
sources <- unique(cd8_data$source)
n_sources <- length(sources)
cor_matrix <- matrix(NA, n_sources, n_sources)
n_genes_matrix <- matrix(NA, n_sources, n_sources)
rownames(cor_matrix) <- colnames(cor_matrix) <- sources
rownames(n_genes_matrix) <- colnames(n_genes_matrix) <- sources

# Fill matrices
for(i in 1:n_sources) {
  for(j in 1:n_sources) {
    source1 <- sources[i]
    source2 <- sources[j]
    
    if(i == j) {
      # Diagonal: correlation = 1, n_genes = total genes in source
      cor_matrix[i,j] <- 1
      n_genes_matrix[i,j] <- cd8_data %>%
        filter(source == source1) %>%
        pull(gene) %>%
        n_distinct()
    } else {
      # Get shared genes between sources
      shared_data <- cd8_data %>%
        filter(source %in% c(source1, source2)) %>%
        group_by(gene) %>%
        filter(n() == 2) %>%
        ungroup() %>%
        select(gene, source, z_score) %>%
        pivot_wider(names_from = source, values_from = z_score)
      
      n_shared <- nrow(shared_data)
      
      if(n_shared > 0) {
        # Calculate correlation if there are shared genes
        cor_val <- cor(shared_data[[source1]], shared_data[[source2]], 
                      method = "spearman", use = "complete.obs")
        cor_matrix[i,j] <- cor_val
        n_genes_matrix[i,j] <- n_shared
      }
    }
  }
}

# Convert to long format for plotting
plot_data <- expand.grid(
  source1 = rownames(cor_matrix),
  source2 = colnames(cor_matrix)
) %>%
  mutate(
    correlation = as.vector(cor_matrix),
    n_genes = as.vector(n_genes_matrix),
    source1 = format_source_names(source1),
    source2 = format_source_names(source2),
    label = case_when(
      source1 == source2 ~ "",  # Empty label on diagonal
      !is.na(correlation) ~ sprintf("%.2f", correlation),
      TRUE ~ "NA"
    )
  )

# Create correlation plot
p_cd8_cor <- ggplot(plot_data, aes(x = source1, y = source2, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), 
            color = ifelse(plot_data$correlation > 0.5 | plot_data$source1 == plot_data$source2, 
                          "white", "black"),
            size = 2.5) +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4",
                      midpoint = 0, limits = c(-1, 1), na.value = "gray90") +
  labs(
    title = "(D) CD8 T Cell Expression Correlation Across Sources",
    subtitle = "Correlation coefficients between sources",
    x = "Data Source",
    y = "Data Source",
    fill = "Correlation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Save CD8 correlation plot separately
ggsave(file.path(plot_dir, "04_cd8_correlation.png"), p_cd8_cor, 
       width = 8, height = 7, dpi = 300, bg = "white")

# Prepare data for intensity distributions
# Function to calculate z-score normalization by source
calculate_zscore_normalization <- function(data) {
  data %>%
    group_by(source) %>%
    mutate(
      log_intensity = log10(intensity),
      z_score = (log_intensity - mean(log_intensity, na.rm = TRUE)) / sd(log_intensity, na.rm = TRUE)
    ) %>%
    ungroup()
}

# Get top cell types for visualization
top_celltypes <- all_results %>%
  group_by(celltype) %>%
  summarise(total_genes = n_distinct(gene)) %>%
  arrange(desc(total_genes)) %>%
  slice_head(n = 8) %>%
  pull(celltype)

# Prepare intensity data with z-score normalization
intensity_sample <- all_results %>%
  filter(celltype %in% top_celltypes, !is.na(intensity), intensity > 0) %>%
  calculate_zscore_normalization() %>%
  group_by(celltype) %>%
  sample_n(size = min(1000, n()), replace = FALSE) %>%  # Sample for performance
  ungroup()

# Plot A: Gene coverage by cell type and data source
p_top_left <- all_results %>%
  group_by(celltype, source) %>%
  summarise(genes_per_source = n_distinct(gene), .groups = "drop") %>%
  mutate(
    celltype_display = format_celltype_names(celltype),
    source_display = format_source_names(source),
    source_display = case_when(
      str_detect(source_display, "^PXD") ~ paste0("", source_display),
      TRUE ~ source_display
    )
  ) %>%
  group_by(celltype_display) %>%
  mutate(total_genes = sum(genes_per_source)) %>%
  ungroup() %>%
  arrange(desc(total_genes), source_display) %>%
  mutate(
    celltype_display = factor(celltype_display, levels = unique(celltype_display))
  ) %>%
  ggplot(aes(x = celltype_display, y = genes_per_source, fill = source_display)) +
  geom_col(position = "stack", alpha = 0.8) +
  geom_text(aes(label = comma(genes_per_source)), 
            position = position_stack(vjust = 0.5),
            size = 6) +
  scale_fill_brewer(type = "qual", palette = "Set2", name = "Data Source") +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
  labs(
    title = "(A) Gene Coverage by Cell Type and Data Source",
    x = "Cell Type",
    y = "Number of Unique Genes"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.5),  # Position legend inside plot
    legend.background = element_rect(fill = "#FFFFFF99", color = NA),  # Semi-transparent white background
    legend.box.margin = margin(0, 0, 0, 0),  # Remove margin around legend box
    legend.margin = margin(5, 5, 5, 5),  # Add small margin inside legend
    plot.title = element_text(size = 28, face = "bold"),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 15.5, 5.5, 5.5)
  )

# Plot B: Intensity distributions
p_bottom_left <- intensity_sample %>%
  mutate(
    celltype_display = format_celltype_names(celltype),
    source_display = format_source_names(source)
  ) %>%
  ggplot(aes(x = source_display, y = z_score, fill = source_display)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  geom_jitter(size = 0.1, alpha = 0.1, width = 0.2) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  labs(
    title = "(B) Intensity Distributions",
    x = "Data Source",
    y = "Z-Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 28, face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 18, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 15.5, 5.5, 5.5)
  ) +
  facet_wrap(~ celltype_display, scales = "free_x", ncol = 4)

# Create correlation analysis for all cell types with multiple sources
# First identify cell types with multiple sources
multi_source_celltypes <- all_results %>%
  group_by(celltype) %>%
  summarise(n_sources = n_distinct(source)) %>%
  filter(n_sources > 1) %>%
  pull(celltype)

# Calculate correlations for each cell type with multiple sources
correlation_data_all <- all_results %>%
  filter(celltype %in% multi_source_celltypes) %>%
  group_by(celltype, gene, source) %>%
  summarise(
    intensity = median(intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Calculate z-scores within each source
  group_by(source) %>%
  mutate(
    z_score = scale(log10(intensity))[,1]
  ) %>%
  ungroup()

# Function to calculate correlations for a specific cell type
calculate_correlations <- function(data, ct) {
  ct_data <- data %>%
    filter(celltype == ct) %>%
    select(gene, source, z_score) %>%
    group_by(gene) %>%
    filter(n() > 1) %>%  # Only keep genes present in multiple sources
    ungroup()
  
  sources <- unique(ct_data$source)
  n_sources <- length(sources)
  
  if(n_sources < 2) return(NULL)
  
  cor_matrix <- matrix(NA, n_sources, n_sources)
  rownames(cor_matrix) <- colnames(cor_matrix) <- sources
  
  for(i in 1:n_sources) {
    for(j in 1:n_sources) {
      source1 <- sources[i]
      source2 <- sources[j]
      
      if(i == j) {
        cor_matrix[i,j] <- 1
      } else {
        shared_data <- ct_data %>%
          filter(source %in% c(source1, source2)) %>%
          group_by(gene) %>%
          filter(n() == 2) %>%
          ungroup() %>%
          select(gene, source, z_score) %>%
          pivot_wider(names_from = source, values_from = z_score)
        
        if(nrow(shared_data) > 0) {
          cor_val <- cor(shared_data[[source1]], shared_data[[source2]], 
                        method = "spearman", use = "complete.obs")
          cor_matrix[i,j] <- cor_val
        }
      }
    }
  }
  
  list(
    celltype = ct,
    correlations = cor_matrix
  )
}

# Calculate correlations for all relevant cell types
correlation_results <- lapply(multi_source_celltypes, function(ct) {
  calculate_correlations(correlation_data_all, ct)
})

# Function to create visualization based on number of sources
plot_correlation_viz <- function(cor_result) {
  if(is.null(cor_result)) return(NULL)
  
  ct <- cor_result$celltype
  sources <- rownames(cor_result$correlations)
  
  # Get all possible pairs of sources
  source_pairs <- combn(sources, 2, simplify = FALSE)
  
  # Create a list to store all plots for this cell type
  plots_list <- lapply(source_pairs, function(pair) {
    source1 <- pair[1]
    source2 <- pair[2]
    
    # Get the data for this pair of sources
    plot_data <- correlation_data_all %>%
      filter(celltype == ct, source %in% c(source1, source2)) %>%
      select(gene, source, z_score) %>%
      pivot_wider(names_from = source, values_from = z_score) %>%
      drop_na()
    
    # Calculate correlation
    cor_val <- cor(plot_data[[source1]], plot_data[[source2]], 
                  method = "spearman", use = "complete.obs")
    
    # Create scatterplot
    ggplot(plot_data, aes(x = .data[[source1]], y = .data[[source2]])) +
      geom_point(alpha = 0.5, size = 1.2) +
      geom_smooth(method = "lm", se = TRUE, color = "#4575b4", alpha = 0.2) +
      labs(
        title = format_celltype_names(ct),
        subtitle = sprintf("ρ = %.2f", cor_val),
        x = format_source_names(source1),
        y = format_source_names(source2)
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)
      )
  })
  
  return(plots_list)
}

# Create correlation plots for all cell types
correlation_plots <- lapply(correlation_results, plot_correlation_viz)

# Flatten the list of lists into a single list of plots
all_plots <- unlist(correlation_plots, recursive = FALSE)

# Filter out NULL plots
valid_plots <- all_plots[!sapply(all_plots, is.null)]

if (length(valid_plots) > 0) {
  # Calculate how many empty plots we need to add to make a complete 2x4 grid
  n_plots <- length(valid_plots)
  n_empty_needed <- 8 - n_plots  # We want 8 spots total (2 rows x 4 columns)
  
  # Create empty plots if needed
  if (n_empty_needed > 0) {
    empty_plot <- ggplot() + theme_void()
    empty_plots <- replicate(n_empty_needed, empty_plot, simplify = FALSE)
    all_plots_with_empty <- c(valid_plots, empty_plots)
  } else {
    all_plots_with_empty <- valid_plots[1:8]  # Take only first 8 if we have more
  }
  
  # Create title plot
  title_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "(C) Correlation across cell types and sources",
             size = 8, fontface = "bold") +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # Combine plots with title
  p_bottom_right <- wrap_plots(
    title_plot,
    wrap_plots(all_plots_with_empty, ncol = 4, nrow = 2),
    ncol = 1,
    heights = c(0.1, 1)
  )
} else {
  # Create a message plot if no valid correlations were found
  p_bottom_right <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = "No cell types with\nmultiple data sources found",
             size = 6) +
    theme_void() +
    labs(title = "(C) Correlation across cell types and sources") +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.margin = margin(20, 20, 20, 20)
    )
}

# Save individual plots in TIFF format
ggsave(file.path(plot_dir, "03_celltype_gene_counts.tiff"), p_top_left, 
       width = 12, height = 8, dpi = 300, compression = "lzw")
ggsave(file.path(plot_dir, "05_intensity_distributions_zscore.tiff"), p_bottom_left, 
       width = 12, height = 8, dpi = 300, compression = "lzw")
ggsave(file.path(plot_dir, "04_cd8_correlation.tiff"), p_bottom_right, 
       width = 8, height = 7, dpi = 300, compression = "lzw")

# Create the comprehensive plot with the layout and adjust the relative heights
library(patchwork)

# Create a layout specification for three plots
layout <- "
AAA
BCC
"

# Combine plots with the layout
p6_comprehensive <- p_top_left + p_bottom_left + p_bottom_right +
  plot_layout(
    design = layout,
    widths = c(4, 1),  # Make Panel B wider compared to Panel C
    heights = c(1, 1.2)  # Keep bottom row slightly taller
  ) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 18, hjust = 0.5),
      plot.caption = element_text(size = 14, hjust = 1)
    )
  )

# Save the comprehensive plot in TIFF format
ggsave(
  file.path(plot_dir, "00_comprehensive_summary.tiff"),
  p6_comprehensive,
  width = 24,
  height = 16,
  dpi = 300,
  compression = "lzw",
  limitsize = FALSE
)

message(sprintf("  TIFF plots saved to: %s", plot_dir))

# Save results
message("\n6. Saving results...")

# Save detailed results
write_csv(all_results, file.path(output_dir, "celltype_protein_data.csv"))

# Save summaries
write_csv(celltype_summary, file.path(output_dir, "celltype_summary_by_source.csv"))
write_csv(overall_summary, file.path(output_dir, "celltype_overall_summary.csv"))

# Print summary to console
message("\n=== Cell Type Analysis Summary ===")
print(overall_summary)

message(sprintf("\nTotal unique genes across all cell types: %d", 
               length(unique(all_results$gene))))

message(sprintf("Results saved to: %s", output_dir))
message(sprintf("TIFF plots saved to: %s", plot_dir))
message("=== Analysis Complete ===") 