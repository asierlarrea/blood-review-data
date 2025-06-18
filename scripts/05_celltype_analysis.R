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
source("scripts/utilities/data_loader.R")
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
  message(sprintf("  Processing %s", basename(file_path)))
  
  # Read raw data
  raw_data <- read_csv(file_path, show_col_types = FALSE)
  
  # Process using standardized function
  processed_data <- process_gpmdb_data(raw_data, force_mapping = force_mapping)
  
  # Extract cell type from filename
  celltype <- case_when(
    str_detect(basename(file_path), "b_cell") ~ "B_cells",
    str_detect(basename(file_path), "cd4") ~ "CD4_T_cells",
    str_detect(basename(file_path), "cd8") ~ "CD8_T_cells",
    str_detect(basename(file_path), "nk") ~ "NK_cells",
    str_detect(basename(file_path), "monocyte") ~ "Monocytes",
    str_detect(basename(file_path), "erythrocyte") ~ "Erythrocytes",
    str_detect(basename(file_path), "platelet") ~ "Platelets",
    TRUE ~ str_extract(basename(file_path), "(?<=gpmdb_)[^.]+")
  )
  
  # Add metadata
  processed_data$celltype <- celltype
  processed_data$source <- "GPMDB"
  processed_data$filename <- basename(file_path)
  
  gpmdb_results[[celltype]] <- processed_data
  message(sprintf("    Final genes for %s: %d", celltype, nrow(processed_data)))
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

# Group related cell subtypes
all_results <- all_results %>%
  mutate(
    celltype = case_when(
      celltype %in% c("Platelets", "Thrombocyte") ~ "Thrombocytes",
      celltype %in% c("Erythrocyte", "Erythrocytes") ~ "Erythrocytes",
      celltype %in% c("B_cells", "B.plasma") ~ "B_cells",
      celltype %in% c("CD4_T_cells", "T4.EMRA", "Th1", "Th2", "Th17") ~ "CD4_T_cells",
      celltype %in% c("CD8_T_cells", "T8.EMRA") ~ "CD8_T_cells",
      celltype %in% c("NK_cells", "NK.bright", "NK.dim") ~ "NK_cells",
      celltype %in% c("Monocytes", "MO.classical", "MO.intermediate", "MO.nonclassical",
                     "LFQ.intensity.Human_Macs_D1", "LFQ.intensity.Human_Macs_D2",
                     "LFQ.intensity.Human_Macs_D5", "LFQ.intensity.Human_Macs_D6",
                     "LFQ.intensity.Human_Macs_D7") ~ "Monocytes",
      celltype %in% c("mDC", "pDC") ~ "Dendritic_cells",
      celltype == "Basophil" ~ "Basophils",
      celltype == "Eosinophil" ~ "Eosinophils",
      celltype == "Neutrophil" ~ "Neutrophils",
      TRUE ~ celltype
    )
  )

# Print summary of cell types before and after grouping
message("  Original cell types (", n_distinct(all_results$celltype), "): ", 
        paste(sort(unique(all_results$celltype)), collapse = ", "))

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
  mutate(
    total_genes = sum(genes_per_source),
    # Only show labels for non-Erythrocytes
    show_label = !str_detect(celltype_display, "Erythrocytes")
  ) %>%
  ungroup() %>%
  arrange(desc(total_genes), source_display) %>%
  mutate(
    celltype_display = factor(celltype_display, levels = unique(celltype_display))
  ) %>%
  ggplot(aes(x = celltype_display, y = genes_per_source, fill = source_display)) +
  geom_col(position = position_stack(reverse = FALSE)) +
  geom_text(
    aes(label = if_else(show_label, format(genes_per_source, big.mark = ","), "")),
    position = position_stack(vjust = 0.5, reverse = FALSE),
    color = "black",
    size = 6
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    x = "Cell Type",
    y = "Number of Genes",
    fill = "Data Source"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Plot B: Distribution of z-scores by cell type and data source
p_bottom_left <- all_results %>%
  group_by(source) %>%
  mutate(
    z_score = scale(log10(intensity))[,1]
  ) %>%
  ungroup() %>%
  mutate(
    celltype_display = format_celltype_names(celltype),
    source_display = format_source_names(source)
  ) %>%
  ggplot(aes(x = celltype_display, y = z_score, fill = source_display)) +
  geom_violin(scale = "width", adjust = 1.5, position = position_dodge(width = 0.7)) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    x = "Cell Type",
    y = "Z-score of log10(Intensity)",
    fill = "Data Source"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Identify cell types with multiple sources
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
    select(gene, source, intensity) %>%
    # Remove any rows where intensity is NA, 0, or negative
    filter(!is.na(intensity), intensity > 0) %>%
    # Calculate z-scores within each source
    group_by(source) %>%
    mutate(
      z_score = scale(log10(intensity))[,1]
    ) %>%
    ungroup() %>%
    # Keep only genes present in multiple sources
    group_by(gene) %>%
    filter(n() > 1) %>%
    ungroup()
  
  sources <- unique(ct_data$source)
  n_sources <- length(sources)
  
  if(n_sources < 2) return(NULL)
  
  cor_matrix <- matrix(NA, n_sources, n_sources)
  n_genes_matrix <- matrix(NA, n_sources, n_sources)
  rownames(cor_matrix) <- colnames(cor_matrix) <- sources
  rownames(n_genes_matrix) <- colnames(n_genes_matrix) <- sources
  
  for(i in 1:n_sources) {
    for(j in 1:n_sources) {
      source1 <- sources[i]
      source2 <- sources[j]
      
      if(i == j) {
        cor_matrix[i,j] <- 1
        next
      }
      
      # Get data for both sources
      data1 <- ct_data %>% filter(source == source1) %>% select(gene, z_score)
      data2 <- ct_data %>% filter(source == source2) %>% select(gene, z_score)
      
      # Find shared genes with non-NA values
      shared_genes <- intersect(data1$gene, data2$gene)
      data1_shared <- data1 %>% filter(gene %in% shared_genes) %>% filter(!is.na(z_score))
      data2_shared <- data2 %>% filter(gene %in% shared_genes) %>% filter(!is.na(z_score))
      final_shared_genes <- intersect(data1_shared$gene, data2_shared$gene)
      
      # Store number of shared genes
      n_genes_matrix[i,j] <- length(final_shared_genes)
      
      # Print diagnostic information
      message(sprintf("\nDiagnostic for %s correlation between %s and %s:", ct, source1, source2))
      message(sprintf("  Total genes in %s: %d", source1, nrow(data1)))
      message(sprintf("  Total genes in %s: %d", source2, nrow(data2)))
      message(sprintf("  Shared genes: %d", length(shared_genes)))
      message(sprintf("  Valid genes in %s: %d", source1, nrow(data1_shared)))
      message(sprintf("  Valid genes in %s: %d", source2, nrow(data2_shared)))
      message(sprintf("  Final shared genes with valid values: %d", length(final_shared_genes)))
      
      if(length(final_shared_genes) < 10) {
        message(sprintf("Warning: Insufficient data for correlation between %s and %s in %s", source1, source2, ct))
        message(sprintf("Only %d genes have valid values in both datasets (minimum 10 required)", length(final_shared_genes)))
        cor_matrix[i,j] <- NA
        next
      }
      
      # Calculate correlation
      data1_final <- data1_shared %>% filter(gene %in% final_shared_genes)
      data2_final <- data2_shared %>% filter(gene %in% final_shared_genes)
      cor_matrix[i,j] <- cor(data1_final$z_score, data2_final$z_score, method = "spearman")
    }
  }
  
  return(list(correlations = cor_matrix, n_genes = n_genes_matrix))
}

# Calculate correlations for all relevant cell types
correlation_results <- lapply(multi_source_celltypes, function(ct) {
  calculate_correlations(correlation_data_all, ct)
})

# Plot correlation visualization for a cell type
plot_correlation_viz <- function(cor_result) {
  if(is.null(cor_result)) return(NULL)
  
  # Extract correlation matrix and number of genes
  cor_matrix <- cor_result$correlations
  n_genes <- cor_result$n_genes
  
  # Check if we have any valid correlations
  if(all(is.na(cor_matrix))) {
    message("No valid correlations found for this cell type")
    return(NULL)
  }
  
  # Convert matrices to long format for plotting
  cor_data <- cor_matrix %>%
    as.data.frame() %>%
    rownames_to_column("source1") %>%
    pivot_longer(-source1, names_to = "source2", values_to = "correlation")
  
  n_genes_data <- n_genes %>%
    as.data.frame() %>%
    rownames_to_column("source1") %>%
    pivot_longer(-source1, names_to = "source2", values_to = "n_genes")
  
  # Combine correlation and n_genes data
  plot_data <- cor_data %>%
    left_join(n_genes_data, by = c("source1", "source2")) %>%
    mutate(
      source1 = format_source_names(source1),
      source2 = format_source_names(source2),
      label = if_else(!is.na(correlation),
                     sprintf("r = %.2f\n(n = %d)", correlation, n_genes),
                     sprintf("No correlation\n(n = %d)", n_genes))
    )
  
  # Create correlation plot
  p <- ggplot(plot_data, aes(x = source1, y = source2, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = label), size = 3) +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0, limit = c(-1,1), na.value = "grey90"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  return(p)
}

# Function to create CD8 T cell correlation plot
plot_cd8_correlations <- function(data) {
  # Filter for CD8 T cells only
  cd8_data <- data %>%
    filter(celltype == "CD8_T_cells") %>%
    select(gene, source, z_score) %>%
    group_by(gene) %>%
    filter(n() > 1) %>%
    ungroup()
  
  # Get unique sources
  sources <- unique(cd8_data$source)
  
  # Create all possible pairs of sources
  source_pairs <- combn(sources, 2, simplify = FALSE)
  
  # Create plots for each pair
  plots <- lapply(source_pairs, function(pair) {
    source1 <- pair[1]
    source2 <- pair[2]
    
    # Get data for this pair
    pair_data <- cd8_data %>%
      filter(source %in% c(source1, source2)) %>%
      pivot_wider(names_from = source, values_from = z_score) %>%
      drop_na()
    
    # Calculate correlation
    cor_val <- cor(pair_data[[source1]], pair_data[[source2]], 
                   method = "spearman", use = "complete.obs")
    
    # Create scatter plot
    ggplot(pair_data, aes(x = .data[[source1]], y = .data[[source2]])) +
      geom_point(alpha = 0.5, size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "#4575b4", alpha = 0.2) +
      labs(
        title = sprintf("CD8 T cells (n = %d)", nrow(pair_data)),
        subtitle = sprintf("ρ = %.2f", cor_val),
        x = format_source_names(source1),
        y = format_source_names(source2)
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)
      )
  })
  
  # Combine plots
  wrap_plots(plots, ncol = 2)
}

# Create correlation plots
correlation_plots <- lapply(correlation_results, plot_correlation_viz)

# Remove NULL plots and combine valid plots
valid_plots <- correlation_plots[!sapply(correlation_plots, is.null)]

if(length(valid_plots) > 0) {
  # Create combined plot with valid plots only
  combined_plot <- wrap_plots(valid_plots, ncol = 2)
  
  # # # Also save as TIFF for manuscript
  # ggsave(file.path(plot_dir, "correlation_plots.tiff"), 
  #        combined_plot, 
  #        width = 12, height = 8, 
  #        device = "tiff", dpi = 300)
} else {
  message("No valid correlation plots were generated.")
}

# Create CD8 T cell correlation plot
p_cd8_cor <- plot_cd8_correlations(correlation_data_all)

# # Save individual plots
# ggsave(file.path(plot_dir, "01_gene_coverage.tiff"), p_top_left, 
#        width = 12, height = 8, dpi = 300, compression = "lzw")

# Define layout
layout <- "
AABB
CCDD
"

# Combine plots with the layout
p6_comprehensive <- p_top_left + p_bottom_left + p_cd8_cor +
  plot_layout(
    design = layout,
    guides = 'collect'
  ) &
  theme(legend.position = 'bottom')

# Also save as TIFF for manuscript
ggsave(file.path(plot_dir, "00_comprehensive_celltypes_analysis_panel.tiff"), 
       p6_comprehensive, 
       width = 16, height = 12, 
       device = "tiff", dpi = 300)

message("TIFF plots saved to: ", plot_dir)
message("=== Analysis Complete ===")

# Process GPMDB data
process_gpmdb_data <- function(data, force_mapping = FALSE) {
  message("Processing GPMDB data...")
  message(sprintf("  Input rows: %d", nrow(data)))
  
  # Map transcript accessions to gene symbols
  message("  Mapping transcript accessions to genes...")
  data <- data %>%
    mutate(
      # Extract gene names from protein column if available
      gene_from_protein = str_extract(protein, "^[^\\s]+"),
      # Clean up the PSM values and convert to numeric
      psm = as.numeric(str_replace_all(PSM, "[^0-9.]", "")),
      # Use PSM count as intensity
      intensity = psm
    )
  
  # Map accessions to gene symbols
  if(force_mapping || !all(str_detect(data$gene_from_protein, "^[A-Z0-9]+$"))) {
    accessions <- unique(data$accession)
    message(sprintf("Converting %d IDs to gene symbols...", length(accessions)))
    gene_mapping <- convert_to_gene_symbol(accessions)
    
    data <- data %>%
      left_join(gene_mapping, by = c("accession" = "protein_id")) %>%
      mutate(
        gene = coalesce(gene_symbol, gene_from_protein)
      )
  } else {
    data <- data %>%
      mutate(gene = gene_from_protein)
  }
  
  message("  Extracting gene names from descriptions...")
  valid_genes <- data %>%
    filter(!is.na(gene)) %>%
    filter(gene != "") %>%
    filter(!is.na(intensity))
  
  message(sprintf("  Valid genes after mapping: %d", nrow(valid_genes)))
  
  # Deduplicate genes by taking median intensity
  dedup_result <- deduplicate_genes(valid_genes, intensity_col = "intensity")
  message(sprintf("  Final unique genes after deduplication: %d", nrow(dedup_result)))
  
  return(dedup_result)
} 