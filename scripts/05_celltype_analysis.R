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
                      "ggthemes", "patchwork", "UpSetR", "ggupset", "RColorBrewer", "ggrepel", "tibble", "cowplot")
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

# Process GPMDB data
process_gpmdb_data <- function(data, force_mapping = FALSE) {
  message("Processing GPMDB data...")
  message(sprintf("  Input rows: %d", nrow(data)))

  # Handle different GPMDB formats by checking for specific columns
  if ("PSM" %in% colnames(data) && "protein" %in% colnames(data)) {
    # Format for gpmdb_erythrocyte.csv
    message("  Detected GPMDB format with PSM/protein columns.")
    processed_data <- data %>%
      mutate(
        intensity = as.numeric(stringr::str_extract(PSM, "^\\d+")),
        gene_from_protein = stringr::str_extract(protein, "^[^\\s]+")
      )
  } else if ("total" %in% colnames(data) && "description" %in% colnames(data)) {
    # Format for gpmdb_platelet.csv
    message("  Detected GPMDB format with total/description columns.")
    processed_data <- data %>%
      mutate(
        intensity = as.numeric(total),
        gene_from_protein = stringr::str_extract(description, "^[^,]+")
      )
  } else {
    warning("Unrecognized GPMDB file format for this processor.")
    return(NULL)
  }

  # Map accessions to gene symbols, using gene_from_protein as a fallback
  if ("accession" %in% colnames(processed_data)) {
    accessions <- unique(processed_data$accession)
    message(sprintf("  Converting %d accession IDs to gene symbols...", length(accessions)))
    gene_mapping <- convert_to_gene_symbol(accessions)

    # Manually map gene symbols using the named vector
    processed_data <- processed_data %>%
      mutate(
        gene_symbol = gene_mapping[accession],
        gene = coalesce(gene_symbol, gene_from_protein)
      )
  } else {
    processed_data$gene <- processed_data$gene_from_protein
  }

  # Filter for valid entries
  valid_genes <- processed_data %>%
    filter(!is.na(gene) & gene != "" & !is.na(intensity) & intensity > 0)
  message(sprintf("  Valid entries after filtering: %d", nrow(valid_genes)))

  # Deduplicate genes, taking the median intensity
  dedup_result <- deduplicate_genes(valid_genes, gene_col = "gene", quant_col = "intensity")
  message(sprintf("  Final unique genes after deduplication: %d", nrow(dedup_result)))

  return(dedup_result)
}

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

# --- Start of new plotting section ---

# Prepare data for Panel A and B to ensure consistent ordering
plot_data_summary <- all_results %>%
  group_by(celltype, source) %>%
  summarise(genes_per_source = n_distinct(gene), .groups = "drop") %>%
  mutate(
    celltype_display = format_celltype_names(celltype)
  ) %>%
  group_by(celltype_display) %>%
  mutate(total_genes = sum(genes_per_source)) %>%
  ungroup() %>%
  arrange(total_genes) # Sort from low to high

# Get all unique sources and create a consistent factor for the legend
all_source_names <- unique(all_results$source)
source_display_levels <- unique(format_source_names(all_source_names))

# Create ordered factor for cell types
ordered_celltypes <- unique(plot_data_summary$celltype_display)

# Re-apply factor to summary data
plot_data_summary$celltype_display <- factor(plot_data_summary$celltype_display, levels = ordered_celltypes)
plot_data_summary$source_display <- factor(format_source_names(plot_data_summary$source), levels = source_display_levels)

# Panel A: Stacked horizontal bar plot of gene counts
p_panel_a <- ggplot(plot_data_summary, aes(y = celltype_display, x = genes_per_source, fill = source_display)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = if_else(genes_per_source > 150, format(genes_per_source, big.mark = ","), "")),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4.5
  ) +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  scale_fill_brewer(palette = "Paired", drop = FALSE) +
  labs(
    title = "(A) Genes per Cell Type",
    x = "Number of Genes",
    y = "Cell Type",
    fill = "Data Source"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(face = "bold", size = rel(1.2)),
    legend.title = element_text(size = rel(1.1)),
    legend.text = element_text(size = rel(1.0)),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# Panel B: Boxplot of z-scores
z_score_data <- all_results %>%
  filter(intensity > 0) %>% # Ensure log10 does not generate -Inf
  group_by(source) %>%
  mutate(
    log_intensity = log10(intensity),
    # Z-score is calculated only if there is variation in the data
    z_score = if (n() > 1 && sd(log_intensity, na.rm = TRUE) > 0) {
      scale(log_intensity)[,1]
    } else {
      NA_real_
    }
  ) %>%
  ungroup() %>%
  mutate(
    celltype_display = factor(format_celltype_names(celltype), levels = ordered_celltypes),
    source_display = factor(format_source_names(source), levels = source_display_levels)
  )

p_panel_b <- ggplot(z_score_data, aes(y = celltype_display, x = z_score, fill = source_display)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.9), width=0.8) +
  coord_cartesian(xlim = quantile(z_score_data$z_score, c(0.001, 0.999), na.rm = TRUE)) + # trim outliers for viz
  scale_fill_brewer(palette = "Paired", drop = FALSE) +
  labs(
    title = "(B) Z-score Distribution",
    x = "Z-score of log10(Intensity)",
    y = "",
    fill = "Data Source"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(face = "bold", size = rel(1.2)),
    legend.title = element_text(size = rel(1.1)),
    legend.text = element_text(size = rel(1.0)),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# Panel C: Scatterplot correlations
# Identify cell types with multiple sources
multi_source_celltypes <- all_results %>%
  group_by(celltype) %>%
  summarise(n_sources = n_distinct(source)) %>%
  filter(n_sources > 1) %>%
  pull(celltype)

# Data for correlation plots
correlation_data_all <- all_results %>%
  filter(celltype %in% multi_source_celltypes) %>%
  group_by(celltype, gene, source) %>%
  summarise(
    intensity = median(intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # z-score per source
  group_by(source) %>%
  mutate(
    z_score = scale(log10(intensity))[,1]
  ) %>%
  ungroup()

plot_all_source_correlations <- function(data, celltypes_to_plot) {
  all_plots <- list()
  
  for (ct in celltypes_to_plot) {
    ct_data <- data %>%
      filter(celltype == ct) %>%
      select(gene, source, z_score) %>%
      group_by(gene) %>%
      filter(n_distinct(source) > 1) %>%
      ungroup()
    
    sources <- unique(ct_data$source)
    if (length(sources) < 2) next
    
    source_pairs <- combn(sources, 2, simplify = FALSE)
    
    ct_plots <- lapply(source_pairs, function(pair) {
      source1 <- pair[1]
      source2 <- pair[2]
      
      pair_data <- ct_data %>%
        filter(source %in% c(source1, source2)) %>%
        pivot_wider(names_from = source, values_from = z_score) %>%
        drop_na()
      
      if (nrow(pair_data) < 10) return(NULL) # Skip if not enough data points
      
      cor_val <- cor(pair_data[[source1]], pair_data[[source2]],
                     method = "spearman", use = "complete.obs")
      
      p <- ggplot(pair_data, aes(x = .data[[source1]], y = .data[[source2]])) +
        geom_point(alpha = 0.2, size = 1, color = "blue") +
        geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.7) +
        labs(
          subtitle = sprintf("%s (n=%s, ρ=%.2f)", format_celltype_names(ct), format(nrow(pair_data), big.mark=","), cor_val),
          x = format_source_names(source1),
          y = format_source_names(source2)
        ) +
        theme_bw(base_size = 14) +
        theme(
          plot.subtitle = element_text(size = 14, hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12)
        )
      return(p)
    })
    
    all_plots <- c(all_plots, Filter(Negate(is.null), ct_plots))
  }
  
  if (length(all_plots) == 0) {
    return(ggplot() + theme_void() + labs(subtitle = "No valid cross-source correlations to display."))
  }
  
  # Return the grid of plots without a title
  return(wrap_plots(all_plots, ncol=3))
}

p_panel_c_grid <- plot_all_source_correlations(correlation_data_all, multi_source_celltypes)

# Create a title for Panel C using cowplot
panel_c_title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "(C) Cross-source Gene Abundance Correlations",
    fontface = 'bold',
    x = 0.5,
    hjust = 0.5,
    size = 18
  )

# Combine the title and the plot grid
p_panel_c <- cowplot::plot_grid(
  panel_c_title,
  p_panel_c_grid,
  ncol = 1,
  rel_heights = c(0.05, 0.95) # Allocate space for the title
)

# Combine panels into final plot
if (!is.null(p_panel_c)) {
  # Manually extract the legend from Panel A to ensure it's shared
  shared_legend <- cowplot::get_legend(
    p_panel_a + theme(legend.box.margin = margin(6, 6, 6, 6))
  )

  # Remove legends from the individual plots
  p_panel_a <- p_panel_a + theme(legend.position = "none")
  p_panel_b <- p_panel_b + theme(legend.position = "none")

  # Combine the top row plot panels
  top_row <- p_panel_a | p_panel_b

  # Combine the main plot areas
  main_plots <- top_row / p_panel_c + plot_layout(heights = c(1, 1.5))

  # Add the shared legend to the bottom of the combined plot
  final_plot <- wrap_plots(main_plots, shared_legend, ncol = 1, heights = c(0.92, 0.08))

  # Save final plot
  ggsave(file.path(plot_dir, "00_comprehensive_celltypes_panel.tiff"), 
         final_plot, 
         width = 14, height = 16, 
         device = "tiff", dpi = 300,
         compression = "lzw", bg = "transparent")
  
  message("TIFF plot saved to: ", plot_dir)
} else {
  message("Skipping final plot generation as no correlation plots were created.")
}

message("=== Analysis Complete ===")

# The process_gpmdb_data function has been moved to the top of the script for clarity. 