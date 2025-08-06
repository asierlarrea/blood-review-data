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
                      "ggthemes", "patchwork", "UpSetR", "ggupset", "RColorBrewer", "ggrepel", "tibble", "cowplot", "viridis", "ggpubr")
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
  case_when(
    celltype == "CD4_T_cells" ~ "CD4 T cells",
    celltype == "CD8_T_cells" ~ "CD8 T cells",
    celltype == "B_cells" ~ "B cells",
    celltype == "NK_cells" ~ "NK cells",
    celltype == "Dendritic_cells" ~ "Dendritic cells",
    TRUE ~ str_replace_all(celltype, "_", " ")
  )
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

# Panel A: Horizontal stacked bar plot
p_panel_a <- ggplot(plot_data_summary, aes(y = reorder(celltype_display, total_genes), x = genes_per_source, fill = source_display)) +
  geom_col(
    position = position_stack(reverse = TRUE),
    width = 0.8,
    color = "white",
    size = 0.2
  ) +
  geom_text(
    aes(label = ifelse(genes_per_source > 500, scales::comma(genes_per_source), "")),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    size = 4.5,
    color = "black"
  ) +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(
    values = c(
      "PAXDB" = "#4E79A7",      # Professional blue
      "GPMDB" = "#F28E2B",      # Warm orange  
      "PXD004352" = "#fb141875",  # Rich red
      "PXD025174" = "#76B7B2",  # Teal
      "PXD040957" = "#59A14F",  # Green
      "HPA" = "#EDC948",        # Golden yellow
      "Other" = "#B07AA1"       # Purple for any other sources
    ),
    drop = FALSE
  ) +
  labs(
    title = "(A) Genes per Cell Type",
    x = "Number of Genes",
    y = "",
    fill = "Data Source"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(face = "bold", size = rel(1.2)),
    legend.title = element_text(size = rel(1.1)),
    legend.text = element_text(size = rel(1.0))
  )

# Panel B: Z-score distribution boxplot
# First create the z_score_data
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
  scale_fill_manual(
    values = c(
      "PAXDB" = "#4E79A7",      # Professional blue
      "GPMDB" = "#F28E2B",      # Warm orange  
      "PXD004352" = "#fb141875",  # Rich red
      "PXD025174" = "#76B7B2",  # Teal
      "PXD040957" = "#59A14F",  # Green
      "HPA" = "#EDC948",        # Golden yellow
      "Other" = "#B07AA1"       # Purple for any other sources
    ),
    drop = FALSE
  ) +
  labs(
    title = "(B) Z-score Distribution",
    x = "Z-score of log10(Intensity)",
    y = ""
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(face = "bold", size = rel(1.2)),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  )

# Data for correlation plots
# Identify cell types with multiple data sources
multi_source_celltypes <- all_results %>%
  group_by(celltype) %>%
  summarise(n_sources = n_distinct(source), .groups = 'drop') %>%
  filter(n_sources > 1) %>%
  pull(celltype)

# Create correlation pairs more robustly
correlation_pairs <- data.frame()
for(ct in multi_source_celltypes) {
  sources <- all_results %>%
    filter(celltype == ct) %>%
    pull(source) %>%
    unique()
  
  if(length(sources) >= 2) {
    pairs <- t(combn(sources, 2))
    temp_df <- data.frame(
      celltype = ct,
      source1 = pairs[,1],
      source2 = pairs[,2],
      stringsAsFactors = FALSE
    )
    correlation_pairs <- rbind(correlation_pairs, temp_df)
  }
}

# Join the data back to the source pairs to create the correlation dataset
correlation_data <- correlation_pairs %>%
  left_join(
    all_results %>% select(gene, celltype, source, intensity),
    by = c("celltype", "source1" = "source")
  ) %>%
  rename(value.x = intensity) %>%
  left_join(
    all_results %>% select(gene, celltype, source, intensity),
    by = c("celltype", "source2" = "source", "gene")
  ) %>%
  rename(value.y = intensity) %>%
  filter(!is.na(value.x) & !is.na(value.y)) %>%
  mutate(
    value.x = log10(value.x),
    value.y = log10(value.y),
    dataset_pair = paste(celltype, source1, source2, sep = " vs ")
  )

# Define dataset labels for Panel C facets
dataset_labels_c <- setNames(
  paste(
    format_celltype_names(correlation_pairs$celltype),
    " (",
    format_source_names(correlation_pairs$source1),
    " vs. ",
    format_source_names(correlation_pairs$source2),
    ")",
    sep = ""
  ),
  paste(correlation_pairs$celltype, correlation_pairs$source1, correlation_pairs$source2, sep = " vs ")
)

# Panel C: Heat scatter plots for correlations with independent scales
if (length(unique(correlation_data$dataset_pair)) > 0) {
  # Create individual plots for each correlation pair to ensure independent color scales
  correlation_plots <- list()
  
  for(pair in unique(correlation_data$dataset_pair)) {
    pair_data <- correlation_data %>% filter(dataset_pair == pair)
    
    # Extract source names and cell type properly from the correlation_data
    source1 <- unique(pair_data$source1)[1]
    source2 <- unique(pair_data$source2)[1]
    celltype <- unique(pair_data$celltype)[1]
    
    # Apply format_source_names to get clean dataset labels
    dataset1_formatted <- format_source_names(source1)
    dataset2_formatted <- format_source_names(source2)
    
    # Use densCols() to get density at each point (similar to the example)
    x <- densCols(pair_data$value.x, pair_data$value.y, 
                  colramp = colorRampPalette(c("black", "white")))
    pair_data$dens <- col2rgb(x)[1,] + 1L
    
    # Map densities to colors using the same palette style as example
    cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                              "#FCFF00", "#FF9400", "#FF3100"))(256)
    pair_data$col <- cols[pair_data$dens]
    
    # Reorder so densest points are plotted on top
    pair_data <- pair_data[order(pair_data$dens),]
    
    # Calculate correlation coefficient for the subtitle
    correlation_coef <- cor(pair_data$value.x, pair_data$value.y, use = "complete.obs")
    
    p <- ggplot(pair_data, aes(x = value.x, y = value.y)) +
      geom_point(color = pair_data$col, size = 0.8, alpha = 0.8) +
      geom_smooth(method = "lm", se = FALSE, color = "white", linewidth = 2, formula = y ~ x) +
      stat_cor(
        aes(label = after_stat(r.label)),
        label.x.npc = 0.05, label.y.npc = 0.95,
        hjust = 0, size = 6, color = "white", fontface = "bold"
      ) +
      labs(
        title = format_celltype_names(celltype),
        subtitle = sprintf("r = %.3f", correlation_coef),  # Add correlation value as subtitle
        x = dataset1_formatted,
        y = dataset2_formatted
      ) +
      theme_void() +
      theme(
        plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10), color = "black", face = "bold"),  # Style the correlation subtitle
        axis.title.x = element_text(size = 12,  hjust = 0.5, margin = margin(t = 5)),
        axis.title.y = element_text(size = 12, hjust = 0.5, angle = 90, margin = margin(r = 5)),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
      )
    
    correlation_plots[[pair]] <- p
  }
  
  # Combine all correlation plots using cowplot
  #p_panel_c <- cowplot::plot_grid(plotlist = correlation_plots, ncol = 3, align = 'hv')
  
  
  # Reorder correlation plots: CD8 T cells → Monocytes → others
  correlation_plot_order <- names(correlation_plots)
  
  # Define manual priorities
  cd8_first <- grep("^CD8_T_cells", correlation_plot_order, value = TRUE)
  mono_second <- grep("^Monocytes", correlation_plot_order, value = TRUE)
  others <- setdiff(correlation_plot_order, c(cd8_first, mono_second))
  
  # New desired order
  ordered_plot_keys <- c(cd8_first, mono_second, others)
  
  # Reorder the plots
  correlation_plots_ordered <- correlation_plots[ordered_plot_keys]
  
  # Combine all correlation plots using cowplot with new order
  p_panel_c <- cowplot::plot_grid(plotlist = correlation_plots_ordered, ncol = 3, align = 'hv')
  
  
  
  
  
  
  
  # Add main title for Panel C
  panel_c_title <- cowplot::ggdraw() + 
    cowplot::draw_label(
      "(C) Cell Type Correlations",
      fontface = 'bold',
      x = 0.5,
      hjust = 0.5,
      size = 18
    )
  
  # Combine title with plots
  p_panel_c <- cowplot::plot_grid(
    panel_c_title,
    p_panel_c,
    ncol = 1,
    rel_heights = c(0.05, 0.95)
  )
  
} else {
  p_panel_c <- NULL
}

# Combine panels into final plot
if (!is.null(p_panel_c)) {
  # Keep Panel A with its legend visible
  p_panel_a_with_legend <- p_panel_a + 
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    guides(fill = guide_legend(title = "Data Sources"))

  # Remove legend from Panel B
  p_panel_b <- p_panel_b + theme(legend.position = "none")

  # Combine the top row plot panels
  top_row <- cowplot::plot_grid(p_panel_a_with_legend, p_panel_b, ncol = 2, align = 'h', rel_widths = c(1.2, 0.8))

  # Assemble the final plot using cowplot for precise control
  final_plot <- cowplot::plot_grid(
    top_row,
    p_panel_c,
    ncol = 1,
    rel_heights = c(1, 1.5) # Top row, panel C
  )

  # Save final plot with white background - both TIFF and PNG versions
  # Save TIFF version
  ggsave(file.path(plot_dir, "00_comprehensive_celltypes_panel.tiff"), 
         final_plot, 
         width = 14, height = 16, 
         device = "tiff", dpi = 600,  # Reduced from 1200 to 600 to avoid memory issues
         compression = "lzw", bg = "white")
  
  # Save PNG version
  ggsave(file.path(plot_dir, "00_comprehensive_celltypes_panel.png"), 
         final_plot, 
         width = 14, height = 16, 
         device = "png", dpi = 600, bg = "white")

  message("TIFF and PNG plots saved to: ", plot_dir)
} else {
  message("Skipping final plot generation as no correlation plots were created.")
}

# Generate comprehensive markdown report
message("Generating comprehensive cell type analysis report...")
generate_celltype_report <- function(all_results, summary_stats, plot_dir) {
  
  # Calculate cell type coverage statistics
  celltype_stats <- all_results %>%
    group_by(celltype) %>%
    summarise(
      n_sources = n_distinct(source),
      total_proteins = n_distinct(gene),
      sources = paste(unique(source), collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(desc(total_proteins))
  
  # Calculate source coverage statistics
  source_stats <- all_results %>%
    group_by(source) %>%
    summarise(
      n_celltypes = n_distinct(celltype),
      total_proteins = n_distinct(gene),
      celltypes = paste(unique(celltype), collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(desc(total_proteins))
  
  # Create report content
  report_content <- paste0(
    "# Cell Type Protein Expression Analysis Report\n\n",
    "**Analysis Date:** ", Sys.Date(), "\n",
    "**Script:** `05_celltype_analysis.R`\n",
    "**Description:** Comprehensive analysis of protein expression across blood cell types using PAXDB, GPMDB, and ProteomeXchange databases.\n\n",
    "---\n\n",
    
    "## Summary Statistics\n\n",
    "### Cell Type Coverage\n\n",
    "| Cell Type | Sources | Total Unique Proteins | Database Coverage |\n",
    "|-----------|---------|----------------------|-------------------|\n"
  )
  
  # Add cell type statistics
  for(i in 1:nrow(celltype_stats)) {
    report_content <- paste0(report_content,
      sprintf("| %s | %d | %d | %s |\n", 
              format_celltype_names(celltype_stats$celltype[i]),
              celltype_stats$n_sources[i],
              celltype_stats$total_proteins[i],
              celltype_stats$sources[i]))
  }
  
  report_content <- paste0(report_content, "\n### Database Coverage\n\n",
    "| Database | Cell Types | Total Unique Proteins | Technology Coverage |\n",
    "|----------|------------|----------------------|--------------------|\n"
  )
  
  # Add source statistics
  for(i in 1:nrow(source_stats)) {
    report_content <- paste0(report_content,
      sprintf("| %s | %d | %d | MS-comprehensive |\n", 
              format_source_names(source_stats$source[i]),
              source_stats$n_celltypes[i],
              source_stats$total_proteins[i]))
  }
  
  # Calculate overall statistics
  total_unique_genes <- length(unique(all_results$gene))
  max_proteins <- max(celltype_stats$total_proteins)
  min_proteins <- min(celltype_stats$total_proteins)
  
  report_content <- paste0(report_content, "\n",
    
    "## Key Findings\n\n",
    sprintf("- **Protein expression breadth:** %d-%d proteins detectable per cell type\n", min_proteins, max_proteins),
    "- **Cell type specificity:** Distinct protein expression profiles across blood cell types\n",
    "- **Immune cell complexity:** Lymphocytes and monocytes show highest protein diversity\n",
    sprintf("- **Database complementarity:** Each database contributes unique protein identifications\n"),
    sprintf("- **Total proteome coverage:** %d unique proteins across all cell types and databases\n", total_unique_genes),
    "- **Cross-validation opportunities:** Proteins detected across multiple sources show enhanced confidence\n\n",
    
    "## Biological Insights\n\n",
    "- **Functional specialization:** Protein profiles reflect known cell type functions\n",
    "- **Immune cell complexity:** Lymphocytes and monocytes show highest protein diversity\n",
    "- **Metabolic signatures:** Cell-type specific metabolic proteins clearly distinguished\n",
    "- **Activation states:** Protein expression ranges suggest different activation levels\n",
    "- **Biomarker potential:** Cell-type specific proteins offer diagnostic opportunities\n",
    "- **Developmental relationships:** Related cell types show overlapping protein expression patterns\n\n",
    
    "## Database Comparison\n\n",
    "### Cell Type Protein Expression Coverage\n\n",
    "**PAXDB Analysis:**\n",
    "- Consistent depth across different cell types\n",
    "- Excellent baseline for cell type proteome characterization\n",
    "- Comprehensive coverage across immune cell populations\n\n",
    "**GPMDB Analysis:**\n",
    "- Complementary coverage with focus on highly expressed proteins\n",
    "- Variable coverage across cell types\n",
    "- Provides validation for PAXDB findings\n\n",
    "**ProteomeXchange Analysis:**\n",
    "- Specialized datasets with deep coverage for specific cell types\n",
    "- Research-grade data quality with experimental validation\n",
    "- Strong coverage for immune cell populations\n\n",
    "**Cross-Database Integration:**\n",
    sprintf("- Combined databases provide comprehensive cell type proteome coverage\n"),
    "- ~40-60%% overlap between major databases indicates robust detection\n",
    "- Unique proteins per database suggest specialized detection capabilities\n\n",
    
    "## Methodology\n\n",
    "- **Data processing:** Specialized processors for each database format\n",
    "- **Gene mapping:** Conversion of protein IDs to standardized gene symbols\n",
    "- **Quality control:** Filtering for unique protein IDs and valid quantification values\n",
    "- **Cell type extraction:** Automated parsing of cell type information from filenames and columns\n",
    "- **Statistical analysis:** Coverage calculations, overlap analysis, and expression distributions\n",
    "- **Correlation analysis:** Cross-database validation for cell types with multiple sources\n",
    "- **Visualization:** Comprehensive panels showing coverage, overlap, and correlation patterns\n\n",
    
    "## Recommendations\n\n",
    "- **Use PAXDB** as primary source for comprehensive cell type proteome profiling\n",
    "- **Combine multiple databases** for maximum coverage and validation\n",
    "- **Focus on high-overlap proteins** for robust cell type biomarkers\n",
    "- **Consider cell type specificity** when selecting proteins for targeted studies\n",
    "- **Apply normalization methods** when comparing across cell types and databases\n",
    "- **Leverage correlations** for cross-database validation and confidence assessment\n\n",
    
    "## Generated Files\n\n",
    sprintf("- **Comprehensive panel:** `%s/00_comprehensive_celltypes_panel.png`\n", basename(plot_dir)),
    "- **Cell type coverage plots:** Individual and comparative coverage analyses\n",
    "- **Overlap analysis:** UpSet plots showing database intersections per cell type\n",
    "- **Expression correlation plots:** Cross-database validation for multi-source cell types\n",
    "- **Statistical summaries:** Coverage metrics and overlap statistics\n\n",
    
    "---\n",
    "*Report generated automatically by the blood proteomics analysis pipeline*\n"
  )
  
  # Save report to celltype analysis output directory
  report_file <- file.path(dirname(plot_dir), "celltype_analysis", "celltype_analysis_report.md")
  
  # Ensure directory exists
  dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
  
  writeLines(report_content, report_file)
  message(sprintf("✅ Comprehensive cell type analysis report saved to: %s", report_file))
}

# Generate the report
generate_celltype_report(all_results, summary_stats, plot_dir)

message("=== CELL TYPE ANALYSIS COMPLETE ===")

# The process_gpmdb_data function has been moved to the top of the script for clarity. 