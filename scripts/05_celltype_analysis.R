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
ensure_output_dirs()

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", 
                      "ggthemes", "patchwork", "UpSetR", "ggupset", "RColorBrewer", "ggrepel")
load_packages(required_packages)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
force_mapping <- "--force-mapping" %in% args

# Load gene mapping and deduplication utilities
source("scripts/data_processing/simple_id_mapping.R")
source("scripts/utilities/gene_deduplication.R")

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
    TRUE ~ source
  )
}

# Define color palette for cell types
n_celltypes <- length(unique(all_results$celltype))
celltype_colors <- RColorBrewer::brewer.pal(min(n_celltypes, 12), "Set3")
if (n_celltypes > 12) {
  celltype_colors <- colorRampPalette(celltype_colors)(n_celltypes)
}

# Plot 1: Cell type gene counts summary
p1 <- overall_summary %>%
  mutate(
    celltype_display = format_celltype_names(celltype),
    celltype_display = reorder(celltype_display, total_genes)
  ) %>%
  ggplot(aes(x = celltype_display, y = total_genes, fill = celltype)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = comma(total_genes)), hjust = -0.1, size = 3) +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
  labs(
    title = "Gene Coverage by Cell Type",
    subtitle = "Total unique genes detected per cell type across all data sources",
    x = "Cell Type",
    y = "Number of Unique Genes",
    caption = "Source: PAXDB (cell types only), GPMDB (cell types only), ProteomeXchange (PXD004352, PXD025174, PXD040957)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_dir, "celltype_gene_counts.png"), p1, 
       width = 12, height = 8, dpi = 300, bg = "white")

# Plot 2: Data source coverage
# First format source names, then aggregate by formatted names
source_summary <- all_results %>%
  mutate(source_display = format_source_names(source)) %>%
  group_by(source_display) %>%
  summarise(
    total_genes = length(unique(gene)),
    cell_types = n_distinct(celltype),
    .groups = "drop"
  ) %>%
  arrange(desc(total_genes))

p2 <- source_summary %>%
  mutate(source_display = reorder(source_display, total_genes)) %>%
  ggplot(aes(x = source_display, y = total_genes)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = paste0(comma(total_genes), "\ngenes\n(", cell_types, " types)")), 
            hjust = -0.1, size = 3, lineheight = 0.8) +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title = "Data Source Coverage",
    subtitle = "Total genes and cell types per data source",
    x = "Data Source",
    y = "Number of Unique Genes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_dir, "data_source_coverage.png"), p2, 
       width = 10, height = 6, dpi = 300, bg = "white")

# Plot 3: Cell type by source heatmap
# First format source names, then aggregate data that might have multiple sources mapping to same display name
celltype_source_formatted <- celltype_summary %>%
  select(celltype, source, gene_count) %>%
  mutate(source_display = format_source_names(source)) %>%
  group_by(celltype, source_display) %>%
  summarise(gene_count = sum(gene_count), .groups = "drop")

# Create matrix for heatmap
celltype_source_matrix <- celltype_source_formatted %>%
  pivot_wider(names_from = source_display, values_from = gene_count, values_fill = 0)

# Convert to long format for ggplot
heatmap_data <- celltype_source_matrix %>%
  pivot_longer(cols = -celltype, names_to = "source_display", values_to = "gene_count") %>%
  mutate(
    has_data = gene_count > 0,
    gene_count_log = ifelse(gene_count > 0, log10(gene_count), 0),
    celltype_display = format_celltype_names(celltype)
  )

p3 <- heatmap_data %>%
  ggplot(aes(x = source_display, y = celltype_display, fill = gene_count_log)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(gene_count > 0, comma(gene_count), "")), 
            size = 2.5, color = "black") +
  scale_fill_gradient(low = "#f7f7f7", high = "#2171b5", 
                     name = "Log10\nGenes", 
                     labels = function(x) ifelse(x == 0, "0", paste0("10^", round(x, 1)))) +
  labs(
    title = "Cell Type × Data Source Matrix",
    subtitle = "Gene counts per cell type and data source",
    x = "Data Source",
    y = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank()
  )

ggsave(file.path(plot_dir, "celltype_source_matrix.png"), p3, 
       width = 12, height = 10, dpi = 300, bg = "white")

# Plot 4: Intensity distributions by cell type (sample of top cell types)
top_celltypes <- overall_summary %>% 
  slice_head(n = 8) %>% 
  pull(celltype)

intensity_sample <- all_results %>%
  filter(celltype %in% top_celltypes, !is.na(intensity), intensity > 0) %>%
  group_by(celltype) %>%
  sample_n(size = min(1000, n()), replace = FALSE) %>%  # Sample for performance
  ungroup()

p4 <- intensity_sample %>%
  mutate(celltype_display = format_celltype_names(celltype)) %>%
  ggplot(aes(x = log10(intensity), fill = celltype)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~celltype_display, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = celltype_colors[1:length(top_celltypes)]) +
  labs(
    title = "Protein Intensity Distributions",
    subtitle = "Log10-transformed intensity values for top 8 cell types",
    x = "Log10(Intensity)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave(file.path(plot_dir, "intensity_distributions.png"), p4, 
       width = 12, height = 10, dpi = 300, bg = "white")

# Plot 5: Technology comparison (PAXDB vs ProteomeXchange)
tech_comparison <- all_results %>%
  mutate(
    technology = case_when(
      str_detect(source, "PAXDB") ~ "PAXDB",
      str_detect(source, "ProteomeXchange") ~ "ProteomeXchange",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(technology, celltype) %>%
  summarise(gene_count = n_distinct(gene), .groups = "drop") %>%
  pivot_wider(names_from = technology, values_from = gene_count, values_fill = 0) %>%
  filter(PAXDB > 0 | ProteomeXchange > 0)

p5 <- tech_comparison %>%
  mutate(celltype_display = format_celltype_names(celltype)) %>%
  ggplot(aes(x = PAXDB, y = ProteomeXchange)) +
  geom_point(size = 3, alpha = 0.7, color = "steelblue") +
  geom_text_repel(aes(label = celltype_display), size = 3, max.overlaps = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(labels = comma_format()) +
  scale_y_continuous(labels = comma_format()) +
  labs(
    title = "PAXDB vs ProteomeXchange Coverage",
    subtitle = "Gene counts per cell type by technology platform",
    x = "PAXDB Gene Count",
    y = "ProteomeXchange Gene Count",
    caption = "Diagonal line represents equal coverage"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggsave(file.path(plot_dir, "technology_comparison.png"), p5, 
       width = 10, height = 8, dpi = 300, bg = "white")

# Plot 6: Comprehensive summary dashboard
p6_top <- (p1 | p4) 
p6_bottom <- p3
p6_comprehensive <- p6_top / p6_bottom + 
  plot_annotation(
    title = "Cell Type Protein Expression Analysis - Comprehensive Summary",
    subtitle = paste0("Analysis of ", length(unique(all_results$gene)), 
                     " unique genes across ", length(unique(all_results$celltype)), 
                     " cell types from ", length(unique(all_results$source)), " data sources"),
    caption = "Generated by: 05_celltype_analysis.R"
  ) &
  theme(plot.title = element_text(size = 16, face = "bold"))

ggsave(file.path(plot_dir, "comprehensive_summary.png"), p6_comprehensive, 
       width = 20, height = 14, dpi = 300, bg = "white")

message(sprintf("  Plots saved to: %s", plot_dir))

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
message(sprintf("Plots saved to: %s", plot_dir))
message("=== Analysis Complete ===") 