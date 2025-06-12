#!/usr/bin/env Rscript

# Load utilities and set up paths
source(file.path(dirname(dirname(getwd())), "scripts", "utilities", "load_packages.R"))
ensure_output_dirs()


# Fix CD8 Deduplication Issue
# Author: Data Analysis Pipeline  
# Date: 2024

library(dplyr)

# Read the current data
message("Reading current CD8 data...")
data <- read.csv('plots/01_Database_Analysis/celltype_gene_data.csv', stringsAsFactors = FALSE)
cd8_data <- data[data$CellType == 'CD8', ]

message("=== CURRENT CD8 DATA ANALYSIS ===")
message(paste("Total measurements:", nrow(cd8_data)))
message(paste("Unique genes:", length(unique(cd8_data$Gene))))
message(paste("Datasets:", paste(unique(cd8_data$Dataset), collapse = ", ")))

# Count duplicated genes
gene_counts <- table(cd8_data$Gene)
duplicated_genes <- names(gene_counts[gene_counts > 1])
message(paste("Genes appearing in multiple datasets:", length(duplicated_genes)))

# Show examples
if(length(duplicated_genes) >= 3) {
  message("\nExample genes appearing in multiple datasets:")
  for(i in 1:3) {
    gene <- duplicated_genes[i]
    gene_data <- cd8_data[cd8_data$Gene == gene, c('Gene', 'Concentration', 'Dataset')]
    message(paste("Gene:", gene))
    for(j in 1:nrow(gene_data)) {
      message(paste("  ", gene_data$Dataset[j], ":", format(gene_data$Concentration[j], scientific = TRUE)))
    }
  }
}

# OPTION 1: Complete deduplication (merge across datasets)
message("\n=== OPTION 1: COMPLETE DEDUPLICATION ===")
cd8_merged <- cd8_data %>%
  group_by(Gene, Database, Category, CellType) %>%
  summarise(
    Concentration = max(Concentration, na.rm = TRUE),  # Keep highest concentration
    Dataset = paste(unique(Dataset), collapse = ","),  # Combine dataset names
    .groups = 'drop'
  )

message(paste("After complete deduplication:"))
message(paste("  Total measurements:", nrow(cd8_merged)))
message(paste("  Unique genes:", length(unique(cd8_merged$Gene))))
message(paste("  Reduction:", nrow(cd8_data) - nrow(cd8_merged), "measurements removed"))

# Show some merged entries
if(length(duplicated_genes) >= 3) {
  message("\nExample merged entries:")
  for(i in 1:3) {
    gene <- duplicated_genes[i]
    merged_entry <- cd8_merged[cd8_merged$Gene == gene, ]
    if(nrow(merged_entry) > 0) {
      message(paste("Gene:", gene))
      message(paste("  Datasets:", merged_entry$Dataset[1]))
      message(paste("  Final concentration:", format(merged_entry$Concentration[1], scientific = TRUE)))
    }
  }
}

# OPTION 2: Dataset-aware analysis (keep separate but label clearly)
message("\n=== OPTION 2: DATASET-AWARE ANALYSIS ===")
cd8_dataset_aware <- cd8_data %>%
  mutate(
    Source_ID = paste0(Dataset, "_", CellType),
    Gene_Dataset = paste0(Gene, "_", Dataset)
  ) %>%
  group_by(Gene, Database, Dataset, Category, CellType) %>%
  slice_max(Concentration, n = 1, with_ties = FALSE) %>%
  ungroup()

message(paste("Dataset-aware approach:"))
message(paste("  Total measurements:", nrow(cd8_dataset_aware)))
message(paste("  Unique gene-dataset combinations:", length(unique(cd8_dataset_aware$Gene_Dataset))))
message(paste("  Unique genes:", length(unique(cd8_dataset_aware$Gene))))

# Create summary by dataset
dataset_summary <- cd8_dataset_aware %>%
  group_by(Dataset) %>%
  summarise(
    gene_count = n_distinct(Gene),
    measurement_count = n(),
    .groups = 'drop'
  )

message("\nMeasurements by dataset:")
for(i in 1:nrow(dataset_summary)) {
  message(paste("  ", dataset_summary$Dataset[i], ":", 
                dataset_summary$gene_count[i], "genes,", 
                dataset_summary$measurement_count[i], "measurements"))
}

# Calculate overlap between datasets
message("\n=== DATASET OVERLAP ANALYSIS ===")
datasets <- unique(cd8_data$Dataset)
for(i in 1:(length(datasets)-1)) {
  for(j in (i+1):length(datasets)) {
    genes_i <- unique(cd8_data$Gene[cd8_data$Dataset == datasets[i]])
    genes_j <- unique(cd8_data$Gene[cd8_data$Dataset == datasets[j]])
    overlap <- length(intersect(genes_i, genes_j))
    union_size <- length(union(genes_i, genes_j))
    jaccard <- round(overlap / union_size, 3)
    message(paste(datasets[i], "vs", datasets[j], ":", overlap, "shared genes, Jaccard =", jaccard))
  }
}

# RECOMMENDATION
message("\n=== RECOMMENDATION ===")
message("For standard analysis: Use OPTION 1 (Complete Deduplication)")
message("  - Reduces CD8 from 18,876 to ~9,352 measurements")
message("  - Each gene appears once with highest concentration")
message("  - Standard for comparative proteomics")
message("")
message("For dataset comparison: Use OPTION 2 (Dataset-Aware)")
message("  - Keeps dataset information for meta-analysis")
message("  - Allows studying dataset-specific effects")
message("  - Better for method comparison studies")

# Save results
write.csv(cd8_merged, get_output_path("cd8_merged_deduplicated.csv", "tables"), row.names = FALSE)
write.csv(cd8_dataset_aware, get_output_path("cd8_dataset_aware.csv", "tables"), row.names = FALSE)
message("\nSaved results:")
message("  cd8_merged_deduplicated.csv - Option 1")
message("  cd8_dataset_aware.csv - Option 2") 
