#!/usr/bin/env Rscript

# Script to add gene names to biomarker file using simple_id_mapping.R functions

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Source the simple_id_mapping.R script
source("scripts/data_processing/simple_id_mapping.R")

# Read the biomarker file
biomarker_file <- "data/metadata/markerdb_biomarker.csv"
biomarkers <- read_csv(biomarker_file, col_names = TRUE, show_col_types = FALSE)

# Get the protein IDs (they are in the first column)
protein_ids <- biomarkers[[1]][-1]  # Skip the header row

# Convert protein IDs to gene symbols
gene_symbols <- convert_to_gene_symbol(protein_ids)

# Create a new dataframe with the gene names
result <- data.frame(
  protein_id = protein_ids,
  gene_name = gene_symbols
)

# Save the result
output_file <- "data/metadata/markerdb_biomarker_with_genes.csv"
write_csv(result, output_file)

message(sprintf("Added gene names to biomarker file. Results saved to: %s", output_file)) 