# Display formatting utilities for consistent plot labels
# Converts underscore formats to parentheses for publication-ready plots

# Format database names for display
format_database_names <- function(db_names) {
  formatted <- db_names
  formatted <- gsub("HPA_MS", "HPA (MS)", formatted)
  formatted <- gsub("HPA_Immunoassay", "HPA (Immunoassay)", formatted)
  formatted <- gsub("HPA_PEA", "HPA (PEA)", formatted)
  formatted <- gsub("HPA_Immuno", "HPA (Immuno)", formatted)
  formatted <- gsub("GPMDB_Plasma", "GPMDB (Plasma)", formatted)
  formatted <- gsub("GPMDB_Serum", "GPMDB (Serum)", formatted)
  formatted <- gsub("GPMDB_Erythrocyte", "GPMDB (Erythrocyte)", formatted)
  formatted <- gsub("GPMDB_Platelet", "GPMDB (Platelet)", formatted)
  formatted <- gsub("PaxDB_Plasma", "PaxDB (Plasma)", formatted)
  formatted <- gsub("PeptideAtlas_Plasma", "PeptideAtlas (Plasma)", formatted)
  
  return(formatted)
}

# Format cell type names for display  
format_cell_type_names <- function(cell_names) {
  formatted <- cell_names
  formatted <- gsub("B_cell", "B cell", formatted)
  formatted <- gsub("_", " ", formatted)
  
  return(formatted)
}

# Format technical terms for display
format_technical_terms <- function(terms) {
  formatted <- terms
  formatted <- gsub("log_abundance", "log abundance", formatted)
  formatted <- gsub("gene_symbol", "gene symbol", formatted)
  formatted <- gsub("abundance_type", "abundance type", formatted)
  formatted <- gsub("dynamic_range", "dynamic range", formatted)
  formatted <- gsub("biomarker_coverage", "biomarker coverage", formatted)
  formatted <- gsub("protein_count", "protein count", formatted)
  formatted <- gsub("cell_type", "cell type", formatted)
  
  return(formatted)
}

# Custom scale functions with proper labeling
scale_color_database <- function(name = "Database", ...) {
  scale_color_discrete(name = name, labels = format_database_names, ...)
}

scale_fill_database <- function(name = "Database", ...) {
  scale_fill_discrete(name = name, labels = format_database_names, ...)
}

# Function to add formatted labels to existing ggplot
add_formatted_database_labels <- function(plot_obj) {
  # Add or modify scale_*_discrete to include proper formatting
  plot_obj + 
    scale_color_discrete(labels = format_database_names) +
    scale_fill_discrete(labels = format_database_names)
}

# Create database color palette with formatted names
create_database_colors <- function() {
  colors <- c(
    "HPA (MS)" = "#1F78B4",
    "HPA (Immunoassay)" = "#33A02C", 
    "HPA (PEA)" = "#FF7F00",
    "HPA (Immuno)" = "#33A02C",
    "GPMDB (Plasma)" = "#E31A1C",
    "GPMDB (Serum)" = "#FB9A99",
    "GPMDB (Erythrocyte)" = "#FDBF6F",
    "GPMDB (Platelet)" = "#FF7F00",
    "PaxDB (Plasma)" = "#6A3D9A",
    "PeptideAtlas (Plasma)" = "#CAB2D6"
  )
  return(colors)
}

# Print formatting examples
show_formatting_examples <- function() {
  cat("Database Name Formatting Examples:\n")
  cat("==================================\n")
  test_names <- c("HPA_MS", "HPA_Immunoassay", "HPA_PEA", "GPMDB_Plasma", 
                 "GPMDB_Serum", "PaxDB_Plasma", "PeptideAtlas_Plasma")
  formatted_names <- format_database_names(test_names)
  
  for (i in 1:length(test_names)) {
    cat(sprintf("%-20s -> %s\n", test_names[i], formatted_names[i]))
  }
  
  cat("\nCell Type Formatting Examples:\n")
  cat("==============================\n")
  cell_types <- c("B_cell", "NK_cell", "T_cell", "cell_type")
  formatted_cells <- format_cell_type_names(cell_types)
  
  for (i in 1:length(cell_types)) {
    cat(sprintf("%-15s -> %s\n", cell_types[i], formatted_cells[i]))
  }
} 