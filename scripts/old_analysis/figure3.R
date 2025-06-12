# ============================================================================

# Load utilities and set up paths
source("scripts/utilities/load_packages.R")
ensure_output_dirs()

# ============================================================================
# Figure 3: Cell Type Proteome Bubble Plot
# Description: Creates a bubble plot showing protein counts by database and cell type
# ============================================================================

  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

# Load required libraries
if (!require(ggplot2, quietly = TRUE)) {
  library(ggplot2)
}
if (!require(scales, quietly = TRUE)) {
  library(scales)
}

# Load data with error handling
data_file <- get_data_path("bubble.csv")
if (!file.exists(data_file)) {
  stop(paste("Error: File bubble.csv not found in data/raw/ directory."))
}

proteomics_data <- read.csv(data_file, stringsAsFactors = FALSE)

# Validate data
if (nrow(proteomics_data) == 0) {
  stop("Error: Data file is empty.")
}

# Check required columns
required_cols <- c("Database", "Cell.type", "Protein.count")
missing_cols <- setdiff(required_cols, names(proteomics_data))
if (length(missing_cols) > 0) {
  stop(paste("Error: Missing columns in data:", 
             paste(missing_cols, collapse = ", ")))
}

# Define cell type order (from most to least specialized)
cell_type_order <- rev(c("CD8", "CD4", "B cell", "NK", 
                        "Platelet", "Erythrocyte", 
                        "Monocyte", "DC", "Macrophage", 
                        "Neutrophil", "Eosinophil", "Basophil"))

# Create the bubble plot
p <- ggplot(proteomics_data, aes(x = Database, 
                                y = factor(`Cell.type`, 
                                          levels = cell_type_order),
                                size = `Protein.count`, 
                                color = Database)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(3, 15), 
                       name = "Protein Count",
                       labels = scales::comma) +
  scale_color_discrete(name = "Database", 
                       labels = function(x) gsub("_", " (", gsub("_([^_]+)$", " (\\1))", x))) +
  theme_minimal() +
  labs(x = "Database", 
       y = "Cell Type", 
       title = "Protein Counts by Database and Cell Type",
       subtitle = "Bubble size represents number of proteins") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "right") +
  guides(color = "none")  # Hide color legend since it's redundant with x-axis

# Create output directory if it doesn't exist
output_dir <- "outputs/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# === Save as TIFF ===
tiff_file <- file.path(output_dir, "cell_types_plot.tiff")
cat("Saving plot to:", tiff_file, "\n")

tiff(tiff_file, width = 10, height = 8, units = "in", res = 600, compression = "lzw")

# Generate the plot for saving
print(p)

# Close TIFF device
dev.off()
cat("Plot saved successfully!\n")

# === Display plot in R console/RStudio ===
print(p)

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total observations:", nrow(proteomics_data), "\n")
cat("Unique databases:", length(unique(proteomics_data$Database)), "\n")
cat("Unique cell types:", length(unique(proteomics_data$Cell.type)), "\n")
cat("Protein count range:", min(proteomics_data$Protein.count, na.rm = TRUE), 
    "to", max(proteomics_data$Protein.count, na.rm = TRUE), "\n")

# Show top 5 highest protein counts
top_5 <- proteomics_data[order(proteomics_data$Protein.count, 
                               decreasing = TRUE), ][1:5, ]
cat("\nTop 5 highest protein counts:\n")
print(top_5[, c("Database", "Cell.type", "Protein.count")])

