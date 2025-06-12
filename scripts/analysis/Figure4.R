# ============================================================================

# Load utilities and set up paths
source("scripts/utilities/load_packages.R")
ensure_output_dirs()

# ============================================================================
# Figure 4: Cell Type Proteome UpSet Plot
# Description: Creates an UpSet plot comparing protein overlaps across cell types
# ============================================================================

# Set CRAN mirror for package installation
if (length(getOption("repos")) == 0 || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

# Load required libraries
if (!require(UpSetR, quietly = TRUE)) {
  cat("Installing UpSetR package...\n")
  install.packages("UpSetR", dependencies = TRUE)
  library(UpSetR)
}
if (!require(ggplot2, quietly = TRUE)) {
  cat("Installing ggplot2 package...\n")
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require(grid, quietly = TRUE)) {
  cat("Installing grid package...\n")
  install.packages("grid", dependencies = TRUE)
  library(grid)
}

# Load data with error handling
data_file <- get_data_path("upset_cell_type.csv")
if (!file.exists(data_file)) {
  stop(paste("Error: File upset_cell_type.csv not found in data/raw/ directory."))
}

plot_df <- read.csv(data_file, stringsAsFactors = FALSE)

# Validate data
if (nrow(plot_df) == 0) {
  stop("Error: Data file is empty.")
}

# Define sets (cell types) to analyze
set_vars <- c("CD8", "CD4", "Platelet", "NK", "Monocyte", "Erythrocyte", 
              "Bcell", "DC", "Macrophage", "Basophil", "Eosinophil", 
              "Neutrophil")

# Verify that all required columns exist
missing_cols <- setdiff(set_vars, names(plot_df))
if (length(missing_cols) > 0) {
  warning(paste("Missing columns in data:", 
                paste(missing_cols, collapse = ", ")))
  # Use only available columns
  set_vars <- intersect(set_vars, names(plot_df))
}

if (length(set_vars) == 0) {
  stop("Error: No valid cell type columns found in data.")
}

# Define aesthetic settings
text_scale <- c(1.3, 1.3, 1, 1, 1.3, 1.2)
main_bar_col <- "violetred4"
sets_bar_col <- "turquoise4"
matrix_col <- "slateblue4"
shade_col <- "wheat4"

# Create output directory if it doesn't exist
output_dir <- "outputs/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# === Save as TIFF ===
tiff_file <- file.path(output_dir, "upset_cell_types_plot.tiff")
cat("Saving plot to:", tiff_file, "\n")

tiff(tiff_file, width = 12, height = 8, units = "in", res = 600, compression = "lzw")

# Generate the UpSet plot for saving
upset(
  plot_df, 
  sets = set_vars,
  mainbar.y.label = "Proteins",
  sets.x.label = "Proteins by cell type",
  order.by = "freq", 
  point.size = 2,
  line.size = 1,
  text.scale = text_scale,
  main.bar.color = main_bar_col,
  sets.bar.color = sets_bar_col,
  matrix.color = matrix_col,
  shade.color = shade_col
)

# Add plot title to saved version
grid.text("Cell Type Proteome Comparison", 
          x = 0.7, y = 0.95, 
          gp = gpar(fontsize = 12, fontface = "bold"))

# Close TIFF device
dev.off()
cat("Plot saved successfully!\n")

# === Display plot in R console/RStudio ===
upset(
  plot_df, 
  sets = set_vars,
  mainbar.y.label = "Proteins",
  sets.x.label = "Proteins by cell type",
  order.by = "freq", 
  point.size = 2,
  line.size = 1,
  text.scale = text_scale,
  main.bar.color = main_bar_col,
  sets.bar.color = sets_bar_col,
  matrix.color = matrix_col,
  shade.color = shade_col
)

# Add plot title to displayed version
grid.text("Cell Type Proteome Comparison", 
          x = 0.7, y = 0.95, 
          gp = gpar(fontsize = 12, fontface = "bold"))

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total rows in dataset:", nrow(plot_df), "\n")
cat("Cell types analyzed:", length(set_vars), "\n")
cat("Available cell types:", paste(set_vars, collapse = ", "), "\n")

for (col in set_vars) {
  if (col %in% names(plot_df)) {
    count <- sum(plot_df[[col]], na.rm = TRUE)
    cat(paste(col, ":", count, "proteins\n"))
  }
}
