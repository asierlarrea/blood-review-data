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
data_file <- "cell_type.csv"
if (!file.exists(data_file)) {
  stop(paste("Error: File", data_file, "not found in current directory."))
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

# Create UpSet plot
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

# Add plot title
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
