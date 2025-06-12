# ============================================================================

# Load utilities and set up paths
source("scripts/utilities/load_packages.R")
source("scripts/utilities/display_formatting.R")
ensure_output_dirs()

# ============================================================================
# Figure 2: UpSet Plot of Plasma Proteins
# Description: Creates an UpSet plot showing protein intersections across databases
# ============================================================================

# Load required packages
required_packages <- c("UpSetR", "ggplot2", "grid")
load_packages(required_packages)

# Load data with error handling
data_file <- "data/raw/plasma_upset.csv"
if (!file.exists(data_file)) {
  stop(paste("Error: File", data_file, "not found in current directory."))
}

plot_df <- read.csv(data_file, stringsAsFactors = FALSE)
plot_df <- as.data.frame(plot_df)

# Validate data
if (nrow(plot_df) == 0) {
  stop("Error: Data file is empty.")
}

# Define variables for plot aesthetics (use original data column names)
set_vars <- c("GPMDB", "PaxDB", "PeptideAtlas", 
              "HPA_Immuno", "HPA_MS", "HPA_PEA")
text_scale <- c(1.3, 1.3, 1, 1, 1.3, 1.2)
main_bar_col <- "violetred4"
sets_bar_col <- "turquoise4"
matrix_col <- "slateblue4"
shade_col <- "wheat4"

# Verify that all required columns exist
missing_cols <- setdiff(set_vars, names(plot_df))
if (length(missing_cols) > 0) {
  stop(paste("Error: Missing columns in data:", 
             paste(missing_cols, collapse = ", ")))
}

# Create output directory if it doesn't exist
output_dir <- "outputs/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# === Save as TIFF ===
tiff_file <- file.path(output_dir, "upset_plot.tiff")
cat("Saving plot to:", tiff_file, "\n")

tiff(tiff_file, width = 12, height = 8, units = "in", res = 600, compression = "lzw")

# Create a copy with formatted column names for display
plot_df_display <- plot_df
names(plot_df_display) <- format_database_names(names(plot_df_display))
formatted_set_vars <- format_database_names(set_vars)

# Generate the plot for saving
upset(
  plot_df_display, 
  sets = formatted_set_vars,
  mainbar.y.label = "Identified proteins",
  sets.x.label = "Proteins by database",
  order.by = "freq", 
  point.size = 2,
  line.size = 1,
  text.scale = text_scale,
  main.bar.color = main_bar_col,
  sets.bar.color = sets_bar_col,
  matrix.color = matrix_col,
  shade.color = shade_col
)

# Add title to saved plot
grid.text("UpSet Plot of Plasma Proteins", 
          x = 0.7, y = 0.95, 
          gp = gpar(fontsize = 12, fontface = "bold"))

# Close TIFF device
dev.off()
cat("Plot saved successfully!\n")

# === Display plot in R console/RStudio ===
upset(
  plot_df_display, 
  sets = formatted_set_vars,
  mainbar.y.label = "Identified proteins",
  sets.x.label = "Proteins by database",
  order.by = "freq", 
  point.size = 2,
  line.size = 1,
  text.scale = text_scale,
  main.bar.color = main_bar_col,
  sets.bar.color = sets_bar_col,
  matrix.color = matrix_col,
  shade.color = shade_col
)

# Add title to displayed plot
grid.text("UpSet Plot of Plasma Proteins", 
          x = 0.7, y = 0.95, 
          gp = gpar(fontsize = 12, fontface = "bold"))

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total rows in dataset:", nrow(plot_df), "\n")
for (i in seq_along(set_vars)) {
  col <- set_vars[i]
  formatted_col <- formatted_set_vars[i]
  if (col %in% names(plot_df)) {
    count <- sum(plot_df[[col]], na.rm = TRUE)
    cat(paste(formatted_col, ":", count, "proteins\n"))
  }
}
