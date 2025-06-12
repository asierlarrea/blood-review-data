# ============================================================================

# Load utilities and set up paths
source("scripts/utilities/load_packages.R")
ensure_output_dirs()

# ============================================================================
# Figure 5: Cell Type Intensity Distribution Boxplot
# Description: Creates boxplots showing protein intensity distributions by cell type
# ============================================================================

  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

# Load required libraries
required_packages <- c("ggplot2", "tidyr", "dplyr", 
                       "readr", "scales")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
  }
}

# Load data with error handling
data_file <- "data/raw/boxplot.csv"
if (!file.exists(data_file)) {
  stop(paste("Error: File", data_file, "not found in current directory."))
}

cell_data <- read_csv(data_file)

# Validate data
if (nrow(cell_data) == 0) {
  stop("Error: Data file is empty.")
}

# Convert intensity columns to numeric
intensity_cols <- grep("^Intensity_", names(cell_data), 
                       value = TRUE)
if (length(intensity_cols) == 0) {
  stop("Error: No intensity columns found (should start with 'Intensity_')")
}

cell_data <- cell_data %>% 
  mutate(across(all_of(intensity_cols), ~parse_number(as.character(.))))

# Define cell type order (hierarchical: adaptive -> innate -> blood cells)
cell_type_levels <- c("CD8", "CD4", "B cell", "NK", 
                     "Platelet", "Erythrocyte", 
                     "Monocyte", "DC", "Macrophage", 
                     "Neutrophil", "Eosinophil", "Basophil")

# Pivot to long format
long_data <- cell_data %>% 
  pivot_longer(cols = all_of(intensity_cols), 
               names_to = "cell_type", 
               values_to = "intensity") %>%
  mutate(
    cell_type = gsub("Intensity_", "", cell_type),
    # Format cell type names properly
    cell_type = gsub("_", " ", cell_type),
    cell_type = factor(cell_type, levels = cell_type_levels)
  ) %>%
  filter(!is.na(intensity))  # Remove NA values

# Calculate summary statistics
summary_stats <- long_data %>%
  group_by(cell_type) %>%
  summarise(
    protein_count = n(),
    mean_intensity = mean(intensity, na.rm = TRUE),
    median_intensity = median(intensity, na.rm = TRUE),
    .groups = 'drop'
  )

# Count number of databases per cell type (based on available data)
database_count <- long_data %>%
  group_by(cell_type) %>%
  summarise(
    unique_proteins = n_distinct(intensity[!is.na(intensity)]),
    .groups = 'drop'
  ) %>%
  mutate(
    db_label = case_when(
      cell_type %in% c("CD8") ~ "Databases: 4",
      cell_type %in% c("CD4") ~ "Databases: 3", 
      cell_type %in% c("B cell", "NK", "Platelet", 
                       "Erythrocyte", "Monocyte") ~ "Databases: 2",
      TRUE ~ "Databases: 1"
    )
  )

# Create the boxplot
p <- ggplot(long_data, aes(x = cell_type, y = intensity, fill = cell_type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  
  # Add protein counts as text labels
  geom_text(data = summary_stats, 
            aes(x = cell_type, y = 10^1, label = protein_count, fill = NULL),
            color = "black", size = 3.5, fontface = "bold") +
  
  # Add database count labels
  geom_text(data = database_count,
            aes(x = cell_type, y = 10^0.5, label = db_label, fill = NULL),
            color = "black", size = 3.5) +
  
  # Formatting
  scale_y_log10() +
  scale_fill_discrete(name = "Cell Type") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Cell Type", 
    y = "Log₁₀(Z-score)", 
    title = "Protein Intensity Distribution by Cell Type",
    subtitle = paste("Numbers show protein counts;", 
                     "labels indicate database coverage")
  )

# Create output directory if it doesn't exist
output_dir <- "outputs/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# === Save as TIFF ===
tiff_file <- file.path(output_dir, "boxplot_celltypes.tiff")
cat("Saving plot to:", tiff_file, "\n")

tiff(tiff_file, width = 14, height = 8, units = "in", res = 600, compression = "lzw")

# Generate the plot for saving
print(p)

# Close TIFF device
dev.off()
cat("Plot saved successfully!\n")

# === Display plot in R console/RStudio ===
print(p)

# Print detailed summary statistics
cat("\n=== Detailed Summary Statistics ===\n")
print(summary_stats)

cat("\n=== Database Coverage ===\n")
print(database_count)

cat("\n=== Overall Statistics ===\n")
cat("Total data points:", nrow(long_data), "\n")
cat("Cell types analyzed:", length(unique(long_data$cell_type)), "\n")
cat("Intensity range: 10^", 
    round(log10(min(long_data$intensity, na.rm = TRUE)), 2), 
    " to 10^", round(log10(max(long_data$intensity, na.rm = TRUE)), 2), "\n")
