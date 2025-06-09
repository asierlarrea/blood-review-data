# ============================================================================
# Figure 3: Cell Type Proteome Bubble Plot
# Description: Creates a bubble plot showing protein counts by database and cell type
# ============================================================================

# Set CRAN mirror for package installation
if (length(getOption("repos")) == 0 || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

# Load required libraries
if (!require(ggplot2, quietly = TRUE)) {
  cat("Installing ggplot2 package...\n")
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require(scales, quietly = TRUE)) {
  cat("Installing scales package...\n")
  install.packages("scales", dependencies = TRUE)
  library(scales)
}

# Load data with error handling
data_file <- "cell_types.csv"
if (!file.exists(data_file)) {
  stop(paste("Error: File", data_file, "not found in current directory."))
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
  scale_color_discrete(name = "Database") +
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

# Display the plot
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

