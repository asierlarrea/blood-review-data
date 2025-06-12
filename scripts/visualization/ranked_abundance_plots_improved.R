#!/usr/bin/env Rscript

# Simplified Ranked Abundance Plots
# Creates basic ranked abundance visualization

# Load required packages
source("scripts/utilities/load_packages.R")
required_packages <- c("ggplot2", "dplyr", "gridExtra")
load_packages(required_packages)

# Create output directory
output_dir <- "outputs/plots/05_Ranked_Abundance_Improved"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

message("Creating simplified ranked abundance plots...")

# Generate synthetic data for demonstration
set.seed(42)
n_proteins <- 200

demo_data <- data.frame(
  gene_symbol = paste0("GENE", 1:n_proteins),
  abundance = exp(rnorm(n_proteins, mean = 5, sd = 2)),
  source = "Demo Data",
  abundance_type = "Arbitrary Units"
) %>%
  arrange(desc(abundance)) %>%
  mutate(
    rank = row_number(),
    log_abundance = log10(abundance),
    # Highlight some "biomarkers"
    is_biomarker = gene_symbol %in% paste0("GENE", c(5, 12, 23, 45, 67, 89, 123, 156))
  )

# Create ranked abundance plot
p1 <- ggplot(demo_data, aes(x = rank, y = log_abundance)) +
  geom_area(fill = "lightgrey", alpha = 0.7) +
  geom_point(data = demo_data[demo_data$is_biomarker, ], 
             aes(x = rank, y = log_abundance), 
             color = "red", size = 2, alpha = 0.8) +
  scale_y_continuous(
    name = "Concentration (log10 scale)",
    labels = function(x) parse(text = paste0("10^", round(x, 1)))
  ) +
  scale_x_continuous(
    name = "Protein Rank (highest to lowest)"
  ) +
  theme_minimal() +
  labs(
    title = "Ranked Protein Abundance (Demonstration)",
    subtitle = paste0("n = ", n_proteins, " proteins | Biomarkers highlighted in red")
  )

# Save plot
ggsave(file.path(output_dir, "ranked_abundance_demo.png"), p1, 
       width = 12, height = 8, dpi = 300)

# Create summary statistics
summary_stats <- data.frame(
  source = "Demo Data",
  n_proteins = n_proteins,
  biomarkers_found = sum(demo_data$is_biomarker),
  min_abundance = min(demo_data$abundance),
  max_abundance = max(demo_data$abundance),
  dynamic_range = log10(max(demo_data$abundance) / min(demo_data$abundance))
)

write.csv(summary_stats, file.path(output_dir, "ranked_abundance_summary_demo.csv"), 
          row.names = FALSE)

message("✓ Simplified ranked abundance analysis completed")
message(paste("✓ Output saved to:", output_dir))
message("✓ Created demonstration ranked abundance plot") 