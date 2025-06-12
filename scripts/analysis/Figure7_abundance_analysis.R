#!/usr/bin/env Rscript

# Load utilities and set up paths
source("scripts/utilities/load_packages.R")
ensure_output_dirs()


# Figure 7: Protein Abundance Analysis 
# Working directly with raw database files for quantitative abundance analysis

# Set CRAN mirror
if(is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

# Function to install packages if not available
install_if_missing <- function(packages) {
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing", pkg, "..."))
      install.packages(pkg, dependencies = TRUE)
      if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
        stop(paste("Failed to install package:", pkg))
      }
    }
  }
}

# Install required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "ggridges", "viridis", 
                      "scales", "gridExtra", "stringr", "ggbeeswarm")
install_if_missing(required_packages)

# Create output directory
output_dir <- "outputs/plots/03_Abundance_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create a simple local mapping function instead of external dependencies
message("Using simplified ID mapping system...")
enhanced_fast_map_to_gene_symbol <- function(accessions, descriptions = NULL) {
  if (length(accessions) == 0) return(character(0))
  
  gene_symbols <- character(length(accessions))
  
  for(i in seq_along(accessions)) {
    acc <- as.character(accessions[i])
    desc <- if(!is.null(descriptions)) as.character(descriptions[i]) else ""
    
    # Simple gene symbol extraction
    gene_symbol <- NA
    
    # If it's already a gene symbol
    if (grepl("^[A-Z][A-Z0-9_-]{1,14}$", acc) && !grepl("^ENS|^[OPQ][0-9]|^[A-Z][0-9]", acc)) {
      gene_symbol <- toupper(acc)
    }
    
    # Extract from description using GN= pattern
    if (is.na(gene_symbol) && !is.na(desc) && desc != "") {
      if(grepl("GN=([A-Za-z0-9_-]+)", desc)) {
        gene_symbol <- toupper(gsub(".*GN=([A-Za-z0-9_-]+).*", "\\1", desc))
      }
    }
    
    gene_symbols[i] <- gene_symbol
  }
  
  return(gene_symbols)
}

# Function to parse cell type from column name (same as Figure 6)
parse_cell_type <- function(column_name) {
  col <- as.character(column_name)
  
  # PXD004352 format: Intensity_B.memory_01_activated -> B
  if(grepl("^Intensity_([A-Za-z0-9]+)", col)) {
    cell_type <- gsub("^Intensity_([A-Za-z0-9]+).*", "\\1", col)
    
    # Convert pDC to DC
    if(cell_type == "pDC") cell_type <- "DC"
    
    # Combine CD4+ T cell subsets (T4, Th1, Th2, Th17, mTregs, nTregs)
    if(cell_type %in% c("T4", "Th1", "Th2", "Th17", "mTregs", "nTregs")) {
      cell_type <- "CD4"
    }
    
    # Combine CD8+ T cell subsets (T8)
    if(cell_type == "T8") {
      cell_type <- "CD8"
    }
    
    # Combine dendritic cell subsets (mDC with DC)
    if(cell_type == "mDC") {
      cell_type <- "DC"
    }
    
    # Rename monocytes for clarity  
    if(cell_type == "MO") cell_type <- "MONOCYTES"
    
    return(toupper(cell_type))
  }
  
  # PXD025174 format: Copy_Number_CD4 -> CD4
  if(grepl("^Copy_Number_([A-Za-z0-9]+)", col)) {
    cell_type <- gsub("^Copy_Number_([A-Za-z0-9]+)", "\\1", col)
    return(toupper(cell_type))
  }
  
  # PXD040957_CD8 format: LFQ.intensity.CD8TCellnaive_03 -> CD8
  if(grepl("^LFQ\\.intensity\\.CD8", col)) {
    return("CD8")
  }
  
  # PXD040957_Macrophages format: LFQ.intensity.Human_Macs_D1_M1 -> MACROPHAGES
  if(grepl("^LFQ\\.intensity\\.Human_Macs", col)) {
    return("MACROPHAGES")
  }
  
  return(NA)
}

# Function to read and process database files for abundance analysis
read_abundance_data <- function() {
  message("Reading abundance data from raw database files...")
  
  abundance_data <- list()
  
  # 1. PeptideAtlas - normalized PSMs
  message("  - Reading PeptideAtlas abundance data...")
  peptideatlas <- read.csv(get_data_path("PeptideAtlas.csv"), stringsAsFactors = FALSE)
  
  # Map to gene symbols using enhanced mapping
  gene_symbols <- enhanced_fast_map_to_gene_symbol(peptideatlas$biosequence_accession, 
                                                  peptideatlas$biosequence_desc)
  
  abundance_data$PeptideAtlas <- data.frame(
    gene_symbol = gene_symbols,
    abundance = as.numeric(peptideatlas$norm_PSMs_per_100K),
    log_abundance = log10(as.numeric(peptideatlas$norm_PSMs_per_100K) + 1),
    database = "PeptideAtlas",
    cell_type = "Plasma",
    abundance_type = "Normalized PSMs per 100K",
    stringsAsFactors = FALSE
  )
  
  # 2. PaxDb files - protein abundance in ppm
  message("  - Reading PaxDb abundance data...")
  paxdb_files <- list.files(pattern = "^PaxDb_.*\\.csv$")
  for(file in paxdb_files) {
    cell_type <- gsub("PaxDb_|\\.csv", "", file)
    message(paste("    Reading", file, "for", cell_type))
    
    lines <- readLines(file)
    data_start <- which(grepl("^string_external_id,abundance", lines))
    if(length(data_start) > 0) {
      # Read from the header line directly
      temp_data <- read.csv(file, skip = data_start - 1, stringsAsFactors = FALSE, header = TRUE)
      
      # Verify we have the expected columns
      if(!"string_external_id" %in% colnames(temp_data) || !"abundance" %in% colnames(temp_data)) {
        message(paste("      Warning: Missing expected columns in", file))
        next
      }
      
      # Clean protein IDs
      protein_ids <- gsub("^9606\\.", "", temp_data$string_external_id)
      
      # Only keep rows with valid protein IDs and abundance values
      valid_rows <- !is.na(protein_ids) & protein_ids != "" & 
                   !is.na(temp_data$abundance) & temp_data$abundance > 0
      
      if(sum(valid_rows) > 0) {
        # Map to gene symbols using enhanced mapping
        gene_symbols <- enhanced_fast_map_to_gene_symbol(protein_ids[valid_rows])
        
        abundance_data[[paste0("PaxDb_", cell_type)]] <- data.frame(
          gene_symbol = gene_symbols,
          abundance = temp_data$abundance[valid_rows],
          log_abundance = log10(temp_data$abundance[valid_rows] + 1),
          database = "PaxDb",
          cell_type = cell_type,
          abundance_type = "ppm",
          stringsAsFactors = FALSE
        )
      } else {
        message(paste("      Warning: No valid data found in", file))
      }
    } else {
      message(paste("      Warning: Could not find data header in", file))
    }
  }
  
  # 3. HPA concentration data
  message("  - Reading HPA concentration data...")
  hpa_files <- list.files(pattern = "^HPA_.*\\.csv$")
  for(file in hpa_files) {
    technique <- gsub("HPA_|\\.csv", "", file)
    message(paste("    Reading", file, "for", technique))
    
    hpa_data <- read.csv(file, stringsAsFactors = FALSE)
    
    # The HPA file has Mean.concentration column, filter valid rows
    if("Mean.concentration" %in% colnames(hpa_data)) {
      concentration_col <- "Mean.concentration"
    } else if("Concentration" %in% colnames(hpa_data)) {
      concentration_col <- "Concentration"
    } else {
      message(paste("      Warning: No concentration column found in", file))
      next
    }
    
    # Filter rows with valid gene names and concentrations
    valid_rows <- !is.na(hpa_data$Gene) & hpa_data$Gene != "" & 
                 !is.na(hpa_data[[concentration_col]]) & hpa_data[[concentration_col]] != ""
    
    if(sum(valid_rows) > 0) {
      valid_data <- hpa_data[valid_rows, ]
      
      # Convert concentration to numeric (remove units)
      # Handle encoding issues by converting to character and cleaning
      conc_strings <- as.character(valid_data[[concentration_col]])
      # Convert to UTF-8 and handle invalid characters
      conc_strings <- iconv(conc_strings, to = "UTF-8", sub = "")
      # Extract numeric values
      concentration_values <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", conc_strings)))
      
      # Remove rows where concentration conversion failed
      numeric_valid <- !is.na(concentration_values) & concentration_values > 0
      
      if(sum(numeric_valid) > 0) {
        # HPA already provides gene symbols
        gene_symbols <- toupper(valid_data$Gene[numeric_valid])
        
        abundance_data[[paste0("HPA_", technique)]] <- data.frame(
          gene_symbol = gene_symbols,
          abundance = concentration_values[numeric_valid],
          log_abundance = log10(concentration_values[numeric_valid] + 0.001),  # Add small value for log transform
          database = "HPA",
          cell_type = technique,
          abundance_type = "mg/L",
          stringsAsFactors = FALSE
        )
      } else {
        message(paste("      Warning: No valid concentration data in", file))
      }
    } else {
      message(paste("      Warning: No valid data found in", file))
    }
  }
  
  # 4. GPMDB spectral counts
  message("  - Reading GPMDB spectral count data...")
  gpmdb_files <- list.files(pattern = "^GPMDB_.*\\.csv$")
  for(file in gpmdb_files) {
    cell_type <- gsub("GPMDB_|\\.csv", "", file)
    message(paste("    Reading", file, "for", cell_type))
    
    lines <- readLines(file)
    header_line <- which(grepl("^#,accession,total", lines))
    if(length(header_line) > 0) {
      # Skip to the line after the header
      gpmdb_data <- read.csv(file, skip = header_line, stringsAsFactors = FALSE, header = FALSE)
      
      # Set column names manually
      colnames(gpmdb_data) <- c("rank", "accession", "total", "log_e", "EC", "description")
      
      # Filter valid rows
      valid_rows <- !is.na(gpmdb_data$accession) & gpmdb_data$accession != "" & 
                   !is.na(gpmdb_data$total) & gpmdb_data$total > 0
      
      if(sum(valid_rows) > 0) {
        # Map GPMDB protein accessions to gene symbols using enhanced mapping
        gene_symbols <- enhanced_fast_map_to_gene_symbol(gpmdb_data$accession[valid_rows], 
                                                        gpmdb_data$description[valid_rows])
        
        abundance_data[[paste0("GPMDB_", cell_type)]] <- data.frame(
          gene_symbol = gene_symbols,
          abundance = gpmdb_data$total[valid_rows],
          log_abundance = log10(gpmdb_data$total[valid_rows] + 1),
          database = "GPMDB",
          cell_type = cell_type,
          abundance_type = "Spectral Counts",
          stringsAsFactors = FALSE
        )
      } else {
        message(paste("      Warning: No valid data found in", file))
      }
    } else {
      message(paste("      Warning: Could not find data header in", file))
    }
  }
  
  # 5. ProteomeXchange files - with proper cell type parsing
  message("  - Reading ProteomeXchange data with cell type parsing...")
  pxd_files <- list.files(pattern = "^PXD.*\\.csv$")
  for(file in pxd_files) {
    message(paste("    Processing", file, "for cell type-specific abundance"))
    
    # Read data with header (not skip = 1)
    pxd_data <- read.csv(file, stringsAsFactors = FALSE)
    
    # Find gene names column
    gene_col <- NULL
    if("Gene names" %in% colnames(pxd_data)) {
      gene_col <- "Gene names"
    } else if("Gene.names" %in% colnames(pxd_data)) {
      gene_col <- "Gene.names"
    } else if("Gene_names" %in% colnames(pxd_data)) {
      gene_col <- "Gene_names"
    } else if("Gene.Name" %in% colnames(pxd_data)) {
      gene_col <- "Gene.Name"
    }
    
    if(!is.null(gene_col)) {
      # Find intensity/abundance columns and their cell types
      intensity_cols <- grep("(LFQ\\.intensity|Intensity_|Copy_Number_)", colnames(pxd_data), value = TRUE)
      
      if(length(intensity_cols) > 0) {
        # Parse cell types from column names
        cell_type_map <- data.frame(
          column = intensity_cols,
          cell_type = sapply(intensity_cols, parse_cell_type),
          stringsAsFactors = FALSE
        )
        
        # Remove columns with unknown cell types
        cell_type_map <- cell_type_map[!is.na(cell_type_map$cell_type), ]
        
        if(nrow(cell_type_map) > 0) {
          # Group by cell type and average intensities for same gene
          unique_cell_types <- unique(cell_type_map$cell_type)
          
          for(cell_type in unique_cell_types) {
            cols_for_celltype <- cell_type_map$column[cell_type_map$cell_type == cell_type]
            
            if(length(cols_for_celltype) > 0) {
              # Calculate mean intensity across all columns for this cell type
              if(length(cols_for_celltype) == 1) {
                mean_intensity <- pxd_data[[cols_for_celltype]]
              } else {
                mean_intensity <- rowMeans(pxd_data[, cols_for_celltype, drop = FALSE], na.rm = TRUE)
              }
              
              valid_genes <- character()
              valid_intensities <- numeric()
              
              for(i in 1:nrow(pxd_data)) {
                gene_str <- as.character(pxd_data[[gene_col]][i])
                intensity <- mean_intensity[i]
                
                # Check for single gene and valid intensity
                if(!is.na(gene_str) && gene_str != "" && !grepl("[;,]", gene_str) &&
                   !is.na(intensity) && intensity > 0) {
                  valid_genes <- c(valid_genes, toupper(trimws(gene_str)))
                  valid_intensities <- c(valid_intensities, intensity)
                }
              }
              
              if(length(valid_genes) > 0) {
                abundance_data[[paste0("PXD_", cell_type, "_", gsub("\\.csv", "", file))]] <- data.frame(
                  gene_symbol = valid_genes,
                  abundance = valid_intensities,
                  log_abundance = log10(valid_intensities + 1),
                  database = "ProteomeXchange",
                  cell_type = cell_type,
                  abundance_type = "LFQ Intensity",
                  stringsAsFactors = FALSE
                )
                message(paste("      ", cell_type, "- Kept", length(valid_genes), "gene abundance measurements"))
              }
            }
          }
        }
      }
    }
  }
  
  # Combine all data
  all_abundance <- do.call(rbind, abundance_data)
  
  # Aggregate by gene symbol (sum abundances for same gene from same database/cell_type)
  aggregated_abundance <- all_abundance %>%
    group_by(gene_symbol, database, cell_type, abundance_type) %>%
    summarise(
      abundance = sum(abundance, na.rm = TRUE),
      log_abundance = log10(sum(abundance, na.rm = TRUE) + 0.001),
      .groups = 'drop'
    ) %>%
    filter(!is.na(gene_symbol) & gene_symbol != "" & 
           abundance > 0 & is.finite(abundance))
  
  message(paste("Total gene abundance measurements:", nrow(aggregated_abundance)))
  message(paste("Unique genes:", length(unique(aggregated_abundance$gene_symbol))))
  return(aggregated_abundance)
}

# Function to calculate abundance statistics
calculate_abundance_stats <- function(abundance_data) {
  message("Calculating abundance statistics...")
  
  stats_summary <- abundance_data %>%
    group_by(database, cell_type, abundance_type) %>%
    summarise(
      n_proteins = n(),
      min_abundance = min(abundance, na.rm = TRUE),
      max_abundance = max(abundance, na.rm = TRUE),
      median_abundance = median(abundance, na.rm = TRUE),
      mean_abundance = mean(abundance, na.rm = TRUE),
      sd_abundance = sd(abundance, na.rm = TRUE),
      q25 = quantile(abundance, 0.25, na.rm = TRUE),
      q75 = quantile(abundance, 0.75, na.rm = TRUE),
      dynamic_range = log10(max_abundance / min_abundance),
      .groups = 'drop'
    )
  
  return(stats_summary)
}

# Read abundance data
abundance_data <- read_abundance_data()

# Calculate statistics
abundance_stats <- calculate_abundance_stats(abundance_data)

# Plot 1: Abundance Distribution Ridges by Database
message("Creating abundance distribution ridges plot...")
tiff(file.path(output_dir, "abundance_distributions_by_database.tiff"), 
     width = 14, height = 10, units = "in", res = 600, compression = "lzw")

ggplot(abundance_data, aes(x = log_abundance, y = database, fill = database)) +
  geom_density_ridges(alpha = 0.7, scale = 2) +
  scale_fill_viridis_d(name = "Database") +
  theme_minimal() +
  labs(title = "Protein Abundance Distributions by Database",
       subtitle = "Log10-transformed abundance values showing dynamic range differences",
       x = "Log10(Abundance + offset)", y = "Database") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12))

dev.off()

# Plot 2: Abundance Distributions by Cell Type (PaxDb only for comparison)
message("Creating PaxDb cell type abundance comparison...")
paxdb_data <- abundance_data %>% filter(database == "PaxDb")

tiff(file.path(output_dir, "paxdb_celltype_abundance_comparison.tiff"), 
     width = 14, height = 10, units = "in", res = 600, compression = "lzw")

ggplot(paxdb_data, aes(x = log_abundance, y = reorder(cell_type, log_abundance, median), 
                       fill = cell_type)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  scale_fill_viridis_d(name = "Cell Type") +
  theme_minimal() +
  labs(title = "PaxDb Protein Abundance Across Cell Types",
       subtitle = "Comparing protein abundance profiles between blood cell types",
       x = "Log10(Abundance in ppm)", y = "Cell Type") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"))

dev.off()

# Plot 3: HPA Concentration Analysis
message("Creating HPA concentration analysis...")
hpa_data <- abundance_data %>% filter(database == "HPA")

if(nrow(hpa_data) > 0) {
  tiff(file.path(output_dir, "hpa_concentration_analysis.tiff"), 
       width = 12, height = 8, units = "in", res = 600, compression = "lzw")
  
  p1 <- ggplot(hpa_data, aes(x = cell_type, y = abundance, fill = cell_type)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
    scale_y_log10(labels = label_number(suffix = " mg/L")) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "HPA Protein Concentrations by Method",
         x = "Detection Method", y = "Concentration (mg/L, log scale)")
  
  p2 <- ggplot(hpa_data, aes(x = abundance, fill = cell_type)) +
    geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
    scale_x_log10(labels = label_number(suffix = " mg/L")) +
    scale_fill_brewer(type = "qual", palette = "Set2", name = "Method") +
    theme_minimal() +
    labs(title = "Distribution of HPA Concentrations",
         x = "Concentration (mg/L, log scale)", y = "Number of Proteins")
  
  grid.arrange(p1, p2, nrow = 2, 
               top = "Human Protein Atlas (HPA) Concentration Analysis")
  
  dev.off()
}

# Plot 4: Dynamic Range Comparison
message("Creating dynamic range comparison...")
tiff(file.path(output_dir, "abundance_dynamic_range_comparison.tiff"), 
     width = 12, height = 8, units = "in", res = 600, compression = "lzw")

abundance_stats_viz <- abundance_stats %>%
  mutate(db_celltype = paste0(database, " (", cell_type, ")"))

ggplot(abundance_stats_viz, aes(x = reorder(db_celltype, dynamic_range), 
                                y = dynamic_range, fill = database)) +
  geom_col(alpha = 0.8) +
  coord_flip() +
  scale_fill_viridis_d(name = "Database") +
  theme_minimal() +
  labs(title = "Protein Abundance Dynamic Range by Database/Cell Type",
       subtitle = "Higher values indicate wider range of protein abundances",
       x = "Database - Cell Type", y = "Dynamic Range (log10 orders of magnitude)") +
  geom_text(aes(label = round(dynamic_range, 1)), 
            hjust = -0.1, size = 3)

dev.off()

# Plot 5: Abundance vs Coverage Analysis
message("Creating abundance vs coverage analysis...")
coverage_abundance <- abundance_data %>%
  group_by(database, cell_type) %>%
  summarise(
    n_proteins = n(),
    median_abundance = median(abundance, na.rm = TRUE),
    mean_log_abundance = mean(log_abundance, na.rm = TRUE),
    abundance_range = max(abundance, na.rm = TRUE) - min(abundance, na.rm = TRUE),
    .groups = 'drop'
  )

tiff(file.path(output_dir, "abundance_vs_coverage.tiff"), 
     width = 10, height = 8, units = "in", res = 600, compression = "lzw")

ggplot(coverage_abundance, aes(x = n_proteins, y = mean_log_abundance, 
                              color = database, size = abundance_range)) +
  geom_point(alpha = 0.8) +
  scale_color_viridis_d(name = "Database") +
  scale_size_continuous(name = "Abundance\nRange", 
                       labels = label_scientific()) +
  theme_minimal() +
  labs(title = "Protein Coverage vs Average Abundance",
       subtitle = "Relationship between number of proteins and their average abundance",
       x = "Number of Proteins Detected", 
       y = "Mean Log10(Abundance)") +
  geom_text(aes(label = cell_type), 
            hjust = 0, vjust = 0, size = 3, show.legend = FALSE)

dev.off()

# Plot 6: High Abundance Proteins Analysis
message("Creating high abundance proteins analysis...")
# Identify top 10% most abundant proteins in each dataset
high_abundance <- abundance_data %>%
  group_by(database, cell_type) %>%
  mutate(abundance_percentile = percent_rank(abundance)) %>%
  filter(abundance_percentile >= 0.9) %>%
  ungroup()

# Count how often proteins appear in high abundance across datasets
high_abundance_counts <- high_abundance %>%
  count(gene_symbol, name = "high_abundance_count") %>%
  arrange(desc(high_abundance_count))

# Top consistently high abundance proteins
top_proteins <- high_abundance_counts %>%
  head(20) %>%
  left_join(high_abundance %>% 
             select(gene_symbol, abundance) %>%
             group_by(gene_symbol) %>%
             summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = 'drop'),
           by = "gene_symbol")

tiff(file.path(output_dir, "consistently_high_abundance_proteins.tiff"), 
     width = 12, height = 8, units = "in", res = 600, compression = "lzw")

  ggplot(top_proteins, aes(x = reorder(gene_symbol, high_abundance_count), 
                        y = high_abundance_count)) +
  geom_col(fill = "darkred", alpha = 0.7) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Consistently High Abundance Proteins",
       subtitle = "Proteins appearing in top 10% abundance across multiple datasets",
       x = "Protein ID", y = "Number of Datasets in Top 10%") +
  geom_text(aes(label = high_abundance_count), hjust = -0.1, size = 3)

dev.off()

# Save detailed statistics
write.csv(abundance_stats, file.path(output_dir, "abundance_statistics_summary.csv"), 
          row.names = FALSE)

write.csv(high_abundance_counts, file.path(output_dir, "high_abundance_proteins.csv"), 
          row.names = FALSE)

# Write summary report
sink(file.path(output_dir, "abundance_analysis_summary.txt"))
cat("Protein Abundance Analysis Summary\n")
cat("==================================\n\n")

cat("Data Summary:\n")
cat(paste("Total abundance measurements:", nrow(abundance_data), "\n"))
cat(paste("Unique genes with abundance data:", length(unique(abundance_data$gene_symbol)), "\n"))
cat(paste("Databases analyzed:", length(unique(abundance_data$database)), "\n"))
cat(paste("Cell types/conditions analyzed:", length(unique(abundance_data$cell_type)), "\n\n"))

cat("Abundance Ranges by Database:\n")
for(db in unique(abundance_stats$database)) {
  db_stats <- abundance_stats[abundance_stats$database == db, ]
  cat(paste("  ", db, ":\n"))
  cat(paste("    Dynamic range:", round(mean(db_stats$dynamic_range), 2), "orders of magnitude\n"))
  cat(paste("    Median proteins per dataset:", round(median(db_stats$n_proteins)), "\n"))
  cat(paste("    Abundance type:", unique(db_stats$abundance_type)[1], "\n\n"))
}

cat("Top 10 Most Consistently High Abundance Proteins:\n")
for(i in 1:min(10, nrow(high_abundance_counts))) {
  gene <- high_abundance_counts$gene_symbol[i]
  count <- high_abundance_counts$high_abundance_count[i]
      cat(paste("  ", i, ".", gene, "- appears in", count, "datasets\n"))
}

sink()

message("Abundance analysis completed!")
message(paste("Results saved to:", output_dir)) 
