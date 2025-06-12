#!/usr/bin/env Rscript

# Load utilities and set up paths
source("scripts/utilities/load_packages.R")
ensure_output_dirs()


# Figure 9: Comparative Analysis Across Cell Types and Databases
# Working directly with raw database files for cross-database protein comparison

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
required_packages <- c("ggplot2", "dplyr", "tidyr", "VennDiagram", "UpSetR", 
                      "pheatmap", "scales", "gridExtra", "stringr")
install_if_missing(required_packages)

# Create output directory
output_dir <- "outputs/plots/08_Database_Comparison"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to map protein IDs to gene symbols
map_to_gene_symbol <- function(protein_ids, descriptions = NULL) {
  gene_symbols <- character(length(protein_ids))
  
  for(i in seq_along(protein_ids)) {
    protein_id <- as.character(protein_ids[i])
    description <- if(!is.null(descriptions)) as.character(descriptions[i]) else ""
    
    # Initialize with the protein ID as fallback
    gene_symbol <- protein_id
    
    # Extract gene symbol from description if available
    if(!is.null(descriptions) && !is.na(description) && description != "") {
      # Look for gene symbols in various formats in descriptions
      # Pattern 1: "GENE_NAME," or "GENE_NAME " at the beginning
      if(grepl("^[A-Z][A-Z0-9_-]+[,\\s]", description)) {
        gene_symbol <- gsub("^([A-Z][A-Z0-9_-]+)[,\\s].*", "\\1", description)
      }
      # Pattern 2: Gene symbols in brackets or after Source:HGNC Symbol
      else if(grepl("Source:HGNC Symbol;Acc:(HGNC:)?([A-Z0-9_-]+)", description)) {
        gene_symbol <- gsub(".*Source:HGNC Symbol;Acc:(?:HGNC:)?([A-Z0-9_-]+).*", "\\1", description, perl = TRUE)
      }
      # Pattern 3: Gene symbols at the start before comma
      else if(grepl("^[A-Z0-9_-]+,", description)) {
        gene_symbol <- gsub("^([A-Z0-9_-]+),.*", "\\1", description)
      }
    }
    
    # If it's still a protein accession, try to extract gene from common patterns
    if(grepl("^(ENSP|P\\d|Q\\d|O\\d)", gene_symbol)) {
      # Common gene name patterns from the protein ID itself
      gene_symbol <- paste0("GENE_", substr(protein_id, 1, 8))
    }
    
    # Clean the gene symbol
    gene_symbol <- gsub("[^A-Z0-9_-]", "", toupper(gene_symbol))
    
    # Ensure we have a valid gene symbol
    if(gene_symbol == "" || nchar(gene_symbol) < 2) {
      gene_symbol <- paste0("UNKNOWN_", substr(protein_id, 1, 8))
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

# Function to read protein data for comparative analysis
read_comparative_data <- function() {
  message("Reading protein data for comparative analysis...")
  
  protein_data <- list()
  
  # 1. PeptideAtlas
  message("  - Reading PeptideAtlas...")
  peptideatlas <- read.csv(get_data_path("PeptideAtlas.csv"), stringsAsFactors = FALSE)
  
  # Map to gene symbols using descriptions
  gene_symbols <- map_to_gene_symbol(peptideatlas$biosequence_accession, 
                                    peptideatlas$biosequence_desc)
  
  protein_data$PeptideAtlas <- data.frame(
    gene_symbol = gene_symbols,
    database = "PeptideAtlas",
    cell_type = "Plasma",
    abundance = as.numeric(peptideatlas$norm_PSMs_per_100K),
    stringsAsFactors = FALSE
  )
  
  # 2. PaxDb files - multiple cell types
  message("  - Reading PaxDb files...")
  paxdb_files <- list.files(pattern = "^PaxDb_.*\\.csv$")
  for(file in paxdb_files) {
    cell_type <- gsub("PaxDb_|\\.csv", "", file)
    message(paste("    Reading", file, "for", cell_type))
    
    lines <- readLines(file)
    data_start <- which(grepl("^#string_external_id,abundance", lines))
    if(length(data_start) > 0) {
      # Read from the line after the header
      temp_data <- read.csv(file, skip = data_start, stringsAsFactors = FALSE, header = FALSE)
      
      # Set column names based on actual number of columns
      if(ncol(temp_data) == 2) {
        colnames(temp_data) <- c("string_external_id", "abundance")
      } else if(ncol(temp_data) == 3) {
        colnames(temp_data) <- c("string_external_id", "abundance", "raw_spectral_count")
      } else {
        message(paste("      Warning: Unexpected number of columns in", file, ":", ncol(temp_data)))
        next
      }
      
      # Clean protein IDs
      gene_symbols <- gsub("^9606\\.", "", temp_data$string_external_id)
      
      # Only keep rows with valid protein IDs and abundance values
      valid_rows <- !is.na(gene_symbols) & gene_symbols != "" & 
                   !is.na(temp_data$abundance) & temp_data$abundance > 0
      
      if(sum(valid_rows) > 0) {
        protein_data[[paste0("PaxDb_", cell_type)]] <- data.frame(
          gene_symbol = gene_symbols[valid_rows],
          database = "PaxDb",
          cell_type = cell_type,
          abundance = temp_data$abundance[valid_rows],
          stringsAsFactors = FALSE
        )
      } else {
        message(paste("      Warning: No valid data found in", file))
      }
    } else {
      message(paste("      Warning: Could not find data header in", file))
    }
  }
  
  # 3. HPA files
  message("  - Reading HPA files...")
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
        protein_data[[paste0("HPA_", technique)]] <- data.frame(
          gene_symbol = valid_data$Gene[numeric_valid],
          database = "HPA",
          cell_type = technique,
          abundance = concentration_values[numeric_valid],
          stringsAsFactors = FALSE
        )
      } else {
        message(paste("      Warning: No valid concentration data in", file))
      }
    } else {
      message(paste("      Warning: No valid data found in", file))
    }
  }
  
  # 4. GPMDB files - multiple cell types
  message("  - Reading GPMDB files...")
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
        protein_data[[paste0("GPMDB_", cell_type)]] <- data.frame(
          gene_symbol = gpmdb_data$accession[valid_rows],
          database = "GPMDB",
          cell_type = cell_type,
          abundance = gpmdb_data$total[valid_rows],
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
    message(paste("    Processing", file, "for cell type-specific comparative analysis"))
    
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
                protein_data[[paste0("PXD_", cell_type, "_", gsub("\\.csv", "", file))]] <- data.frame(
                  gene_symbol = valid_genes,
                  database = "ProteomeXchange",
                  cell_type = cell_type,
                  abundance = valid_intensities,
                  stringsAsFactors = FALSE
                )
                message(paste("      ", cell_type, "- Kept", length(valid_genes), "gene entries"))
              }
            }
          }
        }
      }
    }
  }
  
  # Combine all data
  all_data <- do.call(rbind, protein_data)
  
  # Clean protein IDs
  all_data$gene_symbol <- sapply(all_data$gene_symbol, function(x) {
    first_id <- strsplit(as.character(x), ";")[[1]][1]
    main_id <- gsub("-\\d+$", "", first_id)
    return(main_id)
  })
  
  # Remove invalid entries
  all_data <- all_data[!is.na(all_data$gene_symbol) & 
                      all_data$gene_symbol != "" &
                      !is.na(all_data$abundance) &
                      is.finite(all_data$abundance) &
                      all_data$abundance > 0, ]
  
  # Create cell type groups
  all_data$cell_group <- case_when(
    grepl("CD4|CD8", all_data$cell_type, ignore.case = TRUE) ~ "T Cells",
    grepl("B_Cell|B-Cell", all_data$cell_type, ignore.case = TRUE) ~ "B Cells",
    grepl("NK", all_data$cell_type, ignore.case = TRUE) ~ "NK Cells",
    grepl("Monocyte|Macrophage", all_data$cell_type, ignore.case = TRUE) ~ "Myeloid",
    grepl("Neutrophil|Eosinophil|Basophil", all_data$cell_type, ignore.case = TRUE) ~ "Granulocytes",
    grepl("Platelet", all_data$cell_type, ignore.case = TRUE) ~ "Platelets",
    grepl("Erythrocyte", all_data$cell_type, ignore.case = TRUE) ~ "Erythrocytes",
    grepl("Plasma", all_data$cell_type, ignore.case = TRUE) ~ "Plasma",
    TRUE ~ "Other"
  )
  
  # Mark immune vs non-immune cells
  all_data$cell_category <- ifelse(
    all_data$cell_group %in% c("T Cells", "B Cells", "NK Cells", "Myeloid", "Granulocytes"),
    "Immune", "Non-Immune"
  )
  
  message(paste("Total protein measurements loaded:", nrow(all_data)))
  message(paste("Unique proteins:", length(unique(all_data$gene_symbol))))
  message(paste("Cell types:", length(unique(all_data$cell_type))))
  
  return(all_data)
}

# Function to create protein overlap matrix for UpSet plot
create_upset_data <- function(protein_data, by_variable = "database") {
  message(paste("Creating UpSet data by", by_variable))
  
  if(by_variable == "database") {
    # Group by database
    overlap_data <- protein_data %>%
      group_by(database) %>%
      summarise(proteins = list(unique(gene_symbol)), .groups = 'drop')
    
    group_names <- overlap_data$database
    protein_lists <- overlap_data$proteins
  } else {
    # Group by cell group
    overlap_data <- protein_data %>%
      filter(cell_group != "Other") %>%
      group_by(cell_group) %>%
      summarise(proteins = list(unique(gene_symbol)), .groups = 'drop')
    
    group_names <- overlap_data$cell_group
    protein_lists <- overlap_data$proteins
  }
  
  # Create binary matrix for UpSet
  all_proteins <- unique(unlist(protein_lists))
  upset_matrix <- matrix(0, nrow = length(all_proteins), ncol = length(group_names))
  rownames(upset_matrix) <- all_proteins
  colnames(upset_matrix) <- group_names
  
  for(i in seq_along(protein_lists)) {
    proteins_in_group <- protein_lists[[i]]
    upset_matrix[proteins_in_group, i] <- 1
  }
  
  return(upset_matrix)
}

# Function to calculate Jaccard similarity
calculate_jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# Read comparative data
protein_data <- read_comparative_data()

# Plot 1: Database Protein Overlap UpSet Plot
message("Creating database overlap UpSet plot...")
db_upset_matrix <- create_upset_data(protein_data, "database")

# Validate upset matrix dimensions
if(ncol(db_upset_matrix) >= 2 && nrow(db_upset_matrix) > 0) {
  tiff(file.path(output_dir, "database_protein_overlap_upset.tiff"), 
       width = 12, height = 8, units = "in", res = 600, compression = "lzw")
  
  upset(as.data.frame(db_upset_matrix), 
        sets = colnames(db_upset_matrix),
        sets.bar.color = "#56B4E9",
        main.bar.color = "#E69F00",
        matrix.color = "#CC79A7",
        text.scale = 1.2,
        order.by = "freq")
  
  dev.off()
  message("Database UpSet plot created successfully")
} else {
  message("Warning: Insufficient data for database UpSet plot")
  message(paste("Matrix dimensions:", nrow(db_upset_matrix), "x", ncol(db_upset_matrix)))
}

# Plot 2: Cell Type Protein Overlap UpSet Plot
message("Creating cell type overlap UpSet plot...")
celltype_upset_matrix <- create_upset_data(protein_data, "cell_group")

# Validate upset matrix dimensions
if(ncol(celltype_upset_matrix) >= 2 && nrow(celltype_upset_matrix) > 0) {
  tiff(file.path(output_dir, "celltype_protein_overlap_upset.tiff"), 
       width = 12, height = 8, units = "in", res = 600, compression = "lzw")
  
  upset(as.data.frame(celltype_upset_matrix), 
        sets = colnames(celltype_upset_matrix),
        sets.bar.color = "#56B4E9",
        main.bar.color = "#E69F00", 
        matrix.color = "#CC79A7",
        text.scale = 1.2,
        order.by = "freq")
  
  dev.off()
  message("Cell type UpSet plot created successfully")
} else {
  message("Warning: Insufficient data for cell type UpSet plot")
  message(paste("Matrix dimensions:", nrow(celltype_upset_matrix), "x", ncol(celltype_upset_matrix)))
}

# Plot 3: Immune vs Non-Immune Cell Comparison
message("Creating immune vs non-immune comparison...")
immune_comparison <- protein_data %>%
  filter(cell_category %in% c("Immune", "Non-Immune")) %>%
  group_by(cell_category) %>%
  summarise(
    unique_proteins = n_distinct(gene_symbol),
    total_measurements = n(),
    mean_abundance = mean(abundance, na.rm = TRUE),
    median_abundance = median(abundance, na.rm = TRUE),
    .groups = 'drop'
  )

# Get proteins unique to each category
immune_proteins <- protein_data %>% 
  filter(cell_category == "Immune") %>% 
  pull(gene_symbol) %>% unique()

nonimmune_proteins <- protein_data %>% 
  filter(cell_category == "Non-Immune") %>% 
  pull(gene_symbol) %>% unique()

shared_proteins <- intersect(immune_proteins, nonimmune_proteins)
immune_only <- setdiff(immune_proteins, nonimmune_proteins)
nonimmune_only <- setdiff(nonimmune_proteins, immune_proteins)

# Create comparison data
comparison_data <- data.frame(
  Category = c("Immune Only", "Shared", "Non-Immune Only"),
  Count = c(length(immune_only), length(shared_proteins), length(nonimmune_only)),
  Percentage = c(
    length(immune_only) / length(union(immune_proteins, nonimmune_proteins)) * 100,
    length(shared_proteins) / length(union(immune_proteins, nonimmune_proteins)) * 100,
    length(nonimmune_only) / length(union(immune_proteins, nonimmune_proteins)) * 100
  )
)

tiff(file.path(output_dir, "immune_vs_nonimmune_comparison.tiff"), 
     width = 12, height = 8, units = "in", res = 600, compression = "lzw")

p1 <- ggplot(comparison_data, aes(x = Category, y = Count, fill = Category)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c("Immune Only" = "#E31A1C", 
                              "Shared" = "#33A02C", 
                              "Non-Immune Only" = "#1F78B4")) +
  theme_minimal() +
  labs(title = "Protein Distribution: Immune vs Non-Immune Cells",
       x = "Protein Category", y = "Number of Proteins") +
  geom_text(aes(label = paste0(Count, "\n(", round(Percentage, 1), "%)")), 
            vjust = 1.5, color = "white", size = 4) +
  theme(legend.position = "none")

p2 <- ggplot(immune_comparison, aes(x = cell_category, y = unique_proteins, 
                                   fill = cell_category)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c("Immune" = "#E31A1C", "Non-Immune" = "#1F78B4")) +
  theme_minimal() +
  labs(title = "Total Unique Proteins by Cell Category",
       x = "Cell Category", y = "Unique Proteins") +
  geom_text(aes(label = unique_proteins), vjust = 1.5, color = "white", size = 4) +
  theme(legend.position = "none")

grid.arrange(p1, p2, ncol = 2,
             top = "Comparative Analysis: Immune vs Non-Immune Cells")

dev.off()

# Plot 4: T Cell Subset Comparison (CD4 vs CD8)
message("Creating T cell subset comparison...")
tcell_data <- protein_data %>%
  filter(grepl("CD4|CD8", cell_type, ignore.case = TRUE)) %>%
  mutate(tcell_type = ifelse(grepl("CD4", cell_type, ignore.case = TRUE), "CD4", "CD8"))

if(nrow(tcell_data) > 0) {
  cd4_proteins <- tcell_data %>% filter(tcell_type == "CD4") %>% pull(gene_symbol) %>% unique()
  cd8_proteins <- tcell_data %>% filter(tcell_type == "CD8") %>% pull(gene_symbol) %>% unique()
  
  tcell_shared <- intersect(cd4_proteins, cd8_proteins)
  cd4_only <- setdiff(cd4_proteins, cd8_proteins)
  cd8_only <- setdiff(cd8_proteins, cd4_proteins)
  
  tcell_comparison <- data.frame(
    Category = c("CD4 Only", "Shared", "CD8 Only"),
    Count = c(length(cd4_only), length(tcell_shared), length(cd8_only)),
    Percentage = c(
      length(cd4_only) / length(union(cd4_proteins, cd8_proteins)) * 100,
      length(tcell_shared) / length(union(cd4_proteins, cd8_proteins)) * 100,
      length(cd8_only) / length(union(cd4_proteins, cd8_proteins)) * 100
    )
  )
  
  tiff(file.path(output_dir, "tcell_subset_comparison.tiff"), 
       width = 10, height = 6, units = "in", res = 600, compression = "lzw")
  
  ggplot(tcell_comparison, aes(x = Category, y = Count, fill = Category)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(values = c("CD4 Only" = "#FF7F00", 
                                "Shared" = "#33A02C", 
                                "CD8 Only" = "#6A3D9A")) +
    theme_minimal() +
    labs(title = "T Cell Subset Protein Comparison: CD4 vs CD8",
         subtitle = "Unique and shared proteins between T cell subsets",
         x = "Protein Category", y = "Number of Proteins") +
    geom_text(aes(label = paste0(Count, "\n(", round(Percentage, 1), "%)")), 
              vjust = 1.5, color = "white", size = 4) +
    theme(legend.position = "none")
  
  dev.off()
}

# Plot 5: Database Coverage Efficiency Analysis
message("Creating database coverage efficiency analysis...")
efficiency_data <- protein_data %>%
  group_by(database, cell_type) %>%
  summarise(
    protein_count = n_distinct(gene_symbol),
    total_measurements = n(),
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(database) %>%
  summarise(
    cell_types_covered = n(),
    total_proteins = sum(protein_count),
    unique_proteins = length(unique(protein_data$gene_symbol[protein_data$database == database[1]])),
    efficiency = unique_proteins / total_proteins,  # Ratio of unique to total
    coverage_breadth = cell_types_covered,
    .groups = 'drop'
  )

tiff(file.path(output_dir, "database_coverage_efficiency.tiff"), 
     width = 10, height = 8, units = "in", res = 600, compression = "lzw")

ggplot(efficiency_data, aes(x = coverage_breadth, y = unique_proteins, 
                           size = total_proteins, color = database)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(3, 12), name = "Total\nMeasurements") +
  scale_color_viridis_d(name = "Database") +
  theme_minimal() +
  labs(title = "Database Coverage Efficiency Analysis",
       subtitle = "Unique proteins vs cell type coverage breadth",
       x = "Number of Cell Types Covered", 
       y = "Unique Proteins Detected") +
  geom_text(aes(label = database), hjust = 0, vjust = 0, 
            size = 3, show.legend = FALSE)

dev.off()

# Plot 6: Cross-Database Protein Consistency
message("Creating cross-database consistency analysis...")
# Find proteins that appear in multiple databases
protein_db_counts <- protein_data %>%
  group_by(gene_symbol) %>%
  summarise(
    databases_count = n_distinct(database),
    databases = paste(unique(database), collapse = ", "),
    cell_types_count = n_distinct(cell_type),
    mean_abundance = mean(abundance, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(databases_count))

consistency_summary <- protein_db_counts %>%
  count(databases_count, name = "protein_count") %>%
  mutate(percentage = protein_count / sum(protein_count) * 100)

tiff(file.path(output_dir, "cross_database_consistency.tiff"), 
     width = 12, height = 8, units = "in", res = 600, compression = "lzw")

p1 <- ggplot(consistency_summary, aes(x = factor(databases_count), y = protein_count)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Cross-Database Protein Consistency",
       x = "Number of Databases", y = "Number of Proteins") +
  geom_text(aes(label = paste0(protein_count, "\n(", round(percentage, 1), "%)")), 
            vjust = 1.5, color = "white", size = 3.5)

# Top proteins found in most databases
top_consistent <- protein_db_counts %>%
  filter(databases_count >= 3) %>%
  head(15)

p2 <- ggplot(top_consistent, aes(x = reorder(gene_symbol, databases_count), 
                                y = databases_count)) +
  geom_col(fill = "darkgreen", alpha = 0.8) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Most Consistently Detected Proteins",
       x = "Protein ID", y = "Number of Databases") +
  geom_text(aes(label = databases_count), hjust = -0.1, size = 3)

grid.arrange(p1, p2, ncol = 2,
             top = "Cross-Database Protein Detection Consistency")

dev.off()

# Save analysis results
write.csv(efficiency_data, file.path(output_dir, "database_efficiency_analysis.csv"), 
          row.names = FALSE)

write.csv(protein_db_counts, file.path(output_dir, "cross_database_protein_consistency.csv"), 
          row.names = FALSE)

write.csv(comparison_data, file.path(output_dir, "immune_vs_nonimmune_summary.csv"), 
          row.names = FALSE)

if(exists("tcell_comparison")) {
  write.csv(tcell_comparison, file.path(output_dir, "tcell_subset_comparison.csv"), 
            row.names = FALSE)
}

# Summary report
sink(file.path(output_dir, "comparative_analysis_summary.txt"))
cat("Comparative Analysis Summary\n")
cat("===========================\n\n")

cat("Dataset Overview:\n")
cat(paste("Total unique proteins:", length(unique(protein_data$gene_symbol)), "\n"))
cat(paste("Total databases:", length(unique(protein_data$database)), "\n"))
cat(paste("Total cell types:", length(unique(protein_data$cell_type)), "\n"))
cat(paste("Total measurements:", nrow(protein_data), "\n\n"))

cat("Immune vs Non-Immune Comparison:\n")
cat(paste("Immune-specific proteins:", comparison_data$Count[1], 
          "(", round(comparison_data$Percentage[1], 1), "%)\n"))
cat(paste("Shared proteins:", comparison_data$Count[2], 
          "(", round(comparison_data$Percentage[2], 1), "%)\n"))
cat(paste("Non-immune specific proteins:", comparison_data$Count[3], 
          "(", round(comparison_data$Percentage[3], 1), "%)\n\n"))

if(exists("tcell_comparison")) {
  cat("T Cell Subset Comparison (CD4 vs CD8):\n")
  cat(paste("CD4-specific proteins:", tcell_comparison$Count[1], 
            "(", round(tcell_comparison$Percentage[1], 1), "%)\n"))
  cat(paste("Shared proteins:", tcell_comparison$Count[2], 
            "(", round(tcell_comparison$Percentage[2], 1), "%)\n"))
  cat(paste("CD8-specific proteins:", tcell_comparison$Count[3], 
            "(", round(tcell_comparison$Percentage[3], 1), "%)\n\n"))
}

cat("Database Efficiency Rankings:\n")
for(i in 1:nrow(efficiency_data)) {
  db <- efficiency_data$database[i]
  proteins <- efficiency_data$unique_proteins[i]
  coverage <- efficiency_data$coverage_breadth[i]
  cat(paste("  ", i, ".", db, "- ", proteins, "unique proteins across", 
            coverage, "cell types\n"))
}

cat("\nTop 10 Most Consistently Detected Proteins (across databases):\n")
top_consistent_summary <- head(protein_db_counts, 10)
for(i in 1:nrow(top_consistent_summary)) {
  protein <- top_consistent_summary$gene_symbol[i]
  db_count <- top_consistent_summary$databases_count[i]
  databases <- top_consistent_summary$databases[i]
  cat(paste("  ", i, ".", protein, "- found in", db_count, "databases (", databases, ")\n"))
}

sink()

message("Comparative analysis completed!")
message(paste("Results saved to:", output_dir)) 
