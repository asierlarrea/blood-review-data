#!/usr/bin/env Rscript

# Figure 6: Database Correlation and Coverage Analysis
# Working directly with raw database files for protein-level analysis

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
required_packages <- c("ggplot2", "dplyr", "tidyr", "corrplot", "VennDiagram", 
                      "RColorBrewer", "scales", "gridExtra", "stringr", "biomaRt")
install_if_missing(required_packages)

# Load enhanced ID mapping functions with biomaRt ENSP mappings
source("enhanced_fast_mapping.R")

# Create output directory
output_dir <- "plots/01_Database_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to extract actual gene symbols from descriptions
extract_gene_symbol <- function(description) {
  if(is.na(description) || description == "") {
    return(NA)
  }
  
  desc <- as.character(description)
  
  # Pattern 1: GN=GENE_NAME (UniProt format)
  if(grepl("GN=([A-Za-z0-9_-]+)", desc)) {
    gene <- gsub(".*GN=([A-Za-z0-9_-]+).*", "\\1", desc)
    return(toupper(gene))
  }
  
  # Pattern 2: \GName=GENE_NAME (PeptideAtlas format)
  if(grepl("\\\\GName=([A-Za-z0-9_-]+)", desc)) {
    gene <- gsub(".*\\\\GName=([A-Za-z0-9_-]+).*", "\\1", desc)
    return(toupper(gene))
  }
  
  # Pattern 3: Gene name at the beginning of description
  if(grepl("^([A-Z][A-Z0-9_-]+)\\s", desc)) {
    gene <- gsub("^([A-Z][A-Z0-9_-]+)\\s.*", "\\1", desc)
    if(nchar(gene) >= 2 && nchar(gene) <= 15) {  # Reasonable gene name length
      return(toupper(gene))
    }
  }
  
  # Pattern 4: Look for gene symbols in parentheses or after specific keywords
  if(grepl("gene[[:space:]]*[=:]?[[:space:]]*([A-Z][A-Z0-9_-]+)", desc, ignore.case = TRUE)) {
    gene <- gsub(".*gene[[:space:]]*[=:]?[[:space:]]*([A-Z][A-Z0-9_-]+).*", "\\1", desc, ignore.case = TRUE)
    return(toupper(gene))
  }
  
  return(NA)
}

# Function to parse cell type from column name
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

# Function to read and process database files - extracting real gene symbols
read_database_files <- function() {
  message("Reading and processing database files to extract gene symbols...")
  
  # Initialize data structure
  all_plasma_data <- list()
  all_celltype_data <- list()
  
  # 1. PeptideAtlas - PLASMA
  message("  - Processing PeptideAtlas (Plasma)...")
  if(file.exists("PeptideAtlas.csv")) {
    peptideatlas <- read.csv("PeptideAtlas.csv", stringsAsFactors = FALSE)
    
    gene_symbols <- character(nrow(peptideatlas))
    concentrations <- numeric(nrow(peptideatlas))
    valid_entries <- logical(nrow(peptideatlas))
    
    for(i in 1:nrow(peptideatlas)) {
      gene_symbol <- NA
      
      # First try biosequence_gene_name column
      if(!is.na(peptideatlas$biosequence_gene_name[i]) && peptideatlas$biosequence_gene_name[i] != "") {
        gene_name <- as.character(peptideatlas$biosequence_gene_name[i])
        if(!grepl("[;,]", gene_name)) {  # Single gene only
          gene_symbol <- toupper(trimws(gene_name))
        }
      }
      
      # If not found, try extracting from description
      if(is.na(gene_symbol)) {
        gene_symbol <- extract_gene_symbol(peptideatlas$biosequence_desc[i])
      }
      
      # Only keep valid gene symbols (real HGNC symbols)
      if(!is.na(gene_symbol) && 
         !grepl("^(UNIPROT_|GENE_|UNKNOWN_|ENSP_)", gene_symbol) &&
         nchar(gene_symbol) >= 2 && nchar(gene_symbol) <= 15) {
        
        # Try to get normalized PSMs per 100K value
        conc_value <- suppressWarnings(as.numeric(peptideatlas$norm_PSMs_per_100K[i]))
        if(length(conc_value) > 0 && !is.na(conc_value) && conc_value > 0) {
          gene_symbols[i] <- gene_symbol
          concentrations[i] <- conc_value
          valid_entries[i] <- TRUE
        }
      }
    }
    
    if(sum(valid_entries) > 0) {
      pa_data <- data.frame(
        Gene = gene_symbols[valid_entries],
        Concentration = concentrations[valid_entries],
        Database = "PeptideAtlas",
        Category = "Plasma",
        CellType = "Plasma",
        Dataset = "PeptideAtlas",
        stringsAsFactors = FALSE
      )
      all_plasma_data <- c(all_plasma_data, list(pa_data))
      message(paste("    PeptideAtlas: Kept", sum(valid_entries), "of", nrow(peptideatlas), "entries"))
    }
  }
  
  # 2. GPMDB Plasma
  message("  - Processing GPMDB Plasma...")
  if(file.exists("GPMDB_Plasma.csv")) {
    gpmdb_plasma <- read.csv("GPMDB_Plasma.csv", stringsAsFactors = FALSE)
    
    # Find header line
    header_line <- which(grepl("^#,accession,total", readLines("GPMDB_Plasma.csv")))
    if(length(header_line) > 0) {
      # Skip to the line after the header
      gpmdb_plasma <- read.csv("GPMDB_Plasma.csv", skip = header_line, stringsAsFactors = FALSE, header = FALSE)
      
      # Set column names manually
      colnames(gpmdb_plasma) <- c("rank", "accession", "total", "log_e", "EC", "description")
      
      # Filter valid rows
      valid_rows <- !is.na(gpmdb_plasma$accession) & gpmdb_plasma$accession != "" & 
                   !is.na(gpmdb_plasma$total) & gpmdb_plasma$total > 0
      
      if(sum(valid_rows) > 0) {
        # Use proper biological ID mapping for GPMDB accessions
        accessions <- gpmdb_plasma$accession[valid_rows]
        descriptions <- gpmdb_plasma$description[valid_rows]
        
        # Map accessions to gene symbols using enhanced mapping with biomaRt
        gene_symbols <- enhanced_fast_map_to_gene_symbol(accessions, descriptions)
        
        # Remove entries without valid gene symbols
        valid_genes <- !is.na(gene_symbols) & gene_symbols != "" & 
                      nchar(gene_symbols) >= 2 & nchar(gene_symbols) <= 15 &
                      !grepl("^(ENS|UNIPROT_|GENE_|UNKNOWN_)", gene_symbols)
        
        if(sum(valid_genes) > 0) {
          gpmdb_data <- data.frame(
            Gene = gene_symbols[valid_genes],
            Concentration = gpmdb_plasma$total[valid_rows][valid_genes],
            Database = "GPMDB",
            Category = "Plasma",
            CellType = "Plasma",
            Dataset = "GPMDB_Plasma",
            stringsAsFactors = FALSE
          )
          all_plasma_data <- c(all_plasma_data, list(gpmdb_data))
          message(paste("    GPMDB Plasma: Kept", sum(valid_genes), "entries"))
        }
      }
    }
  }
  
  # 3. PaxDB Plasma
  message("  - Processing PaxDB Plasma...")
  paxdb_files <- c("PaxDB_plasma.csv", "PaxDB_Plasma.csv")
  for(file in paxdb_files) {
    if(file.exists(file)) {
      message(paste("    Processing", file))
      
      paxdb_data <- read.csv(file, stringsAsFactors = FALSE)
      
      if("string_external_id" %in% colnames(paxdb_data) && "abundance" %in% colnames(paxdb_data)) {
        # Extract ENSP IDs and abundances
        valid_rows <- !is.na(paxdb_data$string_external_id) & 
                     !is.na(paxdb_data$abundance) & 
                     paxdb_data$abundance > 0 &
                     grepl("9606\\.ENSP[0-9]+", paxdb_data$string_external_id)
        
        if(sum(valid_rows) > 0) {
          # Extract ENSP IDs
          ensp_ids <- gsub("9606\\.", "", paxdb_data$string_external_id[valid_rows])
          abundances <- paxdb_data$abundance[valid_rows]
          
          # Map ENSP IDs to gene symbols using enhanced mapping with biomaRt
          gene_symbols <- enhanced_fast_map_to_gene_symbol(ensp_ids)
          
          # Filter for valid gene mappings
          valid_mappings <- !is.na(gene_symbols) & gene_symbols != "" &
                           nchar(gene_symbols) >= 2 & nchar(gene_symbols) <= 15
          
          if(sum(valid_mappings) > 0) {
            paxdb_plasma_data <- data.frame(
              Gene = gene_symbols[valid_mappings],
              Concentration = abundances[valid_mappings],
              Database = "PaxDB",
              Category = "Plasma",
              CellType = "Plasma",
              Dataset = "PaxDB_Plasma",
              stringsAsFactors = FALSE
            )
            all_plasma_data <- c(all_plasma_data, list(paxdb_plasma_data))
            message(paste("    PaxDB Plasma: Kept", sum(valid_mappings), "entries from", file))
          }
        }
      }
      break  # Only process one of the PaxDB files
    }
  }
  
  # 4. HPA files - PLASMA
  message("  - Processing HPA files (Plasma)...")
  hpa_files <- list.files(pattern = "^HPA_.*\\.csv$")
  for(file in hpa_files) {
    if(file.exists(file)) {
      message(paste("    Processing", file))
      
      hpa_data <- read.csv(file, stringsAsFactors = FALSE)
      
      # Find concentration column
      concentration_col <- NULL
      if("Mean.concentration" %in% colnames(hpa_data)) {
        concentration_col <- "Mean.concentration"
      } else if("Concentration" %in% colnames(hpa_data)) {
        concentration_col <- "Concentration"
      }
      
      if(!is.null(concentration_col) && "Gene" %in% colnames(hpa_data)) {
        valid_genes <- character()
        valid_concentrations <- numeric()
        
        for(i in 1:nrow(hpa_data)) {
          gene_name <- as.character(hpa_data$Gene[i])
          conc_str <- as.character(hpa_data[[concentration_col]][i])
          
          # Check for single gene and valid concentration
          if(!is.na(gene_name) && gene_name != "" && !grepl("[;,]", gene_name) &&
             !is.na(conc_str) && conc_str != "") {
            
            # Clean concentration string and extract numeric
            tryCatch({
              clean_conc_str <- iconv(conc_str, to = "UTF-8", sub = "")
              conc_numeric <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", clean_conc_str)))
              
              if(!is.na(conc_numeric) && conc_numeric > 0) {
                valid_genes <- c(valid_genes, toupper(trimws(gene_name)))
                valid_concentrations <- c(valid_concentrations, conc_numeric)
              }
            }, error = function(e) {
              # Skip entries with encoding issues
              message(paste("    Skipping entry", i, "due to encoding issue"))
            })
          }
        }
        
        if(length(valid_genes) > 0) {
          # Extract HPA technique from filename (e.g., HPA_MS.csv -> HPA_MS)
          hpa_technique <- gsub("\\.csv$", "", file)
          
          hpa_data <- data.frame(
            Gene = valid_genes,
            Concentration = valid_concentrations,
            Database = "HPA",
            Category = "Plasma",
            CellType = "Plasma",
            Dataset = hpa_technique,
            stringsAsFactors = FALSE
          )
          all_plasma_data <- c(all_plasma_data, list(hpa_data))
          message(paste("    HPA: Kept", length(valid_genes), "entries from", file))
        }
      }
    }
  }
  
  # 5. ProteomeXchange files - CELL TYPES
  message("  - Processing ProteomeXchange files (Cell Types)...")
  pxd_files <- list.files(pattern = "^PXD.*\\.csv$")
  for(file in pxd_files) {
    if(file.exists(file)) {
      message(paste("    Processing", file, "for cell types"))
      
      # Read data with header
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
                  # Extract dataset accession from filename (e.g., PXD004352.csv -> PXD004352, PXD040957_CD8.csv -> PXD040957)
                  dataset_accession <- gsub("\\.csv$", "", file)
                  dataset_accession <- gsub("_CD8$|_Macrophages$", "", dataset_accession)
                  
                  px_data <- data.frame(
                    Gene = valid_genes,
                    Concentration = valid_intensities,
                    Database = "ProteomeXchange",
                    Category = "CellType",
                    CellType = cell_type,
                    Dataset = dataset_accession,
                    stringsAsFactors = FALSE
                  )
                  all_celltype_data <- c(all_celltype_data, list(px_data))
                  message(paste("    ProteomeXchange:", cell_type, "- Kept", length(valid_genes), "entries from", file))
                }
              }
            }
          }
        }
      }
    }
  }
  
  # Combine and process data
  message("  - Final processing...")
  
  final_data <- data.frame()
  
  # Process plasma data
  if(length(all_plasma_data) > 0) {
    plasma_data <- do.call(rbind, all_plasma_data)
    plasma_data <- as.data.frame(plasma_data, stringsAsFactors = FALSE)
    plasma_data$Concentration <- as.numeric(plasma_data$Concentration)
    
    # Remove duplicates: if same gene appears multiple times in same database/dataset, keep highest concentration
    plasma_data <- plasma_data %>%
      group_by(Gene, Database, Dataset, Category, CellType) %>%
      slice_max(Concentration, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    final_data <- rbind(final_data, plasma_data)
  }
  
  # Process cell type data
  if(length(all_celltype_data) > 0) {
    celltype_data <- do.call(rbind, all_celltype_data)
    celltype_data <- as.data.frame(celltype_data, stringsAsFactors = FALSE)
    celltype_data$Concentration <- as.numeric(celltype_data$Concentration)
    
    # Remove duplicates: if same gene appears multiple times in same database/dataset/celltype, keep highest concentration
    celltype_data <- celltype_data %>%
      group_by(Gene, Database, Dataset, Category, CellType) %>%
      slice_max(Concentration, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    final_data <- rbind(final_data, celltype_data)
  }
  
  if(nrow(final_data) > 0) {
    message(paste("Total unique gene measurements:", nrow(final_data)))
    message(paste("Unique genes across all sources:", length(unique(final_data$Gene))))
    
    # Summary by category and database
    summary_stats <- final_data %>%
      group_by(Category, Database, CellType) %>%
      summarise(
        gene_count = n_distinct(Gene),
        measurement_count = n(),
        .groups = 'drop'
      )
    
    message("\nSummary by category:")
    for(i in 1:nrow(summary_stats)) {
      message(paste("  ", summary_stats$Category[i], "-", summary_stats$Database[i], 
                    "(", summary_stats$CellType[i], "):", 
                    summary_stats$gene_count[i], "unique genes,", 
                    summary_stats$measurement_count[i], "measurements"))
    }
    
    return(final_data)
  } else {
    message("No valid data found!")
    return(data.frame())
  }
}

# Function to analyze gene overlaps by category
analyze_gene_overlaps_by_category <- function(gene_data, category_filter = NULL) {
  if(!is.null(category_filter)) {
    gene_data <- gene_data %>% filter(Category == category_filter)
    message(paste("Analyzing gene overlaps for category:", category_filter))
  } else {
    message("Analyzing gene overlaps across all categories...")
  }
  
  # Get genes per database-celltype combination
  db_genes <- gene_data %>%
    mutate(DB_CellType = paste0(Database, " (", CellType, ")")) %>%
    group_by(DB_CellType) %>%
    summarise(genes = list(unique(Gene)), .groups = 'drop')
  
  databases <- db_genes$DB_CellType
  n_db <- length(databases)
  
  if(n_db < 2) {
    message("Not enough data sources for overlap analysis")
    return(NULL)
  }
  
  # Calculate pairwise overlaps
  overlap_matrix <- matrix(0, nrow = n_db, ncol = n_db)
  rownames(overlap_matrix) <- databases
  colnames(overlap_matrix) <- databases
  
  for(i in 1:n_db) {
    for(j in 1:n_db) {
      if(i == j) {
        overlap_matrix[i, j] <- length(db_genes$genes[[i]])
      } else {
        overlap_matrix[i, j] <- length(intersect(db_genes$genes[[i]], 
                                                db_genes$genes[[j]]))
      }
    }
  }
  
  # Calculate correlation matrix (Jaccard similarity)
  correlation_matrix <- matrix(0, nrow = n_db, ncol = n_db)
  rownames(correlation_matrix) <- databases
  colnames(correlation_matrix) <- databases
  
  for(i in 1:n_db) {
    for(j in 1:n_db) {
      if(i == j) {
        correlation_matrix[i, j] <- 1.0
      } else {
        intersection_size <- overlap_matrix[i, j]
        union_size <- length(db_genes$genes[[i]]) + length(db_genes$genes[[j]]) - intersection_size
        if(union_size > 0) {
          correlation_matrix[i, j] <- intersection_size / union_size
        } else {
          correlation_matrix[i, j] <- 0
        }
      }
    }
  }
  
  # Gene distribution across data sources
  gene_counts <- gene_data %>%
    mutate(DB_CellType = paste0(Database, " (", CellType, ")")) %>%
    count(DB_CellType, name = "gene_count")
  
  genes_per_source_count <- gene_data %>%
    mutate(DB_CellType = paste0(Database, " (", CellType, ")")) %>%
    group_by(Gene) %>%
    summarise(source_count = n_distinct(DB_CellType), .groups = 'drop') %>%
    count(source_count) %>%
    arrange(source_count)
  
  message("Gene distribution across data sources:")
  for(i in 1:nrow(genes_per_source_count)) {
    message(paste("  ", genes_per_source_count$n[i], "genes found in", genes_per_source_count$source_count[i], "source(s)"))
  }
  
  return(list(
    overlap_matrix = overlap_matrix,
    correlation_matrix = correlation_matrix,
    genes_per_source = gene_counts,
    gene_distribution = genes_per_source_count,
    category = ifelse(is.null(category_filter), "All", category_filter)
  ))
}

# Read all database files and extract gene symbols
gene_data <- read_database_files()

if(nrow(gene_data) == 0) {
  stop("No valid gene data found. Please check your data files.")
}

# Analyze overlaps by category
plasma_analysis <- analyze_gene_overlaps_by_category(gene_data, "Plasma")
celltype_analysis <- analyze_gene_overlaps_by_category(gene_data, "CellType")
overall_analysis <- analyze_gene_overlaps_by_category(gene_data, NULL)

# Plot functions
create_correlation_plot <- function(analysis_result, category, filename) {
  if(is.null(analysis_result)) {
    message(paste("Skipping correlation plot for", category, "- insufficient data"))
    return()
  }
  
  message(paste("Creating gene overlap correlation matrix for", category, "..."))
  
  tryCatch({
    tiff(file.path(output_dir, filename), 
         width = 10, height = 8, units = "in", res = 600)
  }, error = function(e) {
    message(paste("Error creating TIFF file:", e$message))
    png(file.path(output_dir, gsub("\\.tiff$", ".png", filename)), 
        width = 10, height = 8, units = "in", res = 300)
  })
  
  # Convert correlation matrix to long format for ggplot (without reshape2)
  correlation_matrix <- analysis_result$correlation_matrix
  databases <- rownames(correlation_matrix)
  correlation_long <- data.frame()
  for(i in 1:nrow(correlation_matrix)) {
    for(j in 1:ncol(correlation_matrix)) {
      correlation_long <- rbind(correlation_long, data.frame(
        Database1 = databases[i],
        Database2 = databases[j],
        Correlation = correlation_matrix[i, j]
      ))
    }
  }
  
  p <- ggplot(correlation_long, aes(x = Database1, y = Database2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue", 
                         midpoint = 0.5, name = "Jaccard\nSimilarity", 
                         limits = c(0, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          panel.grid = element_blank(),
          plot.title = element_text(size = 14, face = "bold")) +
    labs(title = paste("Gene Overlap Between Data Sources -", category),
         subtitle = "Jaccard similarity coefficient of gene sets",
         x = "Data Source", y = "Data Source") +
    geom_text(aes(label = sprintf("%.3f", Correlation)), 
              color = ifelse(correlation_long$Correlation > 0.5, "white", "black"), 
              size = 3)
  
  print(p)
  dev.off()
}

create_distribution_plot <- function(analysis_result, category, filename) {
  if(is.null(analysis_result)) {
    message(paste("Skipping distribution plot for", category, "- insufficient data"))
    return()
  }
  
  message(paste("Creating gene distribution plot for", category, "..."))
  
  tiff(file.path(output_dir, filename), 
       width = 10, height = 6, units = "in", res = 600)
  
  p <- ggplot(analysis_result$gene_distribution, aes(x = factor(source_count), y = n)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    theme_minimal() +
    labs(title = paste("Gene Distribution Across Data Sources -", category),
         subtitle = "Number of genes found in X data sources",
         x = "Number of Data Sources", y = "Number of Genes") +
    geom_text(aes(label = n), vjust = -0.5, size = 4) +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  print(p)
  dev.off()
}

create_source_counts_plot <- function(analysis_result, category, filename) {
  if(is.null(analysis_result)) {
    message(paste("Skipping source counts plot for", category, "- insufficient data"))
    return()
  }
  
  message(paste("Creating data source counts plot for", category, "..."))
  
  tiff(file.path(output_dir, filename), 
       width = 12, height = 6, units = "in", res = 600)
  
  p <- ggplot(analysis_result$genes_per_source, aes(x = reorder(DB_CellType, gene_count), y = gene_count)) +
    geom_col(fill = "darkgreen", alpha = 0.8) +
    coord_flip() +
    theme_minimal() +
    labs(title = paste("Number of Genes per Data Source -", category),
         subtitle = "Total unique genes with valid HGNC symbols",
         x = "Data Source", y = "Number of Genes") +
    geom_text(aes(label = gene_count), hjust = -0.1, size = 4) +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  print(p)
  dev.off()
}

# Function to create stacked bar plot by dataset accession
create_stacked_dataset_plot <- function(gene_data, category_filter = NULL, filename) {
  if(!is.null(category_filter)) {
    plot_data <- gene_data %>% filter(Category == category_filter)
    plot_title <- paste("Gene Counts by Cell Type and Dataset -", category_filter)
  } else {
    plot_data <- gene_data
    plot_title <- "Gene Counts by Cell Type and Dataset - All Categories"
  }
  
  if(nrow(plot_data) == 0) {
    message(paste("Skipping stacked plot for", ifelse(is.null(category_filter), "All", category_filter), "- no data"))
    return()
  }
  
  # Calculate gene counts by CellType and Dataset
  stacked_data <- plot_data %>%
    group_by(CellType, Dataset) %>%
    summarise(gene_count = n_distinct(Gene), .groups = 'drop') %>%
    arrange(desc(gene_count))
  
  # Create custom color palette for datasets
  n_datasets <- length(unique(stacked_data$Dataset))
  if(n_datasets <= 12) {
    dataset_colors <- RColorBrewer::brewer.pal(min(n_datasets, 12), "Set3")
  } else {
    dataset_colors <- rainbow(n_datasets)
  }
  
  message(paste("Creating stacked dataset plot for", ifelse(is.null(category_filter), "All categories", category_filter), "..."))
  
  tiff(file.path(output_dir, filename), 
       width = 14, height = 8, units = "in", res = 600)
  
  p <- ggplot(stacked_data, aes(x = reorder(CellType, gene_count, sum), y = gene_count, fill = Dataset)) +
    geom_col(position = "stack", alpha = 0.8, color = "white", linewidth = 0.2) +
    coord_flip() +
    scale_fill_manual(values = dataset_colors, name = "Dataset\nAccession") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = plot_title,
      subtitle = "Stacked bars show contribution of each dataset to cell type protein counts",
      x = "Cell Type", 
      y = "Number of Unique Genes"
    ) +
    geom_text(
      data = stacked_data %>% 
        group_by(CellType) %>% 
        summarise(total_genes = sum(gene_count), .groups = 'drop'),
      aes(x = CellType, y = total_genes, label = total_genes, fill = NULL),
      hjust = -0.1, size = 4, fontface = "bold"
    )
  
  print(p)
  dev.off()
}

# Create plots for each category
create_correlation_plot(plasma_analysis, "Plasma", "plasma_gene_overlap_correlation.tiff")
create_distribution_plot(plasma_analysis, "Plasma", "plasma_gene_distribution.tiff")
create_source_counts_plot(plasma_analysis, "Plasma", "plasma_source_counts.tiff")

create_correlation_plot(celltype_analysis, "Cell Types", "celltype_gene_overlap_correlation.tiff")
create_distribution_plot(celltype_analysis, "Cell Types", "celltype_gene_distribution.tiff")
create_source_counts_plot(celltype_analysis, "Cell Types", "celltype_source_counts.tiff")
create_stacked_dataset_plot(gene_data, "CellType", "celltype_stacked_by_dataset.tiff")

create_correlation_plot(overall_analysis, "Overall", "overall_gene_overlap_correlation.tiff")
create_distribution_plot(overall_analysis, "Overall", "overall_gene_distribution.tiff")
create_source_counts_plot(overall_analysis, "Overall", "overall_source_counts.tiff")

# Save all data and results to CSV files
message("Saving data and results to CSV files...")

# 1. Save the final processed gene data
write.csv(gene_data, file.path(output_dir, "processed_gene_data_with_categories.csv"), row.names = FALSE)

# 2. Save plasma-specific data
plasma_data <- gene_data %>% filter(Category == "Plasma")
if(nrow(plasma_data) > 0) {
  write.csv(plasma_data, file.path(output_dir, "plasma_gene_data.csv"), row.names = FALSE)
}

# 3. Save cell type-specific data
celltype_data <- gene_data %>% filter(Category == "CellType")
if(nrow(celltype_data) > 0) {
  write.csv(celltype_data, file.path(output_dir, "celltype_gene_data.csv"), row.names = FALSE)
}

# Function to save analysis results
save_analysis_results <- function(analysis_result, prefix) {
  if(is.null(analysis_result)) return()
  
  # Save overlap matrix
  write.csv(analysis_result$overlap_matrix, file.path(output_dir, paste0(prefix, "_gene_overlap_matrix.csv")))
  
  # Save correlation matrix
  write.csv(analysis_result$correlation_matrix, file.path(output_dir, paste0(prefix, "_jaccard_similarity_matrix.csv")))
  
  # Save gene distribution
  write.csv(analysis_result$gene_distribution, file.path(output_dir, paste0(prefix, "_gene_distribution.csv")), row.names = FALSE)
  
  # Save source counts
  write.csv(analysis_result$genes_per_source, file.path(output_dir, paste0(prefix, "_source_counts.csv")), row.names = FALSE)
  
  # Create pairwise results if applicable
  if(!is.null(analysis_result$overlap_matrix) && nrow(analysis_result$overlap_matrix) > 1) {
    databases <- rownames(analysis_result$overlap_matrix)
    pairwise_results <- data.frame()
    for(i in 1:(length(databases)-1)) {
      for(j in (i+1):length(databases)) {
        pairwise_results <- rbind(pairwise_results, data.frame(
          Source1 = databases[i],
          Source2 = databases[j],
          Shared_Genes = analysis_result$overlap_matrix[i, j],
          Source1_Total = analysis_result$overlap_matrix[i, i],
          Source2_Total = analysis_result$overlap_matrix[j, j],
          Jaccard_Similarity = analysis_result$correlation_matrix[i, j],
          stringsAsFactors = FALSE
        ))
      }
    }
    write.csv(pairwise_results, file.path(output_dir, paste0(prefix, "_pairwise_overlaps.csv")), row.names = FALSE)
  }
}

# Save analysis results for each category
save_analysis_results(plasma_analysis, "plasma")
save_analysis_results(celltype_analysis, "celltype")
save_analysis_results(overall_analysis, "overall")

# Create detailed gene overlap information by category
create_gene_details <- function(data_subset, filename) {
  if(nrow(data_subset) == 0) return()
  
  gene_overlap_detail <- data_subset %>%
    mutate(Source = paste0(Database, " (", CellType, ")")) %>%
    group_by(Gene) %>%
    summarise(
      Sources = paste(unique(Source), collapse = ";"),
      Source_Count = n_distinct(Source),
      CellTypes = paste(unique(CellType), collapse = ";"),
      Total_Measurements = n(),
      Avg_Concentration = mean(Concentration, na.rm = TRUE),
      Max_Concentration = max(Concentration, na.rm = TRUE),
      Min_Concentration = min(Concentration, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(Source_Count), Gene)
  
  write.csv(gene_overlap_detail, file.path(output_dir, filename), row.names = FALSE)
}

create_gene_details(plasma_data, "plasma_gene_details.csv")
create_gene_details(celltype_data, "celltype_gene_details.csv")
create_gene_details(gene_data, "overall_gene_details.csv")

# Write comprehensive summary
sink(file.path(output_dir, "comprehensive_gene_analysis_summary.txt"))
cat("Comprehensive Gene-Based Database Analysis Summary\n")
cat("================================================\n\n")

cat("OVERALL SUMMARY:\n")
cat(paste("Total unique genes across all sources:", length(unique(gene_data$Gene)), "\n"))
cat(paste("Total measurements:", nrow(gene_data), "\n\n"))

cat("BY CATEGORY:\n")
if(nrow(plasma_data) > 0) {
  cat("PLASMA:\n")
  cat(paste("  Unique genes:", length(unique(plasma_data$Gene)), "\n"))
  cat(paste("  Measurements:", nrow(plasma_data), "\n"))
  plasma_summary <- plasma_data %>% 
    group_by(Database, CellType) %>% 
    summarise(genes = n_distinct(Gene), .groups = 'drop')
  for(i in 1:nrow(plasma_summary)) {
    cat(paste("  ", plasma_summary$Database[i], ":", plasma_summary$genes[i], "genes\n"))
  }
  cat("\n")
}

if(nrow(celltype_data) > 0) {
  cat("CELL TYPES:\n")
  cat(paste("  Unique genes:", length(unique(celltype_data$Gene)), "\n"))
  cat(paste("  Measurements:", nrow(celltype_data), "\n"))
  celltype_summary <- celltype_data %>% 
    group_by(Database, CellType) %>% 
    summarise(genes = n_distinct(Gene), .groups = 'drop')
  for(i in 1:nrow(celltype_summary)) {
    cat(paste("  ", celltype_summary$Database[i], "-", celltype_summary$CellType[i], ":", celltype_summary$genes[i], "genes\n"))
  }
  cat("\n")
}

if(!is.null(plasma_analysis) && !is.null(plasma_analysis$gene_distribution)) {
  cat("PLASMA GENE DISTRIBUTION:\n")
  for(i in 1:nrow(plasma_analysis$gene_distribution)) {
    cat(paste("  ", plasma_analysis$gene_distribution$n[i], "genes found in", 
              plasma_analysis$gene_distribution$source_count[i], "plasma source(s)\n"))
  }
  cat("\n")
}

if(!is.null(celltype_analysis) && !is.null(celltype_analysis$gene_distribution)) {
  cat("CELL TYPE GENE DISTRIBUTION:\n")
  for(i in 1:nrow(celltype_analysis$gene_distribution)) {
    cat(paste("  ", celltype_analysis$gene_distribution$n[i], "genes found in", 
              celltype_analysis$gene_distribution$source_count[i], "cell type source(s)\n"))
  }
  cat("\n")
}

sink()

message("Comprehensive gene-based database analysis completed!")
message(paste("Results saved to:", output_dir))
message("\nSummary:")
message(paste("- Total unique genes:", length(unique(gene_data$Gene))))
message(paste("- Plasma genes:", if(nrow(plasma_data) > 0) length(unique(plasma_data$Gene)) else 0))
message(paste("- Cell type genes:", if(nrow(celltype_data) > 0) length(unique(celltype_data$Gene)) else 0))

message("\nCSV files generated:")
message("  MAIN DATA:")
message("    - processed_gene_data_with_categories.csv - All data with categories")
message("    - plasma_gene_data.csv - Plasma-specific data")
message("    - celltype_gene_data.csv - Cell type-specific data")
message("  ANALYSIS RESULTS:")
message("    - plasma_* - Plasma analysis results")
message("    - celltype_* - Cell type analysis results") 
message("    - overall_* - Combined analysis results")
message("  PLOTS:")
message("    - *_correlation.tiff - Overlap correlation matrices")
message("    - *_distribution.tiff - Gene distribution plots")
message("    - *_counts.tiff - Source count plots") 