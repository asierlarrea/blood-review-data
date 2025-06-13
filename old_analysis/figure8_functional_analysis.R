#!/usr/bin/env Rscript

# Load utilities and set up paths
source("scripts/utilities/load_packages.R")
ensure_output_dirs()


# Figure 8: Functional Analysis of Blood Proteins
# Working directly with raw database files for functional categorization

# Set CRAN mirror
if(is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "pheatmap", "RColorBrewer", 
                      "viridis", "scales", "gridExtra", "stringr")
load_packages(required_packages)



# Create output directory
output_dir <- "outputs/plots/05_Functional_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to map protein IDs to gene symbols based on database type
map_to_gene_symbol <- function(protein_ids, descriptions = NULL, database_type = "general") {
  gene_symbols <- character(length(protein_ids))
  
  for(i in seq_along(protein_ids)) {
    protein_id <- as.character(protein_ids[i])
    description <- if(!is.null(descriptions)) as.character(descriptions[i]) else ""
    
    # Initialize with fallback
    gene_symbol <- protein_id
    
    # Database-specific handling
    if(database_type == "paxdb") {
      # PaxDB format: 9606.ENSP00000357727 - remove organism prefix
      clean_id <- gsub("^9606\\.", "", protein_id)
      gene_symbol <- paste0("GENE_", substr(clean_id, 1, 8))
    } else if(database_type == "gpmdb") {
      # GPMDB format: ENSP00000426179 - Ensembl protein IDs
      gene_symbol <- paste0("GENE_", substr(protein_id, 1, 8))
    } else if(database_type == "peptideatlas") {
      # PeptideAtlas format: P27144 - UniProt accessions
      # Try to extract from description first
      if(!is.null(descriptions) && !is.na(description) && description != "") {
        # Look for gene symbols in various formats in descriptions
        if(grepl("GN=([A-Z0-9_-]+)", description)) {
          gene_symbol <- gsub(".*GN=([A-Z0-9_-]+).*", "\\1", description)
        } else if(grepl("^([A-Z0-9_-]+)\\s", description)) {
          gene_symbol <- gsub("^([A-Z0-9_-]+)\\s.*", "\\1", description)
        } else {
          gene_symbol <- paste0("UNIPROT_", substr(protein_id, 1, 6))
        }
      } else {
        gene_symbol <- paste0("UNIPROT_", substr(protein_id, 1, 6))
      }
    } else if(database_type == "hpa") {
      # HPA already has gene symbols - just clean them
      gene_symbol <- toupper(protein_id)
    } else {
      # General case - try to extract from description
      if(!is.null(descriptions) && !is.na(description) && description != "") {
        if(grepl("GN=([A-Z0-9_-]+)", description)) {
          gene_symbol <- gsub(".*GN=([A-Z0-9_-]+).*", "\\1", description)
        } else if(grepl("^([A-Z0-9_-]+)[,\\s]", description)) {
          gene_symbol <- gsub("^([A-Z0-9_-]+)[,\\s].*", "\\1", description)
        } else {
          gene_symbol <- paste0("GENE_", substr(protein_id, 1, 8))
        }
      } else {
        gene_symbol <- paste0("GENE_", substr(protein_id, 1, 8))
      }
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

# Function to read and extract protein data from all databases
read_protein_data <- function() {
  message("Reading protein data from raw database files...")
  
  protein_data <- list()
  
  # 1. PeptideAtlas
  message("  - Reading PeptideAtlas...")
  peptideatlas <- read.csv(get_data_path("peptideatlas.csv"), stringsAsFactors = FALSE)
  
  # Map UniProt accessions to gene symbols
  gene_symbols <- map_to_gene_symbol(peptideatlas$biosequence_accession, 
                                    peptideatlas$biosequence_desc, 
                                    database_type = "peptideatlas")
  
  protein_data$PeptideAtlas <- data.frame(
    gene_symbol = gene_symbols,
    description = peptideatlas$biosequence_desc,
    abundance = as.numeric(peptideatlas$norm_PSMs_per_100K),
    database = "PeptideAtlas",
    cell_type = "Plasma",
    stringsAsFactors = FALSE
  )
  
  # 2. PaxDb files
  message("  - Reading PaxDb files...")
  paxdb_files <- list.files(pattern = "^paxdb_.*\\.csv$")
  for(file in paxdb_files) {
    cell_type <- gsub("paxdb_|\\.csv", "", file)
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
        protein_ids <- gsub("^9606\\.", "", temp_data$string_external_id)
        
        # Only keep rows with valid protein IDs and abundance values
        valid_rows <- !is.na(protein_ids) & protein_ids != "" & 
                     !is.na(temp_data$abundance) & temp_data$abundance > 0
        
        if(sum(valid_rows) > 0) {
          # Map PaxDb protein accessions to gene symbols
          gene_symbols <- map_to_gene_symbol(protein_ids[valid_rows], database_type = "paxdb")
        
        protein_data[[paste0("PaxDb_", cell_type)]] <- data.frame(
          gene_symbol = gene_symbols,
          description = gene_symbols,  # PaxDb doesn't have preferred_name column
          abundance = temp_data$abundance[valid_rows],
          database = "PaxDb",
          cell_type = cell_type,
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
  hpa_files <- list.files(pattern = "^hpa_.*\\.csv$")
  for(file in hpa_files) {
    technique <- gsub("hpa_|\\.csv", "", file)
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
      concentration_values <- suppressWarnings(as.numeric(conc_strings))
      
      # Remove rows where concentration conversion failed
      numeric_valid <- !is.na(concentration_values) & concentration_values > 0
      
      if(sum(numeric_valid) > 0) {
        # Get description column if available
        if("Protein.name" %in% colnames(valid_data)) {
          descriptions <- valid_data$Protein.name[numeric_valid]
        } else if("Gene.description" %in% colnames(valid_data)) {
          descriptions <- valid_data$Gene.description[numeric_valid]
        } else {
          descriptions <- valid_data$Gene[numeric_valid]  # Use gene name as description
        }
        
        # HPA already has gene symbols
        gene_symbols <- map_to_gene_symbol(valid_data$Gene[numeric_valid], database_type = "hpa")
        
        protein_data[[paste0("HPA_", technique)]] <- data.frame(
          gene_symbol = gene_symbols,
          description = descriptions,
          abundance = concentration_values[numeric_valid],
          database = "HPA",
          cell_type = technique,
          stringsAsFactors = FALSE
        )
      } else {
        message(paste("      Warning: No valid concentration data in", file))
      }
    } else {
      message(paste("      Warning: No valid data found in", file))
    }
  }
  
  # 4. GPMDB files
  message("  - Reading GPMDB files...")
  gpmdb_files <- list.files(pattern = "^gpmdb_.*\\.csv$")
  for(file in gpmdb_files) {
    cell_type <- gsub("gpmdb_|\\.csv", "", file)
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
        # Map GPMDB protein accessions (ENSP format) to gene symbols
        gene_symbols <- map_to_gene_symbol(gpmdb_data$accession[valid_rows], 
                                          gpmdb_data$description[valid_rows], 
                                          database_type = "gpmdb")
        
        protein_data[[paste0("GPMDB_", cell_type)]] <- data.frame(
          gene_symbol = gene_symbols,
          description = gpmdb_data$description[valid_rows],
          abundance = gpmdb_data$total[valid_rows],
          database = "GPMDB",
          cell_type = cell_type,
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
    message(paste("    Processing", file, "for cell type-specific data"))
    
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
              valid_descriptions <- character()
              
              for(i in 1:nrow(pxd_data)) {
                gene_str <- as.character(pxd_data[[gene_col]][i])
                intensity <- mean_intensity[i]
                
                # Get description if Protein.names column exists, otherwise use gene name
                if("Protein.names" %in% colnames(pxd_data)) {
                  description <- as.character(pxd_data$Protein.names[i])
                } else {
                  description <- gene_str  # Use gene name as fallback
                }
                
                # Check for single gene and valid intensity
                if(!is.na(gene_str) && gene_str != "" && !grepl("[;,]", gene_str) &&
                   !is.na(intensity) && intensity > 0) {
                  valid_genes <- c(valid_genes, toupper(trimws(gene_str)))
                  valid_intensities <- c(valid_intensities, intensity)
                  valid_descriptions <- c(valid_descriptions, if(is.na(description)) gene_str else description)
                }
              }
              
              if(length(valid_genes) > 0) {
                protein_data[[paste0("PXD_", cell_type, "_", gsub("\\.csv", "", file))]] <- data.frame(
                  gene_symbol = valid_genes,
                  description = valid_descriptions,
                  abundance = valid_intensities,
                  database = "ProteomeXchange",
                  cell_type = cell_type,
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
                      !is.na(all_data$description) &
                      all_data$description != "", ]
  
  message(paste("Total protein entries loaded:", nrow(all_data)))
  return(all_data)
}

# Function to classify proteins by functional categories
classify_proteins <- function(protein_data) {
  message("Classifying proteins by functional categories...")
  
  # Create functional classifications based on description keywords
  protein_data$functional_category <- "Other"
  
  # Immune system proteins
  immune_keywords <- c("immunoglobulin", "antibody", "complement", "interferon", 
                       "interleukin", "cytokine", "chemokine", "CD\\d+", "TCR", "BCR",
                       "HLA-", "immune", "lymphocyte", "T cell", "B cell")
  immune_pattern <- paste(immune_keywords, collapse = "|")
  protein_data$functional_category[grepl(immune_pattern, protein_data$description, ignore.case = TRUE)] <- "Immune System"
  
  # Blood coagulation and hemostasis
  coagulation_keywords <- c("fibrinogen", "factor", "thrombin", "plasmin", "coagulation",
                           "hemostasis", "platelet", "von Willebrand", "prothrombin")
  coag_pattern <- paste(coagulation_keywords, collapse = "|")
  protein_data$functional_category[grepl(coag_pattern, protein_data$description, ignore.case = TRUE)] <- "Coagulation/Hemostasis"
  
  # Transport proteins
  transport_keywords <- c("albumin", "transferrin", "hemoglobin", "haptoglobin", 
                         "ceruloplasmin", "transporter", "carrier", "binding protein")
  transport_pattern <- paste(transport_keywords, collapse = "|")
  protein_data$functional_category[grepl(transport_pattern, protein_data$description, ignore.case = TRUE)] <- "Transport"
  
  # Enzymes
  enzyme_keywords <- c("kinase", "phosphatase", "dehydrogenase", "oxidase", "reductase",
                      "protease", "peptidase", "synthetase", "transferase", "ligase",
                      "lyase", "isomerase", "hydrolase")
  enzyme_pattern <- paste(enzyme_keywords, collapse = "|")
  protein_data$functional_category[grepl(enzyme_pattern, protein_data$description, ignore.case = TRUE)] <- "Enzymes"
  
  # Structural proteins
  structural_keywords <- c("collagen", "elastin", "keratin", "actin", "myosin", 
                          "tubulin", "spectrin", "fibronectin", "laminin")
  structural_pattern <- paste(structural_keywords, collapse = "|")
  protein_data$functional_category[grepl(structural_pattern, protein_data$description, ignore.case = TRUE)] <- "Structural"
  
  # Signaling proteins
  signaling_keywords <- c("receptor", "growth factor", "hormone", "signal", "pathway",
                         "cascade", "activation", "inhibitor", "regulator")
  signaling_pattern <- paste(signaling_keywords, collapse = "|")
  protein_data$functional_category[grepl(signaling_pattern, protein_data$description, ignore.case = TRUE)] <- "Signaling"
  
  # Metabolic proteins
  metabolic_keywords <- c("metabolism", "metabolic", "glycolysis", "gluconeogenesis",
                         "citric acid", "electron transport", "ATP synthase", "respiratory")
  metabolic_pattern <- paste(metabolic_keywords, collapse = "|")
  protein_data$functional_category[grepl(metabolic_pattern, protein_data$description, ignore.case = TRUE)] <- "Metabolism"
  
  return(protein_data)
}

# Function to identify immunoglobulin proteins
identify_immunoglobulins <- function(protein_data) {
  message("Identifying immunoglobulin proteins...")
  
  # Patterns for different immunoglobulin types
  ig_patterns <- list(
    "IgG" = "immunoglobulin.*heavy.*gamma|IgG|IGHG",
    "IgA" = "immunoglobulin.*heavy.*alpha|IgA|IGHA",
    "IgM" = "immunoglobulin.*heavy.*mu|IgM|IGHM",
    "IgE" = "immunoglobulin.*heavy.*epsilon|IgE|IGHE",
    "IgD" = "immunoglobulin.*heavy.*delta|IgD|IGHD",
    "Light_chain_kappa" = "immunoglobulin.*kappa|IGK|light.*kappa",
    "Light_chain_lambda" = "immunoglobulin.*lambda|IGL|light.*lambda"
  )
  
  protein_data$ig_type <- "Non-Ig"
  
  for(ig_type in names(ig_patterns)) {
    matches <- grepl(ig_patterns[[ig_type]], protein_data$description, ignore.case = TRUE)
    protein_data$ig_type[matches] <- ig_type
  }
  
  return(protein_data)
}

# Read protein data
protein_data <- read_protein_data()

# Classify proteins functionally
protein_data <- classify_proteins(protein_data)

# Identify immunoglobulins
protein_data <- identify_immunoglobulins(protein_data)

# Plot 1: Functional Category Distribution
message("Creating functional category distribution plot...")
func_summary <- protein_data %>%
  group_by(functional_category, database) %>%
  summarise(count = n_distinct(gene_symbol), .groups = 'drop')

tiff(file.path(output_dir, "functional_category_distribution.tiff"), 
     width = 12, height = 8, units = "in", res = 600, compression = "lzw")

ggplot(func_summary, aes(x = reorder(functional_category, count, sum), 
                        y = count, fill = database)) +
  geom_col(position = "stack", alpha = 0.8) +
  coord_flip() +
  scale_fill_brewer(type = "qual", palette = "Set3", name = "Database") +
  theme_minimal() +
  labs(title = "Distribution of Proteins by Functional Category",
       subtitle = "Stacked by database contribution",
       x = "Functional Category", y = "Number of Unique Proteins") +
  theme(plot.title = element_text(size = 14, face = "bold"))

dev.off()

# Plot 2: Immunoglobulin Analysis
message("Creating immunoglobulin analysis...")
ig_data <- protein_data %>% filter(ig_type != "Non-Ig")

if(nrow(ig_data) > 0) {
  ig_summary <- ig_data %>%
    group_by(ig_type, database, cell_type) %>%
    summarise(count = n_distinct(gene_symbol), 
              mean_abundance = mean(abundance, na.rm = TRUE),
              .groups = 'drop')
  
  tiff(file.path(output_dir, "immunoglobulin_analysis.tiff"), 
       width = 14, height = 10, units = "in", res = 600, compression = "lzw")
  
  p1 <- ggplot(ig_summary, aes(x = ig_type, y = count, fill = database)) +
    geom_col(position = "dodge", alpha = 0.8) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Immunoglobulin Types by Database",
         x = "Immunoglobulin Type", y = "Number of Proteins") +
    scale_fill_brewer(type = "qual", palette = "Set2", name = "Database")
  
  p2 <- ggplot(ig_summary, aes(x = ig_type, y = mean_abundance, 
                              color = database, size = count)) +
    geom_point(alpha = 0.7) +
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Immunoglobulin Abundance by Type",
         x = "Immunoglobulin Type", y = "Mean Abundance (log scale)") +
    scale_color_brewer(type = "qual", palette = "Set2", name = "Database") +
    scale_size_continuous(name = "Count")
  
  grid.arrange(p1, p2, nrow = 2,
               top = "Immunoglobulin Protein Analysis")
  
  dev.off()
}

# Plot 3: Functional Word Cloud
message("Creating functional word cloud...")
# Extract keywords from protein descriptions
all_descriptions <- paste(protein_data$description, collapse = " ")
# Clean text
all_descriptions <- gsub("[^A-Za-z\\s]", " ", all_descriptions)
all_descriptions <- tolower(all_descriptions)

# Split into words and count
words <- unlist(strsplit(all_descriptions, "\\s+"))
words <- words[nchar(words) > 3]  # Keep only words longer than 3 characters

# Remove common stop words
stop_words <- c("protein", "human", "isoform", "chain", "subunit", "family", 
                "member", "type", "form", "like", "related", "domain", "containing")
words <- words[!words %in% stop_words]

word_freq <- table(words)
word_freq <- sort(word_freq, decreasing = TRUE)
word_freq <- word_freq[1:min(100, length(word_freq))]  # Top 100 words

tiff(file.path(output_dir, "functional_wordcloud.tiff"), 
     width = 12, height = 8, units = "in", res = 600, compression = "lzw")

wordcloud(names(word_freq), as.numeric(word_freq), 
          max.words = 100, random.order = FALSE, rot.per = 0.35,
          colors = brewer.pal(8, "Dark2"), scale = c(4, 0.5))
title("Most Common Terms in Protein Descriptions", line = -1)

dev.off()

# Plot 4: Database Specialization by Function
message("Creating database functional specialization...")
db_func_spec <- protein_data %>%
  group_by(database, functional_category) %>%
  summarise(count = n_distinct(gene_symbol), .groups = 'drop') %>%
  group_by(database) %>%
  mutate(total = sum(count),
         percentage = count / total * 100) %>%
  ungroup()

tiff(file.path(output_dir, "database_functional_specialization.tiff"), 
     width = 14, height = 10, units = "in", res = 600, compression = "lzw")

ggplot(db_func_spec, aes(x = database, y = percentage, fill = functional_category)) +
  geom_col(position = "fill", alpha = 0.8) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(type = "qual", palette = "Set3", name = "Functional\nCategory") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Database Functional Specialization",
       subtitle = "Relative proportion of functional categories per database",
       x = "Database", y = "Percentage of Proteins")

dev.off()

# Plot 5: Cell Type Functional Profiles
message("Creating cell type functional profiles...")
celltype_func <- protein_data %>%
  filter(database %in% c("PaxDb", "GPMDB")) %>%  # Use databases with multiple cell types
  group_by(cell_type, functional_category) %>%
  summarise(count = n_distinct(gene_symbol), .groups = 'drop') %>%
  group_by(cell_type) %>%
  mutate(total = sum(count),
         percentage = count / total * 100) %>%
  ungroup() %>%
  filter(total >= 100)  # Only cell types with at least 100 proteins

if(nrow(celltype_func) > 0) {
  tiff(file.path(output_dir, "celltype_functional_profiles.tiff"), 
       width = 14, height = 10, units = "in", res = 600, compression = "lzw")
  
  ggplot(celltype_func, aes(x = reorder(cell_type, total), y = percentage, 
                           fill = functional_category)) +
    geom_col(position = "fill", alpha = 0.8) +
    coord_flip() +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_brewer(type = "qual", palette = "Set3", name = "Functional\nCategory") +
    theme_minimal() +
    labs(title = "Cell Type Functional Profiles",
         subtitle = "Relative proportion of functional categories per cell type",
         x = "Cell Type", y = "Percentage of Proteins")
  
  dev.off()
}

# Plot 6: Protein Class Treemap
message("Creating protein class treemap...")
# Simplify categories for treemap
treemap_data <- protein_data %>%
  group_by(functional_category) %>%
  summarise(
    count = n_distinct(gene_symbol),
    avg_abundance = mean(abundance, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(count >= 10)  # Only categories with at least 10 proteins

tiff(file.path(output_dir, "protein_class_treemap.tiff"), 
     width = 12, height = 8, units = "in", res = 600, compression = "lzw")

treemap(treemap_data,
        index = "functional_category",
        vSize = "count",
        vColor = "avg_abundance",
        type = "value",
        palette = "RdYlBu",
        title = "Protein Classes by Count and Average Abundance",
        fontsize.title = 14,
        fontsize.labels = 10)

dev.off()

# Save analysis results
write.csv(func_summary, file.path(output_dir, "functional_category_summary.csv"), 
          row.names = FALSE)

if(nrow(ig_data) > 0) {
  write.csv(ig_summary, file.path(output_dir, "immunoglobulin_summary.csv"), 
            row.names = FALSE)
}

write.csv(head(word_freq, 50), file.path(output_dir, "top_protein_keywords.csv"))

# Summary report
sink(file.path(output_dir, "functional_analysis_summary.txt"))
cat("Functional Analysis Summary\n")
cat("==========================\n\n")

cat("Dataset Overview:\n")
cat(paste("Total unique proteins analyzed:", length(unique(protein_data$gene_symbol)), "\n"))
cat(paste("Databases included:", length(unique(protein_data$database)), "\n"))
cat(paste("Cell types/conditions:", length(unique(protein_data$cell_type)), "\n\n"))

cat("Functional Category Distribution:\n")
func_dist <- table(protein_data$functional_category)
for(cat in names(sort(func_dist, decreasing = TRUE))) {
  percentage <- round(func_dist[cat] / sum(func_dist) * 100, 1)
  cat(paste("  ", cat, ":", func_dist[cat], "proteins (", percentage, "%)\n"))
}

if(nrow(ig_data) > 0) {
  cat("\nImmunoglobulin Analysis:\n")
  ig_dist <- table(ig_data$ig_type)
  for(ig in names(sort(ig_dist, decreasing = TRUE))) {
    cat(paste("  ", ig, ":", ig_dist[ig], "proteins\n"))
  }
}

cat("\nTop 10 Most Frequent Keywords:\n")
for(i in 1:min(10, length(word_freq))) {
  cat(paste("  ", i, ".", names(word_freq)[i], "- appeared", word_freq[i], "times\n"))
}

sink()

message("Functional analysis completed!")
message(paste("Results saved to:", output_dir)) 
