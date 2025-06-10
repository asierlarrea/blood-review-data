#!/usr/bin/env Rscript

# Fast Local-Only ID Mapping (No biomaRt)
# Author: Data Analysis Pipeline
# Date: 2024

# Function for fast local ID mapping
fast_map_to_gene_symbol <- function(accessions, descriptions = NULL) {
  if (length(accessions) == 0) {
    return(character(0))
  }
  
  message(paste("Fast mapping", length(accessions), "IDs to gene symbols..."))
  gene_symbols <- character(length(accessions))
  
  # Process each accession
  for(i in seq_along(accessions)) {
    acc <- as.character(accessions[i])
    desc <- if(!is.null(descriptions)) as.character(descriptions[i]) else ""
    
    # Identify accession type and map accordingly
    gene_symbols[i] <- fast_map_single_id(acc, desc)
  }
  
  # Count results
  mapped_count <- sum(!is.na(gene_symbols) & gene_symbols != "")
  message(paste("Fast mapped:", mapped_count, "out of", length(accessions), "IDs"))
  
  return(gene_symbols)
}

# Fast mapping for single ID
fast_map_single_id <- function(accession, description = "") {
  acc <- as.character(accession)
  desc <- as.character(description)
  
  # 1. If it's already a gene symbol
  if (grepl("^[A-Z][A-Z0-9_-]{1,14}$", acc) && !grepl("^ENS|^[OPQ][0-9]|^[A-Z][0-9]", acc)) {
    return(toupper(acc))
  }
  
  # 2. Common protein mappings
  common_gene <- get_common_protein_mapping(acc)
  if (!is.na(common_gene)) {
    return(common_gene)
  }
  
  # 3. Parse from description
  if (!is.na(desc) && desc != "") {
    desc_gene <- extract_gene_from_description_fast(desc)
    if (!is.na(desc_gene)) {
      return(desc_gene)
    }
  }
  
  # 4. For ENSP IDs, try simple pattern matching in description
  if (grepl("^ENSP[0-9]+", acc) && !is.na(desc) && desc != "") {
    ensp_gene <- extract_gene_from_ensembl_desc(desc)
    if (!is.na(ensp_gene)) {
      return(ensp_gene)
    }
  }
  
  return(NA)
}

# Comprehensive common protein mappings
get_common_protein_mapping <- function(accession) {
  # Remove isoform suffixes
  clean_acc <- gsub("-[0-9]+$", "", accession)
  
  # Extensive common protein mappings
  mappings <- list(
    # ENSP IDs for top plasma proteins (from PaxDB)
    "ENSP00000295897" = "ALB",      # Albumin (highest abundance: 81110)
    "ENSP00000236850" = "TF",       # Transferrin (40718)
    "ENSP00000431254" = "HP",       # Haptoglobin (38503)
    "ENSP00000360522" = "APOA1",    # Apolipoprotein A-I (29729)
    "ENSP00000356969" = "FGG",      # Fibrinogen gamma chain (29070)
    "ENSP00000465356" = "AHSG",     # Alpha-2-HS-glycoprotein (27221)
    "ENSP00000416066" = "A2M",      # Alpha-2-macroglobulin (26440)
    "ENSP00000237014" = "CP",       # Ceruloplasmin (19275)
    "ENSP00000227667" = "TTR",      # Transthyretin (18751)
    "ENSP00000385834" = "GC",       # Vitamin D-binding protein (18332)
    "ENSP00000259396" = "PLG",      # Plasminogen (17930)
    "ENSP00000348170" = "C3",       # Complement C3 (17194)
    "ENSP00000345179" = "FGB",      # Fibrinogen beta chain (16630)
    "ENSP00000370010" = "FGA",      # Fibrinogen alpha chain (16536)
    "ENSP00000306099" = "APOB",     # Apolipoprotein B (16421)
    "ENSP00000265983" = "IGHG1",    # Immunoglobulin heavy constant gamma 1 (16017)
    "ENSP00000356671" = "APOA2",    # Apolipoprotein A-II (15061)
    "ENSP00000362924" = "CFH",      # Complement factor H (14804)
    "ENSP00000323929" = "SERPINC1", # Antithrombin (14781)
    "ENSP00000264613" = "IGHM",     # Immunoglobulin heavy constant mu (14377)
    
    # Additional ENSP mappings from biomaRt (top 50)
    "ENSP00000245907" = "C3",       # abundance: 13768
    "ENSP00000292401" = "AZGP1",    # abundance: 13577
    "ENSP00000393887" = "AHSG",     # abundance: 13235 (duplicate, but keeping)
    "ENSP00000350425" = "APOA4",    # abundance: 12331
    "ENSP00000336829" = "FGG",      # abundance: 12045
    "ENSP00000252490" = "APOC2",    # abundance: 9694
    "ENSP00000466775" = "APOC2",    # abundance: 9052
    "ENSP00000421725" = "GC",       # abundance: 8342
    "ENSP00000255040" = "APCS",     # abundance: 7739
    "ENSP00000394936" = "ORM2",     # abundance: 7509
    "ENSP00000333994" = "HBB",      # abundance: 7336
    "ENSP00000263100" = "A1BG",     # abundance: 7068
    "ENSP00000252486" = "APOE",     # abundance: 7000
    "ENSP00000321853" = "SERPINF2", # abundance: 6935
    "ENSP00000254722" = "SERPINF1", # abundance: 6806
    "ENSP00000233242" = "APOB",     # abundance: 6693
    "ENSP00000215727" = "SERPIND1", # abundance: 6501
    "ENSP00000302621" = "LRG1",     # abundance: 6393
    "ENSP00000205948" = "APOH",     # abundance: 6167
    "ENSP00000355627" = "AGT",      # abundance: 6153
    "ENSP00000345968" = "PGLYRP2",  # abundance: 6016
    "ENSP00000263408" = "C9",       # abundance: 5524
    "ENSP00000296130" = "CLEC3B",   # abundance: 5394
    "ENSP00000376793" = "SERPINA3", # abundance: 5184
    "ENSP00000265132" = "AMBP",     # abundance: 5100
    "ENSP00000251595" = "HBA2",     # abundance: 5088
    "ENSP00000356399" = "CFH",      # abundance: 5081
    "ENSP00000266041" = "ITIH4",    # abundance: 4913
    "ENSP00000315130" = "CLU",      # abundance: 4809
    "ENSP00000252244" = "KRT1",     # abundance: 4684
    "ENSP00000266718" = "LUM",      # abundance: 4239
    "ENSP00000329374" = "SERPINA7", # abundance: 4188
    "ENSP00000232003" = "HRG",      # abundance: 4175
    "ENSP00000278407" = "SERPING1", # abundance: 4159
    "ENSP00000396688" = "C4A",      # abundance: 4157
    "ENSP00000416561" = "CFB",      # abundance: 3952
    "ENSP00000246662" = "KRT9",     # abundance: 3920
    "ENSP00000450540" = "SERPINA3", # abundance: 3852
    "ENSP00000222381" = "PON1",     # abundance: 3848
    "ENSP00000308938" = "PLG",      # abundance: 3825
    "ENSP00000308541" = "F2",       # abundance: 3820
    "ENSP00000273283" = "ITIH1",    # abundance: 3767
    "ENSP00000351190" = "ITIH2",    # abundance: 3627
    "ENSP00000322061" = "C7",       # abundance: 3549
    "ENSP00000342850" = "SERPINA6", # abundance: 3497
    "ENSP00000226218" = "VTN",      # abundance: 3363
    "ENSP00000278222" = "SAA4",     # abundance: 3230
    "ENSP00000356037" = "C4BPA",    # abundance: 3169
    "ENSP00000384906" = "SAA1",     # abundance: 3155
    "ENSP00000322421" = "HBA1",     # abundance: 2965
    
    # Blood/Plasma proteins (UniProt IDs)
    "P02768" = "ALB",      # Albumin
    "P02787" = "TF",       # Transferrin  
    "P00738" = "HP",       # Haptoglobin
    "P02765" = "AHSG",     # Alpha-2-HS-glycoprotein
    "P02647" = "APOA1",    # Apolipoprotein A-I
    "P04217" = "A1BG",     # Alpha-1B-glycoprotein
    "P01023" = "A2M",      # Alpha-2-macroglobulin
    "P02671" = "FGA",      # Fibrinogen alpha
    "P02675" = "FGB",      # Fibrinogen beta
    "P02679" = "FGG",      # Fibrinogen gamma
    "P00450" = "CP",       # Ceruloplasmin
    "P02766" = "TTR",      # Transthyretin
    "P02774" = "GC",       # Vitamin D-binding protein
    "P00747" = "PLG",      # Plasminogen
    "P01024" = "C3",       # Complement C3
    "P0C0L4" = "C4A",      # Complement C4-A
    "P0C0L5" = "C4B",      # Complement C4-B
    "P01031" = "C5",       # Complement C5
    
    # Immunoglobulins
    "P01857" = "IGHG1",    # Immunoglobulin heavy constant gamma 1
    "P01859" = "IGHG2",    # Immunoglobulin heavy constant gamma 2
    "P01860" = "IGHG3",    # Immunoglobulin heavy constant gamma 3
    "P01861" = "IGHG4",    # Immunoglobulin heavy constant gamma 4
    "P01876" = "IGHA1",    # Immunoglobulin heavy constant alpha 1
    "P01877" = "IGHA2",    # Immunoglobulin heavy constant alpha 2
    "P01871" = "IGHM",     # Immunoglobulin heavy constant mu
    "P01834" = "IGKC",     # Immunoglobulin kappa constant
    "P0DOX5" = "IGLC2",    # Immunoglobulin lambda constant 2
    
    # Metabolic enzymes
    "P04406" = "GAPDH",    # Glyceraldehyde-3-phosphate dehydrogenase
    "P60174" = "TPI1",     # Triosephosphate isomerase
    "P00558" = "PGK1",     # Phosphoglycerate kinase 1
    "P06733" = "ENO1",     # Enolase 1
    "P14618" = "PKM",      # Pyruvate kinase M1/2
    "P53396" = "ACLY",     # ATP citrate lyase
    
    # Cytoskeletal proteins
    "P60709" = "ACTB",     # Actin beta
    "P63261" = "ACTG1",    # Actin gamma 1
    "P68133" = "ACTA1",    # Actin alpha 1
    "P68032" = "ACTC1",    # Actin alpha cardiac muscle 1
    "P07437" = "TUBB",     # Tubulin beta
    "P68363" = "TUBA1A",   # Tubulin alpha 1a
    
    # Ribosomal proteins
    "P62805" = "RPS6",     # Ribosomal protein S6
    "P46781" = "RPS9",     # Ribosomal protein S9
    "P62277" = "RPS13",    # Ribosomal protein S13
    "P62701" = "RPS4X",    # Ribosomal protein S4 X-linked
    "P62081" = "RPS7",     # Ribosomal protein S7
    
    # Heat shock proteins
    "P07900" = "HSP90AA1", # Heat shock protein 90 alpha A1
    "P08238" = "HSP90AB1", # Heat shock protein 90 alpha B1
    "P11142" = "HSPA8",    # Heat shock cognate 71 kDa protein
    "P0DMV8" = "HSPA1A",   # Heat shock 70 kDa protein 1A
    
    # Common enzymes
    "P04040" = "CAT",      # Catalase
    "P00441" = "SOD1",     # Superoxide dismutase 1
    "P04179" = "SOD2",     # Superoxide dismutase 2
    "P35573" = "GPX1",     # Glutathione peroxidase 1
    "P06280" = "AGXT",     # Alanine-glyoxylate aminotransferase
    
    # Transcription factors
    "Q13438" = "OS9",      # Protein OS-9
    "P04150" = "NR3C1",    # Glucocorticoid receptor
    "P03372" = "ESR1",     # Estrogen receptor 1
    
    # Other important proteins
    "P02545" = "LMNA",     # Prelamin-A/C
    "P07195" = "LDHB",     # L-lactate dehydrogenase B chain
    "P00338" = "LDHA",     # L-lactate dehydrogenase A chain
    "P04075" = "ALDOA",    # Fructose-bisphosphate aldolase A
    "P08133" = "ANXA6"     # Annexin A6
  )
  
  result <- mappings[[clean_acc]]
  return(if(is.null(result)) NA else result)
}

# Fast description parsing
extract_gene_from_description_fast <- function(description) {
  if (is.na(description) || description == "") {
    return(NA)
  }
  
  desc <- as.character(description)
  
  # Pattern 1: GN=GENE_NAME (most common in UniProt)
  if (grepl("GN=([A-Za-z0-9_-]+)", desc)) {
    gene <- gsub(".*GN=([A-Za-z0-9_-]+).*", "\\1", desc)
    if (nchar(gene) >= 2 && nchar(gene) <= 15 && !grepl("[^A-Za-z0-9_-]", gene)) {
      return(toupper(gene))
    }
  }
  
  # Pattern 2: Gene names in parentheses
  if (grepl("\\(([A-Z][A-Z0-9_-]{1,14})\\)", desc)) {
    matches <- regmatches(desc, gregexpr("\\(([A-Z][A-Z0-9_-]{1,14})\\)", desc))[[1]]
    for(match in matches) {
      gene <- gsub("[()]", "", match)
      if (!grepl("PROTEIN|FRAGMENT|DOMAIN|FAMILY", gene, ignore.case = TRUE)) {
        return(toupper(gene))
      }
    }
  }
  
  # Pattern 3: Gene names in brackets
  if (grepl("\\[([A-Z][A-Z0-9_-]{1,14})\\]", desc)) {
    matches <- regmatches(desc, gregexpr("\\[([A-Z][A-Z0-9_-]{1,14})\\]", desc))[[1]]
    for(match in matches) {
      gene <- gsub("[\\[\\]]", "", match)
      if (!grepl("PROTEIN|FRAGMENT|DOMAIN|FAMILY", gene, ignore.case = TRUE)) {
        return(toupper(gene))
      }
    }
  }
  
  # Pattern 4: \GName=GENE_NAME
  if (grepl("\\\\GName=([A-Za-z0-9_-]+)", desc)) {
    gene <- gsub(".*\\\\GName=([A-Za-z0-9_-]+).*", "\\1", desc)
    if (nchar(gene) >= 2 && nchar(gene) <= 15) {
      return(toupper(gene))
    }
  }
  
  # Pattern 5: First word if it looks like a gene symbol
  words <- unlist(strsplit(desc, "[\\s,;]+"))
  for(word in words[1:min(3, length(words))]) {
    clean_word <- gsub("[^A-Za-z0-9_-]", "", word)
    if (nchar(clean_word) >= 2 && nchar(clean_word) <= 15 && 
        grepl("^[A-Z][A-Z0-9_-]*$", clean_word) &&
        !grepl("PROTEIN|GENE|FAMILY|DOMAIN|FRAGMENT|UNCHARACTERIZED", clean_word, ignore.case = TRUE)) {
      return(toupper(clean_word))
    }
  }
  
  return(NA)
}

# Extract gene from Ensembl descriptions
extract_gene_from_ensembl_desc <- function(description) {
  if (is.na(description) || description == "") {
    return(NA)
  }
  
  desc <- as.character(description)
  
  # Ensembl descriptions often have gene symbols at the start
  # Example: "ALB albumin [Source:HGNC Symbol;Acc:HGNC:399]"
  words <- unlist(strsplit(desc, "[\\s\\[]+"))
  first_word <- words[1]
  
  if (nchar(first_word) >= 2 && nchar(first_word) <= 15 && 
      grepl("^[A-Z][A-Z0-9_-]*$", first_word) &&
      !grepl("PROTEIN|GENE|FAMILY|DOMAIN|FRAGMENT", first_word, ignore.case = TRUE)) {
    return(toupper(first_word))
  }
  
  # Look for gene symbols in brackets or after "Symbol;"
  if (grepl("Symbol;Acc:", desc)) {
    # Pattern: [Source:HGNC Symbol;Acc:HGNC:399]
    symbol_pattern <- "Symbol;Acc:HGNC:([0-9]+)"
    if (grepl(symbol_pattern, desc)) {
      # This would need HGNC mapping, skip for now
      return(NA)
    }
  }
  
  # Look for obvious gene patterns
  if (grepl("\\[([A-Z][A-Z0-9_-]{1,14})\\]", desc)) {
    gene <- gsub(".*\\[([A-Z][A-Z0-9_-]{1,14})\\].*", "\\1", desc)
    if (!grepl("PROTEIN|FRAGMENT|DOMAIN|SOURCE|HGNC", gene, ignore.case = TRUE)) {
      return(toupper(gene))
    }
  }
  
  return(NA)
}

# Test function
test_fast_mapping <- function() {
  message("Testing fast ID mapping...")
  
  test_ids <- c(
    "P02768",           # UniProt ALB
    "ENSP00000295897",  # Ensembl ALB
    "ALB",              # Gene symbol
    "Q13438",           # UniProt OS9
    "INVALID"           # Invalid
  )
  
  test_descriptions <- c(
    "Serum albumin GN=ALB PE=1 SV=2",
    "ALB albumin [Source:HGNC Symbol;Acc:HGNC:399]", 
    "",
    "Protein OS-9 GN=OS9 PE=1 SV=1",
    ""
  )
  
  results <- fast_map_to_gene_symbol(test_ids, test_descriptions)
  
  for (i in seq_along(test_ids)) {
    message(paste(test_ids[i], "->", 
                  if (is.na(results[i])) "NOT_MAPPED" else results[i]))
  }
  
  return(results)
}

message("Fast ID mapping functions loaded:")
message("- fast_map_to_gene_symbol(accessions, descriptions)")
message("- test_fast_mapping()") 