# Blood Proteomics Database Analysis

A comprehensive analysis pipeline for blood cell proteomes across multiple databases, featuring advanced biological ID mapping, dataset integration, and comparative analysis of plasma proteins and cellular protein distributions.

## 🔬 Overview

This project provides a comprehensive proteomics analysis pipeline that:
- **Integrates 4 major databases** for plasma protein analysis
- **Maps 12,713 unique genes** across blood cell types and plasma
- **Handles complex biological ID mapping** (UniProt, ENSP, gene symbols)
- **Analyzes protein distributions** across 11 different blood cell types
- **Provides robust deduplication** and dataset merging strategies
- **Generates publication-ready visualizations** and correlation matrices

## 🚀 Key Features

### ✅ Advanced Biological ID Mapping
- **Enhanced ENSP mapping**: 96% coverage using biomaRt (7,033/7,328 IDs)
- **Multi-format support**: UniProt, ENSP, gene symbols, descriptions
- **Fast local fallbacks**: 100+ common protein mappings
- **Automatic batch processing**: Handles large datasets efficiently

### ✅ Comprehensive Database Integration
- **4-source plasma analysis**: PeptideAtlas, GPMDB, HPA, PaxDB
- **11 cell types**: CD4/CD8 T cells, B cells, NK, monocytes, neutrophils, etc.
- **Dataset tracking**: Individual accession tracking (PXD004352, PXD025174, etc.)
- **Smart deduplication**: Configurable gene merging across datasets

### ✅ Publication-Ready Analysis
- **Correlation matrices**: Jaccard similarity between databases
- **Overlap analysis**: Gene distribution across multiple sources
- **Stacked visualizations**: Dataset contribution plots
- **Statistical summaries**: Comprehensive coverage reports

## 📁 Project Structure

```
blood-review-data/
├── README.md                           # This documentation
├── Figure6_database_correlation.R      # Main analysis pipeline
├── enhanced_fast_mapping.R             # ID mapping with biomaRt
├── fast_id_mapping.R                   # Local ID mapping functions
├── map_unmapped_ensp.R                 # biomaRt ENSP mapping
├── fix_deduplication.R                 # Deduplication analysis
├── biomart_ensp_mappings.csv           # 7,013 ENSP→gene mappings
│
├── Data Files/
│   ├── PeptideAtlas.csv                # PeptideAtlas plasma data
│   ├── PaxDB_plasma.csv                # PaxDB plasma abundance
│   ├── GPMDB_Plasma.csv                # GPMDB plasma proteins
│   ├── HPA_*.csv                       # Human Protein Atlas (MS, PEA, Immunoassay)
│   ├── PXD*.csv                        # ProteomeXchange datasets
│   └── PaxDb_*.csv                     # PaxDB cell type data
│
└── plots/01_Database_Analysis/         # Analysis outputs
    ├── processed_gene_data_with_categories.csv  # Main results (7.3MB)
    ├── plasma_gene_data.csv             # Plasma proteins (8,066 genes)
    ├── celltype_gene_data.csv           # Cell type proteins (9,578 genes)
    ├── *_correlation.tiff               # Correlation matrices
    ├── *_distribution.tiff              # Gene distribution plots
    └── *_source_counts.tiff             # Database contribution plots
```

## 🎯 Analysis Results

### Database Coverage Summary
| Database | Plasma Genes | Cell Type Genes | Total Mapped |
|----------|-------------|----------------|-------------|
| **PaxDB** | 7,010 | - | 7,010 |
| **PeptideAtlas** | 4,603 | - | 4,603 |
| **HPA** | 4,433 | - | 4,433 |
| **GPMDB** | 2,373 | - | 2,373 |
| **ProteomeXchange** | - | 9,578 | 9,578 |
| **Total Unique** | **8,066** | **9,578** | **12,713** |

### Key Achievements
- **📈 355x improvement** in PaxDB mapping (20 → 7,033 genes)
- **🔗 4-source plasma analysis** with 1,376 genes in all databases
- **⚡ 96% ENSP mapping success** using biomaRt integration
- **🎯 Perfect 1:1 gene ratios** with proper deduplication

## ⚙️ Requirements

### System Requirements
- **R ≥ 4.0.0**
- **Internet connection** (for biomaRt queries)
- **~2GB RAM** (for large dataset processing)
- **~500MB disk space** (for outputs and mappings)

### R Dependencies
All required packages are managed through a global dependency system:

**Core Packages:**
- `ggplot2`, `dplyr`, `tidyr` (data manipulation & visualization)
- `scales`, `viridis`, `RColorBrewer` (plotting aesthetics)
- `gridExtra`, `ggridges`, `ggbeeswarm` (specialized plots)

**Bioconductor Packages:**
- `biomaRt`, `org.Hs.eg.db` (biological ID mapping)
- `GO.db`, `KEGG.db` (functional annotation)

**Specialized Analysis:**
- `UpSetR`, `VennDiagram` (set analysis)
- `corrplot`, `pheatmap` (correlation & heatmaps)

## 🚦 Quick Start

### 1. **Install All Dependencies (Required First Step)**
```bash
# Install all required packages for the entire project
Rscript install_dependencies.R
```

### 2. **Run Complete Analysis**
```r
# Set working directory
setwd("path/to/blood-review-data")

# Run main analysis pipeline
source("Figure6_database_correlation.R")
```

### 3. **Run Specific Analyses**
```r
# Biomarker analysis
source("biomarker_plasma_analysis.R")

# Ranked abundance plots
source("ranked_abundance_plots_improved.R")

# Test enhanced mapping system
source("enhanced_fast_mapping.R")
test_enhanced_mapping()
```

## 📊 Key Analysis Components

### 🧬 Biological ID Mapping Pipeline

**Enhanced Fast Mapping System:**
```r
# Load all 7,013 biomaRt ENSP mappings
source("enhanced_fast_mapping.R")

# Map IDs with multiple strategies
gene_symbols <- enhanced_fast_map_to_gene_symbol(
  accessions = c("P02768", "ENSP00000295897", "ALB"),
  descriptions = c("Albumin GN=ALB", "", "")
)
# Result: "ALB", "ALB", "ALB"
```

**Mapping Strategies (in order):**
1. **Gene symbol recognition** (direct mapping)
2. **biomaRt ENSP lookup** (comprehensive, 7,013 mappings)
3. **Common protein mappings** (local, 100+ proteins)
4. **Description parsing** (multiple patterns)
5. **Ensembl description extraction** (for ENSP IDs)

### 📈 Database Correlation Analysis

**4-Source Plasma Overlap:**
```r
# Analyze gene overlap between databases
plasma_analysis <- analyze_gene_overlaps_by_category(gene_data, "Plasma")

# Results:
# - 1,376 genes in all 4 sources
# - 2,736 genes in 3 sources  
# - 753 genes in 2 sources
# - 3,201 genes in 1 source
```

### 🔄 Deduplication Strategies

**Option 1: Complete Deduplication (Recommended)**
```r
# Merge genes across all datasets (standard approach)
deduplicated_data <- data %>%
  group_by(Gene, Database, Category, CellType) %>%
  summarise(Concentration = max(Concentration), .groups = 'drop')
# Result: 18,876 → 9,352 measurements (1:1 gene ratio)
```

**Option 2: Dataset-Aware Analysis**
```r
# Keep dataset information for meta-analysis
dataset_aware_data <- data %>%
  group_by(Gene, Database, Dataset, Category, CellType) %>%
  slice_max(Concentration, n = 1) %>%
  ungroup()
# Result: Maintains dataset provenance for comparison studies
```

## 🎨 Output Visualizations

### Generated Plots
1. **Correlation Matrices**: Jaccard similarity between databases
2. **Gene Distribution**: Histogram of genes per source count
3. **Source Counts**: Bar plots of database contributions
4. **Stacked Dataset Plots**: Cell type contributions by dataset

### CSV Outputs
- **`processed_gene_data_with_categories.csv`**: Complete dataset (7.3MB)
- **`plasma_gene_data.csv`**: Plasma proteins only (8,066 genes)
- **`celltype_gene_data.csv`**: Cell type proteins (9,578 genes)
- **`biomart_ensp_mappings.csv`**: All ENSP→gene mappings (7,013)

## 📦 Dependency Management

This project uses a centralized dependency management system for clean, maintainable code.

### **Global Installation System**
```bash
# One-time setup: Install all dependencies
Rscript install_dependencies.R
```

**What it does:**
- ✅ Installs all required CRAN packages
- ✅ Handles Bioconductor packages (biomaRt, org.Hs.eg.db)
- ✅ Creates `load_packages.R` for analysis scripts
- ✅ Provides detailed installation summary
- ✅ Handles optional packages gracefully

### **Individual Script Usage**
All analysis scripts now use the streamlined loading system:
```r
# Instead of managing dependencies in each script
source("load_packages.R")
required_packages <- c("ggplot2", "dplyr", "scales")
load_packages(required_packages)
```

### **Benefits of This System**
- 🧹 **Cleaner code**: No dependency management clutter in analysis scripts
- ⚡ **Faster execution**: No repeated installation checks
- 🔄 **Consistency**: Same packages across all scripts
- 🛠️ **Maintainability**: Single point for dependency updates
- 📊 **Better errors**: Clear messages if packages missing

### **File Structure**
```
├── install_dependencies.R      # Global dependency installer
├── load_packages.R            # Generated package loader (auto-created)
├── Figure6_database_correlation.R    # Uses load_packages()
├── biomarker_plasma_analysis.R       # Uses load_packages()
└── ranked_abundance_plots_improved.R # Uses load_packages()
```

## 🔧 Advanced Usage

### Custom ID Mapping
```r
# Map your own accessions
source("enhanced_fast_mapping.R")
my_ids <- c("P02768", "ENSP00000295897", "Q13438")
results <- enhanced_fast_map_to_gene_symbol(my_ids)
```

### Deduplication Analysis
```r
# Analyze gene duplication patterns
source("fix_deduplication.R")
# Generates: cd8_merged_deduplicated.csv, cd8_dataset_aware.csv
```

### Generate New biomaRt Mappings
```r
# Update ENSP mappings with latest biomaRt
source("map_unmapped_ensp.R")
# Creates: biomart_ensp_mappings.csv, biomart_ensp_mappings.R
```

## 📋 Data Format Requirements

### Input Data Structure
```csv
# PaxDB format (organism prefix handling)
string_external_id,abundance
9606.ENSP00000295897,81110
9606.ENSP00000236850,40718

# GPMDB format (description parsing)
accession,description,total
P02768,"Serum albumin GN=ALB PE=1 SV=2",12500

# ProteomeXchange format (gene names)
Gene names,LFQ.intensity.CD8_1,LFQ.intensity.CD4_1
ALB,1.5e+08,2.3e+07
```

### Output Data Consistency
- **Gene symbols**: Uppercase, validated (ALB, TF, HP)
- **Concentrations**: Numeric, highest value per gene
- **Datasets**: Clean accessions (PXD004352, PXD025174, PXD040957)
- **Categories**: "Plasma" or "CellType"

## 🐛 Troubleshooting

### Common Issues

**1. biomaRt Connection Timeout**
```r
# Switch to enhanced fast mapping (local only)
source("fast_id_mapping.R")
gene_symbols <- fast_map_to_gene_symbol(accessions, descriptions)
```

**2. Memory Issues with Large Datasets**
```r
# Process in smaller batches
batch_size <- 100  # Reduce from default 500
```

**3. ENSP IDs Not Mapping**
```r
# Check if biomart_ensp_mappings.csv exists
file.exists("biomart_ensp_mappings.csv")

# Regenerate if needed
source("map_unmapped_ensp.R")
```

**4. Duplicate Gene Counts**
```r
# Check deduplication strategy
source("fix_deduplication.R")  # Analyze the issue
# Then choose Option 1 (complete) or Option 2 (dataset-aware)
```

### Data Quality Checks
- **File validation**: Automatic existence and format checking
- **ID validation**: Multi-format accession recognition
- **Coverage reporting**: Mapping success rates per dataset
- **Duplicate detection**: Gene count vs. measurement count analysis

## 🤝 Contributing

### Development Workflow
1. **Fork** the repository
2. **Create feature branch**: `git checkout -b feature/new-analysis`
3. **Test changes**: Ensure all scripts run successfully
4. **Update documentation**: Modify README for new features
5. **Submit pull request**: Include test results and examples

### Adding New Databases
1. **Create mapping function** in `enhanced_fast_mapping.R`
2. **Add processing logic** in `Figure6_database_correlation.R`
3. **Update common mappings** if needed
4. **Test integration** with existing analysis

## 📚 Citation

If you use this analysis pipeline in your research, please cite:

```bibtex
@software{blood_proteomics_analysis,
  title = {Blood Proteomics Database Analysis Pipeline},
  author = {Data Analysis Team},
  year = {2024},
  url = {https://github.com/your-repo/blood-review-data},
  note = {Comprehensive proteomics analysis with biomaRt integration}
}
```

## 📄 License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## 🔗 Related Resources

- **biomaRt Documentation**: https://bioconductor.org/packages/biomaRt/
- **PeptideAtlas**: http://www.peptideatlas.org/
- **Human Protein Atlas**: https://www.proteinatlas.org/
- **GPMDB**: http://gpmdb.thegpm.org/
- **PaxDB**: https://pax-db.org/

---

### 📈 Analysis Statistics

- **Total Processing Time**: ~15 minutes (including biomaRt queries)
- **Memory Usage**: ~1.5GB peak
- **Output Size**: ~20MB (plots + CSV files)
- **Gene Mapping Success**: 96% overall (12,713/13,247 total) 