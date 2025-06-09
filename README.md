# Blood Proteomics Analysis

A comprehensive analysis of blood cell proteomes across multiple databases and cell types, focusing on plasma proteins and cellular protein distributions.

## Overview

This project analyzes proteomics data from various databases to understand:
- Protein distributions across different blood cell types
- Database coverage and overlap patterns
- Intensity distributions of proteins in various cellular contexts
- Comparative proteome analysis between immune cell populations

## Project Structure

```
blood-review-data/
├── README.md                       # This file
├── LICENSE                         # Apache 2.0 License
├── .gitignore                      # Git ignore patterns
├── Figure2.R                       # UpSet plot of plasma proteins
├── Figure3.R                       # Cell type bubble plot
├── Figure4.R                       # Cell type UpSet plot
├── Figure5.R                       # Intensity distribution boxplots
├── Data.docx                       # Documentation
└── data/                           # CSV data files
    ├── PeptideAtlas.csv
    ├── PaxDb_*.csv                 # PaxDb database files
    ├── PXD*.csv                    # ProteomeXchange datasets
    ├── HPA_*.csv                   # Human Protein Atlas files
    └── GPMDB_*.csv                 # GPMDB database files
```

## Data Sources

The analysis incorporates data from several major proteomics databases:

- **PeptideAtlas**: Comprehensive peptide repository
- **PaxDb**: Protein abundance database
- **Human Protein Atlas (HPA)**: Multiple assay types (MS, PEA, Immunoassay)
- **GPMDB**: Global Proteome Machine Database
- **ProteomeXchange**: Public proteomics data repository

### Cell Types Analyzed

- **Adaptive Immune**: CD8+ T cells, CD4+ T cells, B cells
- **Innate Immune**: NK cells, Monocytes, Dendritic cells, Macrophages
- **Granulocytes**: Neutrophils, Eosinophils, Basophils
- **Blood Cells**: Platelets, Erythrocytes

## Requirements

### R Dependencies

The analysis requires R (≥ 4.0.0) with the following packages:

```r
# Core packages
install.packages(c(
  "ggplot2",      # Data visualization
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "readr",        # File reading
  "scales",       # Scale functions for ggplot2
  "UpSetR",       # UpSet plot generation
  "grid"          # Grid graphics
))
```

All scripts include automatic package installation, so missing packages will be installed when you run the scripts.

## Usage

### Quick Start

1. **Clone or download** this repository
2. **Set working directory** to the project folder in R:
   ```r
   setwd("path/to/blood-review-data")
   ```
3. **Run individual scripts** to generate figures:
   ```r
   source("Figure2.R")  # Plasma protein UpSet plot
   source("Figure3.R")  # Cell type bubble plot
   source("Figure4.R")  # Cell type UpSet plot
   source("Figure5.R")  # Intensity distribution boxplots
   ```

### Script Details

#### Figure2.R - Plasma Protein Database Overlap
- **Purpose**: Shows protein overlap between databases using UpSet plots
- **Input**: `combined_proteome.csv`
- **Output**: Interactive UpSet plot + summary statistics
- **Key Features**: Database intersection analysis, frequency ordering

#### Figure3.R - Cell Type Coverage by Database
- **Purpose**: Bubble plot showing protein counts by database and cell type
- **Input**: `cell_types.csv`
- **Output**: Bubble plot with protein count visualization
- **Key Features**: Proportional sizing, database comparison

#### Figure4.R - Cell Type Protein Overlap
- **Purpose**: UpSet plot comparing protein overlap across cell types
- **Input**: `cell_type.csv`
- **Output**: Cell type intersection analysis
- **Key Features**: Cross-cell type comparison, frequency analysis

#### Figure5.R - Protein Intensity Distributions
- **Purpose**: Boxplots showing protein intensity distributions by cell type
- **Input**: `cell_type.csv` (with Intensity_ columns)
- **Output**: Log-scale boxplots with statistics
- **Key Features**: Dynamic statistics, database coverage indicators

## Data Format Requirements

### Expected CSV Structure

All scripts expect specific column naming patterns:

- **Database columns**: Named exactly as listed in `set_vars`
- **Intensity columns**: Must start with `Intensity_` prefix
- **Cell type columns**: Should match expected cell type names
- **Required columns**: `Database`, `Cell.type`, `Protein.count` for Figure3.R

### Example Data Structure

```csv
# combined_proteome.csv
GPMDB,PaxDB,PeptideAtlas,HPA_Immuno,HPA_MS,HPA_PEA
1,0,1,0,1,0
0,1,1,1,0,1
...

# cell_types.csv
Database,Cell.type,Protein.count
GPMDB,CD8,1523
PaxDB,CD4,2341
...

# cell_type.csv (for Figure5.R)
Protein_ID,Intensity_CD8,Intensity_CD4,Intensity_B_cell,...
P12345,2.3,1.8,3.1,...
```

## Features

### Error Handling
- Automatic file existence checking
- Data validation and column verification
- Graceful handling of missing data
- Informative error messages

### Visualization Enhancements
- Professional color schemes
- Consistent styling across plots
- Automatic scaling and formatting
- Interactive plot elements where applicable

### Statistical Analysis
- Automatic summary statistics generation
- Database coverage analysis
- Protein count calculations
- Intensity distribution metrics

## Troubleshooting

### Common Issues

1. **File not found errors**:
   - Ensure you're in the correct working directory
   - Check that CSV files exist and have correct names
   - Verify file paths are relative to working directory

2. **Missing columns**:
   - Check CSV headers match expected column names
   - Ensure intensity columns start with `Intensity_`
   - Verify database names match script expectations

3. **Package installation issues**:
   - Update R to latest version
   - Install packages manually if automatic installation fails
   - Check internet connection for package downloads

### Data Quality Checks

The scripts automatically perform several data quality checks:
- Empty file detection
- Missing column validation
- NA value handling
- Data type validation

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add appropriate tests/validation
5. Submit a pull request

## Citation

If you use this analysis in your research, please cite:

```
[Your Citation Information Here]
Blood Proteomics Analysis. (2024). 
Available at: [Repository URL]
```

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or issues, please [open an issue](../../issues) or contact the project maintainers.

---

## Changelog

### Version 1.0.0
- Initial release
- Complete proteomics analysis pipeline
- Four figure generation scripts
- Comprehensive error handling
- Automated statistics generation 