# Blood Proteomics Analysis Pipeline

A comprehensive analysis pipeline for blood proteomics data integration and visualization, focusing on database coverage analysis, protein variability studies, and manuscript-ready visualizations.

## 🗂️ Repository Organization

```
blood-review-data/
├── 📁 data/                    # All data files organized by type
│   ├── raw/                   # Original CSV and data files
│   ├── processed/             # Cleaned and processed data  
│   └── metadata/              # Documentation and metadata files
├── 📁 scripts/                # All analysis scripts organized by function
│   ├── analysis/              # Main Figure*.R and analysis scripts
│   ├── visualization/         # Specialized plotting scripts
│   ├── data_processing/       # Data cleaning and ID mapping
│   └── utilities/             # Helper functions and utilities
├── 📁 outputs/                # All generated outputs
│   ├── plots/                 # Generated visualizations
│   ├── tables/                # Generated tables and summaries
│   ├── reports/               # Analysis reports and documentation
│   └── logs/                  # Script execution logs
├── 📁 manuscript/             # Manuscript files and bibliography
├── 📁 docs/                   # Documentation and README files
├── 📁 config/                 # Configuration and dependency files
├── run_analysis.sh            # Main analysis pipeline runner
└── README.md                  # This file
```

## 🚀 Quick Start

### 1. Install Dependencies
```bash
# Install R packages automatically
Rscript config/install_dependencies.R
```

### 2. Run Complete Analysis
```bash
# Run all analysis scripts with organized outputs
./run_analysis.sh
```

### 3. Check Results
All outputs will be organized in the `outputs/` directory:
- **Plots**: `outputs/plots/` - Publication-ready visualizations
- **Tables**: `outputs/tables/` - Summary statistics and data tables
- **Logs**: `outputs/logs/` - Execution logs for debugging

## 📊 Analysis Components

### Core Analysis Scripts (`scripts/analysis/`)
- **Figure2.R** - UpSet plot of plasma proteins
- **Figure3.R** - Cell type bubble plot
- **Figure4.R** - Cell type UpSet plot
- **Figure5.R** - Intensity distribution boxplots
- **Figure6-9.R** - Extended database and functional analysis

### Visualization Scripts (`scripts/visualization/`)
- **database_coverage_plots.R** - Individual database coverage analysis
- **database_detailed_comparison.R** - Manuscript-specific database statistics
- **variation_distribution_plot.R** - PEA vs MS/MS variability analysis
- **ranked_abundance_plots.R** - Ranked protein abundance visualizations

### Data Processing (`scripts/data_processing/`)
- **ID mapping scripts** - Protein identifier mapping and standardization
- **Data cleaning utilities** - Deduplication and quality control

## 📈 Key Datasets

### Primary Databases (`data/raw/`)
- **PeptideAtlas**: `PeptideAtlas.csv` - 4,608 canonical proteins in plasma/serum
- **Human Protein Atlas**: `HPA_MS.csv`, `HPA_PEA.csv`, `HPA_Immunoassay.csv`
- **GPMDB**: `GPMDB_Plasma.csv`, `GPMDB_Erythrocyte.csv`, `GPMDB_Platelet.csv`
- **PaxDB**: `PaxDb_Plasma.csv`, `PaxDb_CD4.csv`, `PaxDb_CD8.csv`, etc.

### Processed Data
- **Cell type datasets**: Blood cell-specific protein expression data
- **Biomarker lists**: Curated biomarker protein sets
- **Correlation matrices**: Cross-database protein correlation data

## 🎯 Key Features

### Database Integration
- **4 major databases** integrated: PeptideAtlas, HPA, GPMDB, PaxDB
- **12,713 unique genes** mapped across blood components
- **355x improvement** in PaxDB mapping (20 → 7,033 genes)
- **96% ENSP coverage** achieved

### Analysis Scope
- **Blood taxonomy**: Plasma (55%), erythrocytes (43%), leukocytes/platelets (2%)
- **11 blood cell types** analyzed
- **Protein variability**: PEA vs MS/MS coefficient of variation analysis
- **Publication-ready visualizations**: 25+ high-quality plots

### Technical Achievements
- **Dynamic range**: 10-12 orders of magnitude in protein concentrations
- **Technology comparison**: MS proteomics vs affinity-based methods
- **Cross-platform validation**: Correlation analysis across databases

## 📋 Requirements

### R Dependencies
```r
# Core packages (automatically installed)
ggplot2, dplyr, tidyr, viridis, RColorBrewer
scales, gridExtra, stringr, cowplot, ggrepel
UpSetR, corrplot, ggpubr, ComplexHeatmap
```

### System Requirements
- **R version**: ≥ 4.0.0
- **Memory**: ≥ 8GB RAM (for large datasets)
- **Storage**: ≥ 2GB free space for outputs

## 🔧 Customization

### Running Individual Scripts
```bash
# Run specific analysis from project root
Rscript scripts/analysis/Figure2.R
Rscript scripts/visualization/variation_distribution_plot.R
```

### Modifying Output Locations
Edit the utility functions in `scripts/utilities/load_packages.R`:
```r
# Customize output paths
get_output_path("my_plot.png", "plots")    # outputs/plots/my_plot.png
get_data_path("my_data.csv", "raw")        # data/raw/my_data.csv
```

## 📊 Output Organization

### Plot Categories (`outputs/plots/`)
- **01_Database_Analysis** - Database coverage and correlations
- **02_Coverage_Analysis** - Cell type data completeness
- **03_Abundance_Analysis** - Protein concentration distributions
- **04_Cell_Type_Comparison** - Comparative cell type analysis
- **05_Functional_Analysis** - Protein functional categories

### Key Visualizations
1. **Database Coverage**: Individual database protein counts and overlaps
2. **Variability Analysis**: PEA vs MS/MS coefficient of variation comparison
3. **Abundance Distributions**: Protein concentration ranges across platforms
4. **Cell Type Profiles**: Blood cell-specific protein expression patterns

## 🐛 Troubleshooting

### Common Issues
1. **Path errors**: Ensure you're running from the project root directory
2. **Missing packages**: Run `Rscript config/install_dependencies.R`
3. **Memory issues**: Close other applications when processing large datasets
4. **Permission errors**: Check write permissions for `outputs/` directory

### Log Files
Check `outputs/logs/` for detailed error messages:
```bash
# View recent log
ls -la outputs/logs/
tail -20 outputs/logs/[script_name]_[timestamp].log
```

## 📚 Citation

If you use this analysis pipeline, please cite:
```
[Your manuscript citation when published]
```

## 🤝 Contributing

1. Fork the repository
2. Create feature branch: `git checkout -b feature-name`
3. Add scripts to appropriate `scripts/` subdirectory
4. Update documentation
5. Submit pull request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙋 Support

For questions or issues:
1. Check the troubleshooting section above
2. Review log files in `outputs/logs/`
3. Open an issue on GitHub with error details

---

**Last updated**: December 2024  
**Pipeline version**: 2.0 (Organized Structure)  
**R compatibility**: R ≥ 4.0.0 