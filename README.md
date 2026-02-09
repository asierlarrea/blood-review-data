# Blood Proteomics Analysis Pipeline

A comprehensive analysis pipeline for blood proteomics data integration and visualization, focusing on cross-database comparison, biomarker detection, and quantile-normalized abundance analysis across plasma, serum, and cellular components of blood.

## 🔬 Analysis Overview

This pipeline integrates **6 major proteomics databases** to create a unified, gene-centric atlas of the human blood proteome, analyzing:
- **🩸 Plasma proteins**: 8,665 unique genes across 6 data sources
- **🧪 Serum proteins**: 4,163 unique genes across 3 data sources  
- **🔬 Cell-specific proteins**: >13,000 proteins from 11 blood cell types
- **🎯 Clinical biomarkers**: 129 curated biomarkers with detection analysis

## 🗂️ Repository Structure

```
blood-review-data/
├── 📁 data/                          # Organized data files
│   ├── raw/                         # Original source data
│   │   ├── gpmdb/                   # GPMDB spectral count data
│   │   ├── hpa/                     # Human Protein Atlas (MS, PEA, Immunoassay)
│   │   ├── paxdb/                   # PaxDB abundance data (plasma, serum, cells)
│   │   ├── peptideatlas/            # PeptideAtlas canonical proteins
│   │   ├── quantms/                 # QuantMS reprocessed data
│   │   ├── proteomexchange/         # ProteomeXchange cell-specific data
│   │   └── markerdb_biomarker/      # Clinical biomarker annotations
│   ├── processed/                   # Cleaned and mapped data
│   ├── metadata/                    # Biomarker lists and documentation
│   └── cache/                       # ID mapping cache for performance
├── 📁 scripts/                       # Analysis pipeline
│   ├── 01_plasma_protein_analysis.R           # Main plasma proteome analysis
│   ├── 02_peptideatlas_quantification_analysis.R  # PeptideAtlas method comparison
│   ├── 03_biomarker_plasma_analysis.R         # Clinical biomarker detection
│   ├── 04_serum_protein_analysis.R            # Serum proteome analysis
│   ├── 05_celltype_analysis.R                 # Blood cell proteomics
│   ├── utilities/                   # Core functions and themes
│   ├── data_processing/            # ID mapping and cleaning
│   ├── config/                     # Analysis configuration
│   └── visualization/              # Specialized plotting functions
├── 📁 outputs/                       # Generated results
│   ├── plots/                      # Publication-ready visualizations
│   │   ├── 01_plasma_protein_analysis/     # Comprehensive plasma analysis
│   │   ├── 02_peptideatlas_quantification_analysis/  # PeptideAtlas methods
│   │   ├── 03_biomarker_plasma_analysis/   # Biomarker detection plots
│   │   ├── 04_serum_protein_analysis/      # Serum analysis
│   │   └── 05_celltype_analysis/           # Cell-type specific plots
│   ├── tables/                     # Summary statistics and data tables
│   ├── reports/                    # Auto-generated analysis reports
│   └── logs/                       # Execution logs for debugging
├── 📁 manuscript/                    # Manuscript files and bibliography
├── run_analysis.R                   # Complete pipeline runner
├── install_dependencies.R          # Automated R package installation
├── PROJECT_RESULTS_SUMMARY.md      # Comprehensive analysis summary
└── README.md                       # This documentation
```

## 🚀 Quick Start

### 1. Install Dependencies
```bash
# Automatically install all required R packages
Rscript install_dependencies.R
```

### 2. Run Complete Analysis Pipeline
```bash
# Execute all analysis scripts with default settings
Rscript run_analysis.R

# With options for ID re-mapping and plot formats
Rscript run_analysis.R --force-mapping --formats "png,svg,tiff"
```

### 3. Explore Results
```bash
# Check comprehensive outputs
ls outputs/plots/*/                    # All generated visualizations
ls outputs/tables/*/                   # Summary statistics tables
cat outputs/reports/*/                 # Auto-generated analysis reports
```

## 📊 Core Analysis Components

### **Script 01: Plasma Protein Analysis** 
**Comprehensive 6-panel visualization of plasma proteome**
- **Panel A**: Database coverage comparison
- **Panel B**: Protein overlap analysis (UpSet plot)
- **Panel C**: Abundance distribution shapes  
- **Panel D**: Cross-database correlation analysis
- **Panel E**: QuantMS sample presence analysis
- **Panel F**: Complementary platform strengths

**Key Outputs**: `00_comprehensive_plasma_analysis_panel.png` (30×30 inch publication figure)

### **Script 02: PeptideAtlas Quantification Analysis**
**Comparison of PeptideAtlas quantification methods**
- Log-transformed vs Z-score normalized correlations
- Distribution analysis of n_observations vs norm_PSMs_per_100K
- Method validation and correlation assessment

### **Script 03: Biomarker Plasma Analysis** 
**Clinical biomarker detection across databases**
- **Panel A**: Biomarker detection heatmap across all databases
- **Panel B**: Database intersection analysis (UpSet plot)  
- **Panel C**: Biomarker abundance distributions by technology
- Detection rates: PeptideAtlas (81.4%), HPA MS (78.3%), PAXDB (78.3%)

### **Script 04: Serum Protein Analysis**
**Serum-specific proteome characterization**
- Cross-database comparison (GPMDB, PAXDB, HPA Immunoassay, QuantMS)
- Technology-specific analysis and correlation assessment
- 4,163 unique genes identified

### **Script 05: Cell-Type Analysis**
**Blood cellular component proteomics**
- 11 blood cell types analyzed (CD4/CD8 T cells, B cells, NK cells, etc.)
- Integration of PAXDB, ProteomeXchange, and GPMDB cell data
- >13,000 unique proteins across all cell types

## 🔧 Key Technical Features

### **Quantile-to-Normal Normalization**
All analyses use robust quantile-to-normal transformation for cross-database comparison:
```r
# Within each database
rank_quantile = rank(log_abundance) / (n() + 1)
z_score = qnorm(rank_quantile)
```

### **Abundance vs Expression Terminology**
- **Abundance metrics**: GPMDB (spectral counts), PAXDB (ppm), HPA MS (mg/L), PeptideAtlas (PSMs/100K)
- **Expression metrics**: HPA PEA (NPX), QuantMS (iBAQ normalized)

### **Automated ID Mapping**
- Protein accessions → HGNC gene symbols via cached mapping system
- 96% mapping success rate with intelligent fallback strategies
- Cached results for performance: `data/cache/protein_to_gene_mappings.csv`

### **Gene-Level Deduplication**
- Multiple protein isoforms aggregated to gene level using median values
- Handles complex cases with multiple accessions per gene

## 📈 Major Data Sources

| Database | Technology | Plasma/Serum Proteins | Key Files |
|----------|------------|----------------------|-----------|
| **PeptideAtlas** | MS | 4,603 | `peptideatlas.csv` |
| **PAXDB** | MS | 7,021 | `paxdb_plasma.csv`, `paxdb_serum.csv` |
| **HPA MS** | MS | 4,294 | `hpa_ms.csv` |
| **HPA PEA** | PEA | 1,436 | `hpa_pea.csv` |
| **GPMDB** | MS | 2,266 | `gpmdb_plasma.csv`, `gpmdb_serum.csv` |
| **QuantMS** | MS | 2,799 | `quantms/plasma/` (multiple files) |

### **Blood Cell Types** (Script 05)
- **CD8 T cells**: 11,379 proteins (highest coverage)
- **B cells**: 9,802 proteins  
- **NK cells**: 9,481 proteins
- **Platelets**: 8,495 proteins
- Plus: CD4 T cells, monocytes, dendritic cells, erythrocytes

## Availability of Data and Materials

- **PeptideAtlas**: Deutsch EW, Omenn GS, Sun Z, Maes M, Pernemalm M, Palaniappan KK, Letunica N, Vandenbrouck Y, Brun V, Tao SC, et al: Advances and Utility of the Human Plasma Proteome. https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProteins?atlas_build_id=603&organism_id=2&redundancy_constraint=4&presence_level_constraint=1&action=QUERY (2025).
- **Human Protein Atlas**: Alvez MB, Bergstrom S, Kenrick J, Johansson E, Aberg M, Akyildiz M, Altay O, Skold H, Antonopoulos K, Apostolakis E, et al: A human pan-disease blood atlas of the circulating proteome. https://www.proteinatlas.org/about/download (2025).
- **PaxDb**: Huang Q, Szklarczyk D, Oehninger J, von Mering C: PaxDb v6.0: reprocessed, LLM-selected, curated protein abundance data across organisms. https://pax-db.org/species/9606 (2025).
- **GPMDB**: Craig R, Cortens JC, Fenyö D, Beavis RC. Using annotated peptide mass spectrum libraries for protein identification. https://www.thegpm.org/lists/index.html#201507081 (2025).
- **quantms**: Dai C, Pfeuffer J, Wang H, Zheng P, Kall L, Sachsenberg T, Demichev V, Bai M, Kohlbacher O, Perez-Riverol Y: quantms: a cloud-based pipeline for quantitative proteomics enables the reanalysis of public proteomics data. https://quantms.org/datasets (2025).
- **Blood Proteoform Atlas**: Melani RD, Gerbasi VR, Anderson LC, Sikora JW, Toby TK, Hutton JE, Butcher DS, Negrao F, Seckler HS, Srzentic K, et al: The Blood Proteoform Atlas: A reference map of proteoforms in human hematopoietic cells. https://blood-proteoform-atlas.org/proteins (2025).
- **PXD040957**: Canale F, Neumann J, Renesse J et al. Proteomics of immune cells from liver tumors reveals immunotherapy targets. PRIDE. 10.1016/j.xgen.2023.100331 (2023). 
- **PXD004352**: Rieckmann J, Geiger R, Hornburg D et al. Social network architecture of human immune cells unveiled by quantitative proteomics. PRIDE. 10.1038/ni.3693 (2016). 

## 🎯 Key Analysis Outcomes

### **Database Complementarity**
- **Total unique genes**: 8,665 across all plasma databases
- **Core overlap**: Small intersection emphasizes need for integration
- **Technology bias**: MS provides broad coverage, PEA/Immunoassay target specific proteins

### **Biomarker Detection Performance**
- **PeptideAtlas**: 81.4% biomarker coverage (best single source)
- **HPA MS & PAXDB**: 78.3% each  
- **Cross-platform validation**: 60+ biomarkers detected in multiple databases

### **Normalization Strategy Impact**
- **Traditional z-score**: Preserves distribution differences between databases
- **Quantile-to-normal**: Enables robust cross-database comparisons
- **Panel C**: Shows actual distribution shapes using traditional normalization
- **All other panels**: Use quantile-normalized values for fair comparison

## 🛠️ System Requirements

### **R Environment**
```r
# Core analysis packages
ggplot2, dplyr, tidyr, readr, stringr, scales
patchwork, ggpubr, UpSetR, ggupset, tibble

# Specialized packages  
ggridges, RColorBrewer, viridis, png
```

### **Hardware Recommendations**
- **RAM**: ≥8GB (16GB recommended for large datasets)
- **Storage**: ≥5GB free space (plots + intermediate files)
- **CPU**: Multi-core recommended for QuantMS processing

## 🔍 Quality Control Features

### **Data Validation**
- Automated missing file detection before analysis
- Sample count thresholds for QuantMS (≥10 samples per protein)
- Abundance value filtering (positive values only)

### **Reproducibility**
- All analysis parameters centralized in `scripts/config/analysis_config.R`
- Comprehensive logging to `outputs/logs/`
- Version-controlled plotting themes and color schemes

### **Error Handling**
- Graceful handling of missing data sources
- Intelligent fallbacks for mapping failures
- Detailed error reporting with troubleshooting guidance

## 📊 Publication-Ready Outputs

### **Comprehensive Panels**
- **Script 01**: 30×30 inch 6-panel plasma analysis figure
- **Script 03**: 28×22 inch 3-panel biomarker analysis figure
- All panels optimized for publication at 600 DPI

### **Format Support**
- **PNG**: High-resolution for presentations
- **TIFF**: Publication quality with LZW compression  
- **SVG**: Vector graphics for manuscripts (when applicable)

### **Automated Reports**
- Markdown reports with analysis statistics
- Summary tables with cross-database comparisons
- Key findings and methodology documentation

## 🔧 Customization Guide

### **Running Individual Scripts**
```bash
# Run specific analysis with custom parameters
Rscript scripts/01_plasma_protein_analysis.R --force-mapping
Rscript scripts/03_biomarker_plasma_analysis.R --formats "png,tiff"
```

### **Modifying Analysis Parameters**
Edit `scripts/config/analysis_config.R`:
```r
# Customize thresholds
QUANTMS_MIN_SAMPLES <- 10        # Minimum sample presence
PLOT_DPI <- 600                  # Output resolution
DEFAULT_PLOT_WIDTH <- 20         # Figure width (inches)
```

### **Custom Color Schemes**
Modify `scripts/utilities/plot_themes.R`:
```r
# Technology colors
tech_colors <- c("MS" = "#2E86AB", "PEA" = "#F18F01", "Immunoassay" = "#A23B72")
```

## 📚 Citation & Licensing

### **Citation**
If you use this analysis pipeline, please cite:
```bibtex
@software{blood_proteomics_pipeline,
  title = {Blood Proteomics Analysis Pipeline},
  author = {Pérez, Y.},
  year = {2025},
  url = {https://github.com/yourusername/blood-review-data}
}
```

### **License**
This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## 🐛 Troubleshooting

### **Common Issues**
```bash
# Missing packages
Rscript install_dependencies.R

# Permission errors
chmod +x run_analysis.R
sudo chown -R $USER:$USER outputs/

# Memory issues for large datasets
# Edit scripts to process data in chunks
# Or increase available memory: ulimit -v unlimited
```

### **Debug Information**
```bash
# Check execution logs
tail -f outputs/logs/01_plasma_protein_analysis_*.log

# Verify data file integrity
ls -la data/raw/*/

# Test individual components
Rscript -e "source('scripts/utilities/data_loader.R'); test_data_loading()"
```

## 🤝 Contributing

1. **Fork** the repository
2. **Create** feature branch: `git checkout -b feature-analysis-improvement`
3. **Add scripts** to appropriate `scripts/` subdirectory
4. **Update documentation** and add tests
5. **Submit** pull request with clear description

### **Development Guidelines**
- Follow existing code structure and naming conventions
- Add comprehensive documentation for new functions
- Include example usage and expected outputs
- Test with different data subsets before submitting

## 📞 Support

For questions, issues, or contributions:
- **GitHub Issues**: Report bugs and request features
- **Documentation**: Check `PROJECT_RESULTS_SUMMARY.md` for detailed analysis outcomes
- **Logs**: Review `outputs/logs/` for execution details

---

**🩸 Comprehensive blood proteomics analysis pipeline - from raw data to publication-ready visualizations** 🩸
