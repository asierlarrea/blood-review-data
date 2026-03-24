# PeptideAtlas Quantification Methods Comparison Analysis

**Analysis Date:** 2025-12-07
**Script:** `02_peptideatlas_quantification_analysis.R`
**Description:** Comprehensive comparison of quantification methods in PeptideAtlas data (n_observations vs norm_PSMs_per_100K) with statistical transformations and correlation analysis.

---

## Summary Statistics

| Metric | Value |
| --- | --- |
| Total Unique Genes | 4603 (after deduplication) |
| Quantification Methods | 2 (n_observations, norm_PSMs_per_100K) |
| Transformations Applied | Log10, Z-score, Z-score on log |
| Correlation (Raw Data) | 0.707 |
| Correlation (Log-transformed) | 0.966 |
| Correlation (Z-score normalized) | 0.966 |

## Methodology

- **Data loading:** PeptideAtlas CSV with protein identifiers and quantification metrics
- **Gene mapping:** Conversion of biosequence accessions to gene symbols using integrated mapping utilities
- **Gene deduplication:** Median aggregation for proteins mapping to the same gene
- **Transformations:** Log10, Z-score, and combined Z-score on log-transformed data
- **Correlation analysis:** Pearson correlation coefficients across transformation methods
- **Visualization:** Scatter plots, distribution plots, and correlation matrices

## Recommendations

- **Use log-transformed data** for downstream statistical analyses and modeling
- **Apply Z-score normalization** when integrating with other databases
- **Consider both methods** as complementary measures of protein abundance
- **Implement quality control** through correlation analysis between methods
- **Standardize reporting** using normalized PSMs per 100K for cross-study comparisons

## Generated Files

- **Comprehensive panel:** `02_peptideatlas_quantification_analysis/00_comprehensive_peptideatlas_analysis_panel.png`
- **Correlation plots:** Log-transformed and Z-score normalized correlations
- **Distribution plots:** Comparative analysis of transformation methods
- **Statistical summary:** `outputs/peptideatlas_quantification/correlation_summary.csv`
- **Processed data:** `data/processed/peptideatlas_mapped_genes.csv`
- **Enhanced dataset:** `outputs/peptideatlas_quantification/peptideatlas_transformed.csv`

---
*Report generated automatically by the blood proteomics analysis pipeline*

