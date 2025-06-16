# Blood Proteomics Analysis: A Comprehensive Project Summary

**Date:** 2025-06-16  
**Repository:** `blood-review-data`  
**Lead Analyst:** Y. Pérez

---

## 1. Executive Summary
This project integrates six major public proteomics databases and multiple targeted cell-type datasets to create a harmonized, gene-centric atlas of the human blood proteome. By systematically processing **8,374 plasma proteins**, **4,163 serum proteins**, and over **13,000 proteins from 11 circulating cell types**, this work provides a foundational resource for biomarker discovery and systems biology. The analysis emphasizes data source complementarity, evaluates technology biases, and establishes the superior coverage of PeptideAtlas for clinically relevant biomarkers. All code and results are reproducible via the `run_analysis.sh` script.

---

## 2. Global Data Landscape
The study consolidates data from diverse platforms into three primary matrices.

| Matrix              | Data Sources (n) | Technologies                 | Unique Genes | Key Output File                                   |
| ------------------- | ---------------- | ---------------------------- | ------------ | ------------------------------------------------- |
| **Plasma**          | 6                | MS, PEA, Immunoassay         | **8,374**    | `outputs/tables/01.../plasma_protein_summary.csv` |
| **Serum**           | 3                | MS, Immunoassay              | **4,163**    | `outputs/tables/04.../serum_protein_summary.csv`  |
| **Circulating Cells** | 4+ repos         | MS                           | **>13,000**  | `outputs/celltype_analysis/..._summary.csv`      |
| **Biomarker Panel** | 1 (MarkerDB)     | Curated List                 | **129**      | `data/metadata/biomarkers_list.csv`               |

---

## 3. Analysis of the Plasma Proteome
The plasma analysis (`scripts/01_plasma_protein_analysis.R`) reveals significant complementarity between sources.

*   **Source Contribution:** While `PeptideAtlas` (4,603 genes) and `PAXDB` (7,021 genes) are the largest individual contributors, their union is essential for comprehensive coverage, as visualized in the bar chart (`outputs/plots/01.../01_plasma_proteins_by_source.svg`).
*   **Technology Bias:** Mass Spectrometry (MS) overwhelmingly provides the broadest coverage. Affinity-based platforms (PEA, Immunoassay) detect a smaller, distinct set of proteins, often at lower abundance, but cover only ~41% of the clinical biomarker panel.
*   **Protein Overlap:** The intersection of proteins across multiple databases is relatively small, underscoring the value of integration. This is explicitly detailed in the UpSet plot (`outputs/plots/01.../04_upset_protein_overlap.png`), which shows which proteins are unique to or shared between sources.
*   **Quantitative Normalization:** The analysis pipeline systematically evaluates four normalization strategies. Violin plots (`outputs/plots/01.../07-10_dist_violin_*.svg`) demonstrate how cross-database quantile normalization most effectively aligns distributions, a key step for reducing technical batch effects.

---

## 4. Analysis of the Serum Proteome
The serum analysis (`scripts/04_serum_protein_analysis.R`) provides a focused view on this distinct blood fraction.

*   **Dominant Source:** `GPMDB` is the primary source, identifying 3,576 unique genes.
*   **Correlation:** Cross-database correlation heatmaps (`outputs/plots/04.../04_serum_correlation_heatmap.png`) show moderate quantitative agreement between sources after log-transformation and Z-score normalization, highlighting the need for careful integration.
*   **Findings:** The total serum proteome (4,163 genes) is less diverse than plasma's, likely reflecting the removal of coagulation factors.

---

## 5. Circulating Immune & Erythroid Cell-Type Atlas
The cell-type analysis (`scripts/05_celltype_analysis.R`) constructs a deep proteomic inventory of blood's cellular components.

| Cell Type       | Unique Genes | Primary Data Sources                     |
| --------------- | ------------ | ---------------------------------------- |
| **CD8+ T cells**  | **11,379**   | ProteomeXchange (PXD...), PAXDB         |
| **B cells**       | 9,802        | PXD..., PAXDB                            |
| **NK cells**      | 9,481        | PXD..., PAXDB                            |
| **Thrombocytes**  | 8,495        | PXD..., GPMDB                            |

These datasets reveal thousands of intracellular proteins absent in plasma/serum, creating a rich resource for cell-specific biomarker discovery.

---

## 6. Clinical Biomarker Coverage
A key objective was to assess the detection of 129 curated plasma biomarkers (`scripts/03_biomarker_plasma_analysis.R`).

| Repository        | Biomarkers Detected | % of Panel |
| ----------------- | ------------------- | ---------- |
| **PeptideAtlas**    | **105**             | **81.4%**  |
| HPA MS            | 101                 | 78.3%      |
| PAXDB             | 101                 | 78.3%      |
| GPMDB             | 79                  | 61.2%      |

Waterfall plots of abundance (`outputs/plots/03.../combined_waterfall.png`) illustrate that while coverage is high in MS databases, the dynamic range of concentrations spans several orders of magnitude, posing a significant analytical challenge.

---

## 7. Key Methodological Components
1.  **Identifier Harmonization:** A Python script (`add_gene_names.py`) and R utilities map diverse protein identifiers to official HGNC gene symbols.
2.  **Gene-Level Aggregation:** Protein-level quantifications are aggregated to the gene level using the median value to create a robust, unified dataset.
3.  **Reproducibility:** The entire workflow is orchestrated by `run_analysis.sh`, ensuring full reproducibility from raw data to final plots and tables.

---

## 8. Final Conclusions & Future Directions
*   **Conclusion:** This project successfully creates a comprehensive, gene-centric blood proteome atlas. The results quantify the strengths of each major data source and provide a clear rationale for integrative approaches. `PeptideAtlas` is confirmed as the most valuable single resource for clinical biomarker studies.
*   **Future Work:**
    *   Integrate post-translational modification (PTM) data.
    *   Develop an interactive R/Shiny dashboard for data exploration.
    *   Incorporate longitudinal clinical cohorts to validate findings.

---
*This summary is auto-generated. Please cite this work if using the associated data or code.* 