# Blood Proteomics Analysis – Project Summary

**Date:** 2025-06-16  
**Repository:** blood-review-data  
**Scope:** Plasma, serum, circulating cell-type and biomarker analyses (quantile-integration module intentionally excluded)

---

## 1. Study Rationale
Human blood hosts a diverse proteome reflecting systemic physiology and disease. By horizontally integrating six public databases and multiple targeted datasets we provide a harmonised, gene-centric atlas suitable for biomarker discovery and method benchmarking.

---

## 2. Data Landscape
| Matrix | Sources (n) | Technologies | Unique protein-coding genes |
|--------|-------------|--------------|----------------------------|
| Plasma | 6 | MS, PEA, Immunoassay | **8 374** |
| Serum  | 3 | MS, Immunoassay | **4 163** |
| Circulating cells | 11 cell types, 4 main repositories | MS | **13 000+** |
| Biomarker panel | MarkerDB (129 proteins) | – | – |

*See individual CSV tables under `outputs/tables/` for granular statistics.*

---

## 3. Plasma Protein Findings
* `peptideatlas` contributed the largest single set (4 603 genes).
* Affinity platforms (HPA-PEA, Immunoassay) captured low-abundance proteins absent in MS data but covered only ~41 % of the biomarker panel.
* Combined plasma union encompasses ~42 % more genes than any individual source (Fig. P1).

---

## 4. Serum Protein Findings
* `GPMDB` dominated with 3 576 unique genes; `PAXDB` added 823.
* Immunoassay (HPA) provided 336 high-confidence concentration estimates.
* Total serum catalogue: 4 163 genes (Table S1, Fig. S1).

---

## 5. Immune & Erythroid Cell-Type Atlas
Top gene counts:

| Cell type | Genes | Primary datasets |
|-----------|-------|------------------|
| CD8 T cells | **11 379** | PXDs 004352, 025174, 040957; PAXDB |
| B cells | 9 802 | PXD 004352; PAXDB |
| NK cells | 9 481 | PXD 004352; PAXDB |
| CD4 T cells | 8 582 | PXDs 004352, 025174; PAXDB |
| Thrombocytes | 8 495 | PXD 004352; GPMDB |

Complete summaries are in `outputs/celltype_analysis/`.

---

## 6. Clinical Biomarker Coverage
Detection of 129 curated plasma biomarkers:

| Repository | Biomarkers detected | % of panel |
|------------|--------------------|-----------|
| PeptideAtlas | 105 | **81 %** |
| HPA-MS | 101 | 78 % |
| PAXDB | 101 | 78 % |
| GPMDB | 79 | 61 % |
| HPA-PEA | 53 | 41 % |
| HPA Immunoassay | 53 | 41 % |

Waterfall and heat-map visualisations are available in `outputs/plots/03_biomarker_plasma_analysis/`.

---

## 7. Methodological Notes
1. **Identifier harmonisation:** UniProt & Ensembl translation to HGNC gene symbols (`add_gene_names.py`).
2. **Gene-level aggregation:** Median of protein groups per gene.
3. **Transformations:** Source-specific log10 or linear scaling as appropriate.
4. **Quality control:** Removal of missing/zero abundances, duplicate removal (`utilities/gene_deduplication.R`).

*Note:* Cross-database quantile normalisation analysis was removed per project decision; therefore, no integrated abundance matrix is presented here.

---

## 8. Key Takeaways
* The union of six plasma repositories increases coverage by ~2 800 genes relative to the largest single source.
* Serum proteome is less diverse than plasma but still adds 400+ unique genes not observed in plasma.
* Cell-type datasets reveal thousands of intracellular proteins absent from fluid samples—valuable for cell-specific marker discovery.
* PeptideAtlas is the most comprehensive single resource for clinically validated biomarkers.

---

## 9. Repository Structure Reference
```
outputs/
 ├─ tables/         # CSV summaries per analysis
 ├─ plots/          # PNG/PDF visualisations
 ├─ reports/        # Script-level markdown reports
scripts/
 ├─ 01_* – 05_*     # Analysis pipelines (R)
 └─ utilities/      # Helper functions
```

---

## 10. How to Reproduce
```
# Install R dependencies
Rscript install_dependencies.R

# Run the entire workflow (takes ~30 min)
bash run_analysis.sh
```
Individual analyses can be executed via the numbered R scripts.

---

## 11. Future Directions
* Integrate post-translational modifications (PTMs).
* Develop an interactive Shiny dashboard for browsing the atlas.
* Incorporate longitudinal clinical cohorts for biomarker validation.

---

*Please cite this summary if you use the data or code.* 