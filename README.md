# Blood Proteomics Analysis Pipeline

Analysis pipeline for integrating six proteomics databases into a unified, gene-centric atlas of the human blood proteome. Covers plasma, serum, and blood cell compartments with cross-database normalization and biomarker detection analysis.

## Overview

| Compartment | Databases | Unique Genes | Analysis Script |
|---|---|---|---|
| Plasma | 6 (PeptideAtlas, PAXDB, HPA MS, HPA PEA, GPMDB, QuantMS) | ~8,600 | `01_plasma_protein_analysis.R` |
| Serum | 3 (GPMDB, PAXDB, HPA Immunoassay) | ~4,100 | `04_serum_protein_analysis.R` |
| Blood cells | 3 (PAXDB, ProteomeXchange, GPMDB) | >13,000 across 11 cell types | `05_celltype_analysis.R` |
| Biomarkers | 129 curated markers tested across all DBs | -- | `03_biomarker_plasma_analysis.R` |

## Quick Start

```bash
# 1. Install R dependencies
Rscript install_dependencies.R

# 2. Run full pipeline
Rscript run_analysis.R

# 3. Or run individual scripts
Rscript scripts/01_plasma_protein_analysis.R
```

Results go to `outputs/plots/`, `outputs/tables/`, and `outputs/reports/`.

## Repository Structure

```
scripts/
  01_plasma_protein_analysis.R        # 6-panel plasma proteome figure
  02_peptideatlas_quantification_analysis.R  # PeptideAtlas method comparison
  03_biomarker_plasma_analysis.R      # Biomarker detection across DBs
  04_serum_protein_analysis.R         # Serum proteome analysis
  05_celltype_analysis.R              # Blood cell-type proteomics
  06_final_comprehensive_figure.R     # Combined publication figure
  config/analysis_config.R            # Central parameters and paths
  utilities/                          # Shared functions (loading, themes, normalization)
  data_processing/                    # ID mapping and data cleaning

data/
  raw/                                # Source data by database
    peptideatlas/, hpa/, paxdb/, gpmdb/, quantms/, proteomexchange/
  metadata/                           # Biomarker lists
  processed/                          # Intermediate gene-mapped files
  cache/                              # Protein-to-gene mapping cache

outputs/
  plots/                              # Publication-ready figures (PNG, TIFF, SVG)
  tables/                             # Summary statistics per analysis
  reports/                            # Auto-generated markdown reports
```

## Data Sources

| Database | Technology | Proteins | Metric |
|---|---|---|---|
| PeptideAtlas | MS | 4,603 | PSMs/100K |
| PAXDB | MS | 7,021 | ppm |
| HPA MS | MS | 4,294 | mg/L |
| HPA PEA | PEA | 1,436 | NPX |
| GPMDB | MS | 2,266 | spectral counts |
| QuantMS | MS (reanalysis) | 2,799 | iBAQ |

All abundance values are normalized to a common scale using quantile-to-normal transformation for cross-database comparison.

## Availability of Data and Materials

- **PeptideAtlas**: Deutsch EW, Omenn GS, Sun Z, Maes M, Pernemalm M, Palaniappan KK, Letunica N, Vandenbrouck Y, Brun V, Tao SC, et al: Advances and Utility of the Human Plasma Proteome. https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProteins?atlas_build_id=603&organism_id=2&redundancy_constraint=4&presence_level_constraint=1&action=QUERY (2025).
- **Human Protein Atlas**: Alvez MB, Bergstrom S, Kenrick J, Johansson E, Aberg M, Akyildiz M, Altay O, Skold H, Antonopoulos K, Apostolakis E, et al: A human pan-disease blood atlas of the circulating proteome. https://www.proteinatlas.org/about/download (2025).
- **PaxDb**: Huang Q, Szklarczyk D, Oehninger J, von Mering C: PaxDb v6.0: reprocessed, LLM-selected, curated protein abundance data across organisms. https://pax-db.org/species/9606 (2025).
- **GPMDB**: Craig R, Cortens JC, Fenyö D, Beavis RC. Using annotated peptide mass spectrum libraries for protein identification. https://www.thegpm.org/lists/index.html#201507081 (2025).
- **quantms**: Dai C, Pfeuffer J, Wang H, Zheng P, Kall L, Sachsenberg T, Demichev V, Bai M, Kohlbacher O, Perez-Riverol Y: quantms: a cloud-based pipeline for quantitative proteomics enables the reanalysis of public proteomics data. https://quantms.org/datasets (2025).
- **Blood Proteoform Atlas**: Melani RD, Gerbasi VR, Anderson LC, Sikora JW, Toby TK, Hutton JE, Butcher DS, Negrao F, Seckler HS, Srzentic K, et al: The Blood Proteoform Atlas: A reference map of proteoforms in human hematopoietic cells. https://blood-proteoform-atlas.org/proteins (2025).
- **PXD040957**: Canale F, Neumann J, Renesse J et al. Proteomics of immune cells from liver tumors reveals immunotherapy targets. PRIDE. 10.1016/j.xgen.2023.100331 (2023).
- **PXD004352**: Rieckmann J, Geiger R, Hornburg D et al. Social network architecture of human immune cells unveiled by quantitative proteomics. PRIDE. 10.1038/ni.3693 (2016).

## Requirements

**R packages** (installed automatically via `install_dependencies.R`):
- Core: ggplot2, dplyr, tidyr, readr, stringr, scales
- Visualization: patchwork, ggpubr, UpSetR, ggupset, ggridges, viridis

## License

MIT License -- see [LICENSE](LICENSE).
