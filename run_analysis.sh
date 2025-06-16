#!/bin/bash

#===============================================================================
# Blood Review Data Analysis Runner
# 
# Purpose: Execute plasma protein analysis across multiple data sources
# Script: run_analysis.sh
# 
# This script runs the comprehensive analysis of protein quantification 
# in plasma from PeptideAtlas, HPA (MS/PEA/Immunoassay), GPMDB, and PAXDB
#===============================================================================

# Set script to exit on any error
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if command exists
check_command() {
    if ! command -v $1 &> /dev/null; then
        print_error "$1 could not be found. Please install $1 first."
        exit 1
    fi
}

# Parse command line arguments
FORCE_MAPPING=""
PLOT_FORMATS="svg"  # Default to SVG only
while [[ $# -gt 0 ]]; do
    case $1 in
        --force-mapping)
            FORCE_MAPPING="--force-mapping"
            shift
            ;;
        --formats)
            PLOT_FORMATS="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [--force-mapping] [--formats FORMAT_LIST]"
            echo ""
            echo "Options:"
            echo "  --force-mapping  Force re-mapping of protein IDs to gene symbols"
            echo "  --formats LIST   Comma-separated list of output formats (default: svg)"
            echo "                   Supported formats: svg, png, pdf, jpeg, tiff"
            echo "                   Example: svg,png or png,svg,pdf"
            echo "  --help, -h       Show this help message"
            echo ""
            echo "Example usage:"
            echo "  $0                          # Run with default settings (SVG output)"
            echo "  $0 --formats svg,png        # Generate both SVG and PNG formats"
            echo "  $0 --force-mapping --formats png  # Force mapping and PNG output"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            echo "Usage: $0 [--force-mapping] [--formats FORMAT_LIST]"
            echo "  --force-mapping: Force re-mapping of protein IDs to gene symbols"
            echo "  --formats:       Comma-separated list of output formats (default: svg)"
            echo "                   Example: svg,png or png,svg,pdf"
            echo "Use --help for more information."
            exit 1
            ;;
    esac
done

# Print header
echo "==============================================================================="
echo "                    BLOOD REVIEW DATA ANALYSIS"
echo "              Plasma and Serum Protein Quantification"
echo "==============================================================================="
echo ""

# Check prerequisites
print_status "Checking prerequisites..."
check_command "Rscript"

# Check if required directories exist
if [ ! -d "data/raw" ]; then
    print_error "Data directory 'data/raw' not found. Please ensure data is properly organized."
    exit 1
fi

# Check if required data files exist
required_files=(
    "data/raw/peptideatlas/peptideatlas.csv"
    "data/raw/hpa/hpa_ms.csv"
    "data/raw/hpa/hpa_pea.csv"
    "data/raw/hpa/hpa_immunoassay_plasma.csv"
    "data/raw/gpmdb/gpmdb_plasma.csv"
    "data/raw/paxdb/paxdb_plasma.csv"
    "data/raw/gpmdb/gpmdb_serum.csv"
    "data/raw/paxdb/paxdb_serum.csv"
    "data/raw/hpa/hpa_immunoassay_serum.csv"
    "data/raw/paxdb/paxdb_b_cell.csv"
    "data/raw/paxdb/paxdb_cd4.csv"
    "data/raw/paxdb/paxdb_cd8.csv"
    "data/raw/paxdb/paxdb_monocyte.csv"
    "data/raw/paxdb/paxdb_nk.csv"
    "data/raw/proteomexchange/pxd004352.csv"
    "data/raw/proteomexchange/pxd025174.csv"
    "data/raw/proteomexchange/pxd040957_cd8.csv"
    "data/raw/proteomexchange/pxd040957_macrophages.csv"
)

print_status "Verifying required data files..."
missing_files=()
for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        missing_files+=("$file")
    fi
done

if [ ${#missing_files[@]} -ne 0 ]; then
    print_error "Missing required data files:"
    for file in "${missing_files[@]}"; do
        echo "  - $file"
    done
    exit 1
fi

print_success "All required data files found."

# Create output directories
print_status "Setting up output directories..."
mkdir -p outputs
mkdir -p plots

# Run the plasma protein analysis
print_status "Running plasma protein analysis (modular pipeline)..."
echo ""
echo "🚀 Using modular analysis pipeline with:"
echo "  • Centralized configuration"
echo "  • Unified data loading"
echo "  • Standardized plotting themes"
echo "  • Enhanced normalization comparisons"
echo "  • Plot formats: $PLOT_FORMATS"
echo ""
echo "Analyzing protein quantification across:"
echo "  • PeptideAtlas (MS)"
echo "  • HPA MS (Mass Spectrometry)"
echo "  • HPA PEA (Proximity Extension Assay)"
echo "  • HPA Immunoassay"
echo "  • GPMDB (MS)"
echo "  • PAXDB (Expression data)"
echo ""

# Build arguments for R scripts
R_ARGS=""
if [ -n "$FORCE_MAPPING" ]; then
    R_ARGS="$R_ARGS $FORCE_MAPPING"
fi
R_ARGS="$R_ARGS --formats $PLOT_FORMATS"

if Rscript scripts/01_plasma_protein_analysis.R $R_ARGS; then
    print_success "Plasma protein analysis completed successfully!"
else
    print_error "Plasma protein analysis failed!"
    exit 1
fi

# Run the PeptideAtlas quantification analysis
print_status "Running PeptideAtlas quantification comparison analysis..."
echo ""
echo "Comparing PeptideAtlas quantification methods:"
echo "  • n_observations vs norm_PSMs_per_100K"
echo "  • Distribution analysis and correlation"
echo ""

if Rscript scripts/02_peptideatlas_quantification_analysis.R $R_ARGS; then
    print_success "PeptideAtlas quantification analysis completed successfully!"
else
    print_error "PeptideAtlas quantification analysis failed!"
    exit 1
fi

# Run the biomarker plasma analysis
print_status "Running biomarker plasma analysis..."
if Rscript scripts/03_biomarker_plasma_analysis.R $R_ARGS; then
    print_success "Biomarker plasma analysis completed successfully!"
else
    print_error "Biomarker plasma analysis failed!"
    exit 1
fi

# Run the serum protein analysis
print_status "Running serum protein analysis..."
echo ""
echo "Analyzing serum protein quantification across:"
echo "  • GPMDB (MS)"
echo "  • PAXDB (Expression data)"
echo "  • HPA Immunoassay"
echo ""

if Rscript scripts/04_serum_protein_analysis.R $R_ARGS; then
    print_success "Serum protein analysis completed successfully!"
else
    print_error "Serum protein analysis failed!"
    exit 1
fi

# Run the cell type analysis
print_status "Running cell type protein expression analysis..."
echo ""
echo "Analyzing protein expression across blood cell types:"
echo "  • PAXDB: B cells, CD4/CD8 T cells, NK cells, monocytes, plasma, serum"
echo "  • ProteomeXchange pxd004352: 21 immune cell subtypes"
echo "  • ProteomeXchange pxd025174: CD4/CD8 T cell copy numbers"
echo "  • ProteomeXchange pxd040957: CD8 T cells and macrophages"
echo ""

if Rscript scripts/05_celltype_analysis.R $R_ARGS; then
    print_success "Cell type analysis completed successfully!"
else
    print_error "Cell type analysis failed!"
    exit 1
fi


# Display results
echo ""
echo "==============================================================================="
echo "                           ANALYSIS COMPLETE"
echo "==============================================================================="
echo ""
print_success "Generated files:"
echo "  📊 Plots:"
echo "     🚀 Analysis Outputs:"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_by_source.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_by_technology.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_overview.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_abundance_distributions.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_zscore_distributions.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_quantification_distributions_log10.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_quantification_distributions_zscore.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_quantification_density_log10.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_quantification_density_zscore.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_log10_vs_zscore_comparison.png"
echo "     📋 Data Tables:"
echo "     • outputs/tables/01_plasma_protein_analysis/plasma_protein_summary.csv"
echo "     • outputs/tables/01_plasma_protein_analysis/plasma_proteins_normalization_statistics.csv"
echo "     📖 Reports:"
echo "     • outputs/reports/01_plasma_protein_analysis/plasma_protein_analysis_report.md"
echo "     • outputs/plots/02_peptideatlas_quantification_analysis/correlation_plot.png"
echo "     • outputs/plots/02_peptideatlas_quantification_analysis/distribution_plot.png"
echo "     • outputs/plots/02_peptideatlas_quantification_analysis/comprehensive_plot.png"
echo "     • outputs/plots/03_biomarker_plasma_analysis/*_distribution.png"
echo "     • outputs/plots/03_biomarker_plasma_analysis/*_dynamic_range.png"
echo "     • outputs/plots/03_biomarker_plasma_analysis/comprehensive_biomarker_analysis.png"
echo "     • outputs/plots/03_biomarker_plasma_analysis/biomarker_technology_overlap.png"
echo "     • outputs/plots/04_serum_protein_analysis/serum_protein_counts_by_source.png"
echo "     • outputs/plots/04_serum_protein_analysis/serum_protein_counts_combined.png"
echo "     • outputs/plots/04_serum_protein_analysis/serum_protein_overlap_upset.png"
echo "     • outputs/plots/04_serum_protein_analysis/serum_protein_quantification_distributions.png"
echo "     • outputs/plots/04_serum_protein_analysis/serum_protein_comprehensive_summary.png"
echo "     • outputs/plots/05_celltype_analysis/celltype_gene_counts.png"
echo "     • outputs/plots/05_celltype_analysis/data_source_coverage.png"
echo "     • outputs/plots/05_celltype_analysis/celltype_source_matrix.png"
echo "     • outputs/plots/05_celltype_analysis/intensity_distributions.png"
echo "     • outputs/plots/05_celltype_analysis/technology_comparison.png"
echo "     • outputs/plots/05_celltype_analysis/comprehensive_summary.png"
echo "     • outputs/plots/07_plasma_databases_dot_plot/plasma_databases_comparison.svg"
echo ""
echo "  📋 Data & Reports:"
echo "     • outputs/plasma_protein_counts_summary.csv"
echo "     • outputs/03_biomarker_plasma_analysis/biomarker_detection_summary.csv"
echo "     • outputs/serum_protein/serum_protein_counts_summary.csv"
echo "     • outputs/serum_protein/serum_protein_overlap_statistics.csv"
echo "     • outputs/celltype_analysis/celltype_protein_data.csv"
echo "     • outputs/celltype_analysis/celltype_overall_summary.csv"
echo "     • outputs/celltype_analysis/celltype_summary_by_source.csv"
echo ""
echo "  📖 Documentation:"
echo "     • scripts/plasma_protein_analysis_summary.md"
echo ""

# Display quick summary if analysis summary exists
if [ -f "outputs/analysis_summary.txt" ]; then
    print_status "Quick Summary:"
    echo ""
    grep -A 10 "SUMMARY STATISTICS:" outputs/plasma_protein/analysis_summary.txt | head -6
    echo ""
    echo "Technology Classification:"
    echo "• MS-based: PeptideAtlas, HPA MS, GPMDB, PAXDB"
    echo "• PEA: HPA PEA (Proximity Extension Assay)"
    echo "• Immunoassay: HPA Immunoassay"
    echo ""
    echo "PeptideAtlas Quantification Recommendation: norm_PSMs_per_100K"
    echo ""
fi

print_success "Analysis pipeline completed successfully!"
echo ""

echo "🚀 ANALYSIS FEATURES:"
echo "• Modular code architecture with centralized configuration"
echo "• Enhanced normalization comparisons for better cross-database analysis"
echo "• Standardized plotting themes and consistent styling"
echo "• Automatic report generation with comprehensive statistics"
echo "• Improved error handling and data validation"
echo ""

echo "View the generated plots and reports for detailed insights into:"
echo "• Plasma and serum protein quantification across different data sources and technologies"
echo "• Enhanced normalization comparisons for better cross-database analysis"
echo "• Comprehensive statistical summaries of normalization methods"
echo "• Blood cell type protein expression profiles (26 cell types, 14,274 unique genes)" 