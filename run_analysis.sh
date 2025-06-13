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
while [[ $# -gt 0 ]]; do
    case $1 in
        --force-mapping)
            FORCE_MAPPING="--force-mapping"
            shift
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Print header
echo "==============================================================================="
echo "                    BLOOD REVIEW DATA ANALYSIS"
echo "                    Plasma Protein Quantification"
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
print_status "Running plasma protein analysis..."
echo ""
echo "Analyzing protein quantification across:"
echo "  • PeptideAtlas (MS)"
echo "  • HPA MS (Mass Spectrometry)"
echo "  • HPA PEA (Proximity Extension Assay)"
echo "  • HPA Immunoassay"
echo "  • GPMDB (MS)"
echo "  • PAXDB (Expression data)"
echo ""

if [ -n "$FORCE_MAPPING" ]; then
    if Rscript scripts/01_plasma_protein_analysis.R "$FORCE_MAPPING"; then
        print_success "Plasma protein analysis completed successfully!"
    else
        print_error "Plasma protein analysis failed!"
        exit 1
    fi
else
    if Rscript scripts/01_plasma_protein_analysis.R; then
        print_success "Plasma protein analysis completed successfully!"
    else
        print_error "Plasma protein analysis failed!"
        exit 1
    fi
fi

# Run the PeptideAtlas quantification analysis
print_status "Running PeptideAtlas quantification comparison analysis..."
echo ""
echo "Comparing PeptideAtlas quantification methods:"
echo "  • n_observations vs norm_PSMs_per_100K"
echo "  • Distribution analysis and correlation"
echo ""

if [ -n "$FORCE_MAPPING" ]; then
    if Rscript scripts/02_peptideatlas_quantification_analysis.R "$FORCE_MAPPING"; then
        print_success "PeptideAtlas quantification analysis completed successfully!"
    else
        print_error "PeptideAtlas quantification analysis failed!"
        exit 1
    fi
else
    if Rscript scripts/02_peptideatlas_quantification_analysis.R; then
        print_success "PeptideAtlas quantification analysis completed successfully!"
    else
        print_error "PeptideAtlas quantification analysis failed!"
        exit 1
    fi
fi

# Run the biomarker plasma analysis
print_status "Running biomarker plasma analysis..."
if [ -n "$FORCE_MAPPING" ]; then
    if Rscript scripts/biomarker_plasma_analysis.R "$FORCE_MAPPING"; then
        print_success "Biomarker plasma analysis completed successfully!"
    else
        print_error "Biomarker plasma analysis failed!"
        exit 1
    fi
else
    if Rscript scripts/biomarker_plasma_analysis.R; then
        print_success "Biomarker plasma analysis completed successfully!"
    else
        print_error "Biomarker plasma analysis failed!"
        exit 1
    fi
fi

# Display results
echo ""
echo "==============================================================================="
echo "                           ANALYSIS COMPLETE"
echo "==============================================================================="
echo ""
print_success "Generated files:"
echo "  📊 Plots:"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_by_source.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_by_technology.png"
echo "     • outputs/plots/01_plasma_protein_analysis/plasma_proteins_comprehensive.png"
echo "     • outputs/plots/02_peptideatlas_quantification_analysis/quantification_methods_*.png"
echo "     • outputs/plots/02_peptideatlas_quantification_analysis/peptideatlas_quantification_comprehensive.png"
echo "     • outputs/plots/biomarker_plasma_analysis/expression_distribution_*.png"
echo "     • outputs/plots/biomarker_plasma_analysis/combined_biomarker_distributions.png"
echo ""
echo "  📋 Data & Reports:"
echo "     • outputs/plasma_protein_counts_summary.csv"
echo "     • outputs/analysis_summary.txt"
echo "     • outputs/peptideatlas_quantification/quantification_analysis_summary.txt"
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
echo "View the generated plots and reports for detailed insights into"
echo "plasma protein quantification across different data sources and technologies." 