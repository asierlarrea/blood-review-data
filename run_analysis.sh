#!/bin/bash

# ============================================================================
# Blood Proteomics Analysis - Script Runner
# Description: Automated execution of all R analysis scripts
# Usage: ./run_analysis.sh
# ============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
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

# Check if R is installed
if ! command -v Rscript &> /dev/null; then
    print_error "Rscript not found. Please install R first."
    exit 1
fi

# Check if we're in the right directory (look for the new structure)
if [ ! -d "scripts" ] || [ ! -d "data" ] || [ ! -d "outputs" ]; then
    print_error "Required directories (scripts/, data/, outputs/) not found."
    print_error "Please run this script from the organized project directory."
    exit 1
fi

# Create output directories if they don't exist
mkdir -p outputs/{plots,tables,reports,logs}

print_status "Starting Blood Proteomics Analysis Pipeline"
print_status "Working directory: $(pwd)"
print_status "Using organized repository structure"

# Array of scripts to run (with new paths)
scripts=(
    "scripts/analysis/figure2.R" 
    "scripts/analysis/figure3.R" 
    "scripts/analysis/figure4.R" 
    "scripts/analysis/figure5.R"
    "scripts/analysis/figure6_database_correlation.R"
    "scripts/analysis/figure7_abundance_analysis.R"
    "scripts/analysis/figure8_functional_analysis.R"
    "scripts/analysis/figure9_comparative_analysis.R"
    "scripts/visualization/database_coverage_plots.R"
    "scripts/visualization/database_detailed_comparison.R"
    "scripts/visualization/variation_distribution_plot.R"
)

descriptions=(
    "UpSet plot of plasma proteins"
    "Cell type bubble plot"
    "Cell type UpSet plot" 
    "Intensity distribution boxplots"
    "Database coverage and correlation analysis"
    "Protein abundance distribution analysis"
    "Functional categories and immunoglobulin analysis"
    "Comparative biological analysis (immune vs non-immune)"
    "Individual database coverage analysis with descriptive plots"
    "Detailed database comparison based on manuscript statistics"
    "Improved PEA vs MS/MS variability distribution plot with clean aesthetics"
)

# Check for optional analysis scripts
optional_scripts=(
    "scripts/analysis/hpa_pea_msms_variability_correlation.R"
    "scripts/analysis/peptideatlas_quantification_comparison.R"
    "scripts/analysis/biomarker_plasma_analysis.R"
    "scripts/visualization/ranked_abundance_plots_improved.R"
)

optional_descriptions=(
    "HPA PEA vs MS/MS variability correlation analysis"
    "PeptideAtlas quantification comparison"
    "Biomarker plasma analysis"
    "Improved ranked abundance plots"
)

# Run each core script
total_scripts=${#scripts[@]}
successful=0
failed=0
skipped=0

print_status "Running core analysis scripts..."
echo "============================================================================"

for i in "${!scripts[@]}"; do
    script="${scripts[$i]}"
    description="${descriptions[$i]}"
    script_name=$(basename "$script")
    log_file="outputs/logs/${script_name%.R}_$(date +%Y%m%d_%H%M%S).log"
    
    # Check if script exists
    if [ ! -f "$script" ]; then
        print_warning "[$((i+1))/$total_scripts] Skipping $script_name - file not found"
        ((skipped++))
        continue
    fi
    
    print_status "[$((i+1))/$total_scripts] Running $script_name - $description"
    
    # Run the script and capture output
    if Rscript "$script" > "$log_file" 2>&1; then
        print_success "$script_name completed successfully"
        print_status "Log saved to: $log_file"
        ((successful++))
    else
        print_error "$script_name failed! Check log: $log_file"
        ((failed++))
        
        # Show last few lines of error log
        print_warning "Last 10 lines of error log:"
        tail -n 10 "$log_file" | while read line; do
            echo "  $line"
        done
    fi
    
    echo # Empty line for readability
done

# Run optional scripts if they exist
print_status "Checking for optional analysis scripts..."
for i in "${!optional_scripts[@]}"; do
    script="${optional_scripts[$i]}"
    description="${optional_descriptions[$i]}"
    script_name=$(basename "$script")
    log_file="outputs/logs/${script_name%.R}_$(date +%Y%m%d_%H%M%S).log"
    
    if [ -f "$script" ]; then
        print_status "Running optional script: $script_name - $description"
        
        if Rscript "$script" > "$log_file" 2>&1; then
            print_success "$script_name completed successfully"
            ((successful++))
        else
            print_warning "$script_name failed (optional script)"
            ((failed++))
        fi
    fi
done

# Summary
echo "============================================================================"
print_status "Analysis Pipeline Complete"
executed=$((successful + failed))
total_available=$((total_scripts + ${#optional_scripts[@]}))

print_success "Successful: $successful scripts completed"
if [ $skipped -gt 0 ]; then
    print_warning "Skipped: $skipped/$total_scripts core scripts (files not found)"
fi
if [ $failed -gt 0 ]; then
    print_error "Failed: $failed scripts"
    print_warning "Check log files in outputs/logs/ directory for details"
elif [ $successful -eq $executed ] && [ $executed -gt 0 ]; then
    print_success "All available scripts completed successfully!"
fi

# Check for generated plots
plot_files=$(find outputs/plots/ -name "*.png" -o -name "*.pdf" -o -name "*.svg" -o -name "*.tiff" 2>/dev/null | wc -l)
if [ $plot_files -gt 0 ]; then
    print_status "Generated $plot_files plot file(s) in outputs/plots/"
    
    # List organized subdirectories
    if [ -d "outputs/plots" ]; then
        echo "  Plot organization:"
        for dir in outputs/plots/*/; do
            if [ -d "$dir" ]; then
                dir_name=$(basename "$dir")
                file_count=$(find "$dir" -name "*.png" -o -name "*.pdf" -o -name "*.tiff" 2>/dev/null | wc -l)
                if [ $file_count -gt 0 ]; then
                    echo "    • $dir_name: $file_count plots"
                fi
            fi
        done
    fi
fi

# Information about the organized structure
if [ $successful -gt 0 ]; then
    echo
    print_status "Repository Organization:"
    echo "  📁 data/ - All data files organized by type"
    echo "    ├── raw/ - Original CSV and data files"
    echo "    ├── processed/ - Cleaned and processed data"
    echo "    └── metadata/ - Documentation and metadata files"
    echo "  📁 scripts/ - All analysis scripts organized by function"
    echo "    ├── analysis/ - Main Figure*.R and analysis scripts"
    echo "    ├── visualization/ - Specialized plotting scripts"
    echo "    ├── data_processing/ - Data cleaning and ID mapping"
    echo "    └── utilities/ - Helper functions and utilities"
    echo "  📁 outputs/ - All generated outputs"
    echo "    ├── plots/ - Generated visualizations"
    echo "    ├── tables/ - Generated tables and summaries"
    echo "    ├── reports/ - Analysis reports and documentation"
    echo "    └── logs/ - Script execution logs"
    echo "  📁 manuscript/ - Manuscript files and bibliography"
    echo "  📁 docs/ - Documentation and README files"
    echo "  📁 config/ - Configuration and dependency files"
    echo
    print_status "All outputs are now organized in the outputs/ directory!"
fi

echo "============================================================================" 