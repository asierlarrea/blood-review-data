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

# Check if we're in the right directory
if [ ! -f "Figure2.R" ]; then
    print_error "Figure2.R not found. Please run this script from the project directory."
    exit 1
fi

# Check for new analysis scripts
new_scripts=("Figure6_database_correlation.R" "Figure7_abundance_analysis.R" "Figure8_functional_analysis.R" "Figure9_comparative_analysis.R")
missing_new=0
for script in "${new_scripts[@]}"; do
    if [ ! -f "$script" ]; then
        ((missing_new++))
    fi
done

if [ $missing_new -gt 0 ]; then
    print_warning "$missing_new new analysis script(s) not found. Running original analysis only."
    print_status "To get the complete extended analysis, ensure all Figure6-9 scripts are present."
else
    print_status "Extended analysis detected - running complete pipeline with 8 analysis scripts"
fi

print_status "Starting Blood Proteomics Analysis Pipeline"
print_status "Working directory: $(pwd)"

# Create output directory for logs
mkdir -p logs

# Array of scripts to run
scripts=(
    "Figure2.R" 
    "Figure3.R" 
    "Figure4.R" 
    "Figure5.R"
    "Figure6_database_correlation.R"
    "Figure7_abundance_analysis.R"
    "Figure8_functional_analysis.R"
    "Figure9_comparative_analysis.R"
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
)

# Run each script
total_scripts=${#scripts[@]}
successful=0
failed=0
skipped=0

for i in "${!scripts[@]}"; do
    script="${scripts[$i]}"
    description="${descriptions[$i]}"
    log_file="logs/${script%.R}_$(date +%Y%m%d_%H%M%S).log"
    
    # Check if script exists
    if [ ! -f "$script" ]; then
        print_warning "[$((i+1))/$total_scripts] Skipping $script - file not found"
        ((skipped++))
        continue
    fi
    
    print_status "[$((i+1))/$total_scripts] Running $script - $description"
    
    # Run the script and capture output
    if Rscript "$script" > "$log_file" 2>&1; then
        print_success "$script completed successfully"
        print_status "Log saved to: $log_file"
        ((successful++))
    else
        print_error "$script failed! Check log: $log_file"
        ((failed++))
        
        # Show last few lines of error log
        print_warning "Last 10 lines of error log:"
        tail -n 10 "$log_file" | while read line; do
            echo "  $line"
        done
    fi
    
    echo # Empty line for readability
done

# Summary
echo "============================================================================"
print_status "Analysis Pipeline Complete"
executed=$((successful + failed))
print_success "Successful: $successful/$executed scripts executed"
if [ $skipped -gt 0 ]; then
    print_warning "Skipped: $skipped/$total_scripts scripts (files not found)"
fi
if [ $failed -gt 0 ]; then
    print_error "Failed: $failed/$executed scripts"
    print_warning "Check log files in the 'logs/' directory for details"
elif [ $successful -eq $executed ] && [ $executed -gt 0 ]; then
    print_success "All available scripts completed successfully!"
fi

# Check for generated plots (if any)
plot_files=$(find . -name "*.png" -o -name "*.pdf" -o -name "*.svg" -o -name "*.tiff" -o -name "plots/*.tiff" 2>/dev/null | wc -l)
if [ $plot_files -gt 0 ]; then
    print_status "Generated $plot_files plot file(s)"
    
    # Check specifically for plots in organized directories
    plots_dir_files=$(find plots/ -name "*.tiff" 2>/dev/null | wc -l)
    if [ $plots_dir_files -gt 0 ]; then
        print_status "New analysis plots saved in organized directories: $plots_dir_files TIFF files"
        
        # List organized subdirectories if they exist
        if [ -d "plots" ]; then
            echo "  Plot organization:"
            for dir in plots/*/; do
                if [ -d "$dir" ]; then
                    dir_name=$(basename "$dir")
                    file_count=$(find "$dir" -name "*.tiff" 2>/dev/null | wc -l)
                    if [ $file_count -gt 0 ]; then
                        echo "    • $dir_name: $file_count plots"
                    fi
                fi
            done
        fi
    fi
fi

# Additional information for extended analysis
if [ $successful -gt 4 ]; then
    echo
    print_status "Extended analysis completed! Check the following:"
    echo "  • Original plots: Rplots.pdf (if generated)"
    echo "  • New analysis plots: plots/ directory with organized subdirectories"
    echo "  • Analysis recommendations: ANALYSIS_RECOMMENDATIONS.md"
    echo "  • Script logs: logs/ directory"
    echo "  • Total new visualizations: ~18 publication-ready plots"
    echo
    print_status "Plot organization structure:"
    echo "  📁 01_Database_Analysis - Database coverage and correlations"
    echo "  📁 02_Coverage_Analysis - Cell type data completeness"  
    echo "  📁 03_Abundance_Analysis - Protein concentration distributions"
    echo "  📁 04_Cell_Type_Comparison - Cell type comparative analysis"
    echo "  📁 05_Functional_Analysis - Protein functional categories"
    echo "  📁 06_Protein_Classes - Immunoglobulins and protein classes"
    echo "  📁 07_Diversity_Analysis - Cell type diversity metrics"
    echo "  📁 08_Database_Comparison - Database uniqueness analysis"
    echo "  📁 09_Efficiency_Analysis - Coverage efficiency metrics"
fi

echo "============================================================================" 