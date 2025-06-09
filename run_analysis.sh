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

print_status "Starting Blood Proteomics Analysis Pipeline"
print_status "Working directory: $(pwd)"

# Create output directory for logs
mkdir -p logs

# Array of scripts to run
scripts=("Figure2.R" "Figure3.R" "Figure4.R" "Figure5.R")
descriptions=(
    "UpSet plot of plasma proteins"
    "Cell type bubble plot"
    "Cell type UpSet plot" 
    "Intensity distribution boxplots"
)

# Run each script
total_scripts=${#scripts[@]}
successful=0
failed=0

for i in "${!scripts[@]}"; do
    script="${scripts[$i]}"
    description="${descriptions[$i]}"
    log_file="logs/${script%.R}_$(date +%Y%m%d_%H%M%S).log"
    
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
print_success "Successful: $successful/$total_scripts scripts"
if [ $failed -gt 0 ]; then
    print_error "Failed: $failed/$total_scripts scripts"
    print_warning "Check log files in the 'logs/' directory for details"
else
    print_success "All scripts completed successfully!"
fi

# Check for generated plots (if any)
plot_files=$(find . -name "*.png" -o -name "*.pdf" -o -name "*.svg" 2>/dev/null | wc -l)
if [ $plot_files -gt 0 ]; then
    print_status "Generated $plot_files plot file(s)"
fi

echo "============================================================================" 