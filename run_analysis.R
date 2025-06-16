#!/usr/bin/env Rscript

# ===============================================================================
# Blood Review Data Analysis Runner (R Version)
#
# Purpose: Execute the full blood proteomics analysis pipeline.
# This script replaces run_analysis.sh for better platform independence.
# ===============================================================================

# --- 1. Preamble: Load Packages and Define Helpers ---

# Set CRAN mirror for stable package installation in non-interactive sessions
if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}

# Exit on any error
options(error = function() {
  message(crayon::red("\n[ERROR] A critical error occurred. Aborting script."))
  q("no", status = 1, runLast = FALSE)
})

# Load required packages for this runner script
if (!require(argparse)) install.packages("argparse", quiet = TRUE)
if (!require(crayon)) install.packages("crayon", quiet = TRUE)

library(argparse)
library(crayon)

# Helper functions for colored console output
print_status <- function(...) cat(blue("[INFO]"), ..., "\n")
print_success <- function(...) cat(green("[SUCCESS]"), ..., "\n")
print_warning <- function(...) cat(yellow("[WARNING]"), ..., "\n")
print_error <- function(...) cat(red("[ERROR]"), ..., "\n")

# --- 2. Argument Parsing ---

parser <- ArgumentParser(description = "Run the full blood proteomics analysis pipeline.")
parser$add_argument("--force-mapping", action = "store_true", default = FALSE,
                    help = "Force re-mapping of protein IDs to gene symbols.")
parser$add_argument("--formats", type = "character", default = "svg",
                    help = "Comma-separated list of plot formats (e.g., 'svg,png').")
args <- parser$parse_args()

# --- 3. Script Setup and Prerequisite Checks ---

# Print header
cat(bold$cyan(rep("=", 70)), "\n")
cat(bold$cyan("               BLOOD REVIEW DATA ANALYSIS (R RUNNER)"), "\n")
cat(bold$cyan("             Plasma and Serum Protein Quantification"), "\n")
cat(bold$cyan(rep("=", 70)), "\n\n")

# Verify that all required data files exist
print_status("Verifying required data files...")
required_files <- c(
    "data/raw/peptideatlas/peptideatlas.csv", "data/raw/hpa/hpa_ms.csv",
    "data/raw/hpa/hpa_pea.csv", "data/raw/hpa/hpa_immunoassay_plasma.csv",
    "data/raw/gpmdb/gpmdb_plasma.csv", "data/raw/paxdb/paxdb_plasma.csv",
    "data/raw/gpmdb/gpmdb_serum.csv", "data/raw/paxdb/paxdb_serum.csv",
    "data/raw/hpa/hpa_immunoassay_serum.csv", "data/raw/paxdb/paxdb_b_cell.csv",
    "data/raw/paxdb/paxdb_cd4.csv", "data/raw/paxdb/paxdb_cd8.csv",
    "data/raw/paxdb/paxdb_monocyte.csv", "data/raw/paxdb/paxdb_nk.csv",
    "data/raw/proteomexchange/pxd004352.csv", "data/raw/proteomexchange/pxd025174.csv",
    "data/raw/proteomexchange/pxd040957_cd8.csv", "data/raw/proteomexchange/pxd040957_macrophages.csv"
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    print_error("Missing required data files:")
    for (file in missing_files) {
        cat("  - ", file, "\n")
    }
    q(status = 1)
}
print_success("All required data files found.")

# Create output directories if they don't exist
print_status("Setting up output directories...")
dir.create("outputs", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# --- 4. Analysis Execution ---

# Define the sequence of analysis scripts to run
analysis_scripts <- c(
    "01_plasma_protein_analysis.R",
    "02_peptideatlas_quantification_analysis.R",
    "03_biomarker_plasma_analysis.R",
    "04_serum_protein_analysis.R",
    "05_celltype_analysis.R"
)

# Function to run a single analysis script
run_script <- function(script_name, cli_args) {
    script_path <- file.path("scripts", script_name)
    print_status(bold(paste("Running:", script_name)))
    
    # Construct arguments for the script
    rscript_args <- c(script_path)
    if (cli_args$force_mapping) {
        rscript_args <- c(rscript_args, "--force-mapping")
    }
    rscript_args <- c(rscript_args, "--formats", cli_args$formats)
    
    # Execute the script in a separate R process
    result <- system2("Rscript", args = rscript_args, stdout = TRUE, stderr = TRUE)
    
    # Check for errors by looking for non-zero exit status (handled by options(error))
    # or specific error messages in the output.
    if (any(grepl("error", result, ignore.case = TRUE))) {
        print_error(paste(script_name, "failed!"))
        cat("--- Script Output ---\n")
        cat(result, sep = "\n")
        cat("-----------------------\n")
        stop(paste("Execution of", script_name, "halted due to an error."))
    }
    
    print_success(paste(script_name, "completed successfully!"))
    cat("\n")
}

# Loop through and execute each script
for (script in analysis_scripts) {
    run_script(script, args)
}

# --- 5. Final Summary ---

cat(bold$cyan(rep("=", 70)), "\n")
cat(bold$cyan("                           ANALYSIS COMPLETE"), "\n")
cat(bold$cyan(rep("=", 70)), "\n\n")

print_success("All analysis scripts executed successfully.")
cat("  📊 Check the 'outputs/' and 'plots/' directories for results.\n")
cat("  📖 Review the 'PROJECT_RESULTS_SUMMARY.md' for a full overview.\n\n")

cat(bold("🚀 ANALYSIS FEATURES:\n"))
cat("• Modular code architecture with centralized configuration\n")
cat("• Enhanced normalization comparisons for better cross-database analysis\n")
cat("• Standardized plotting themes and consistent styling\n")
cat("• Automatic report generation with comprehensive statistics\n")
cat("• Improved error handling and data validation\n") 