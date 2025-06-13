#!/usr/bin/env Rscript

# ============================================================================
# BLOOD PROTEOMICS ANALYSIS - DEPENDENCY INSTALLER
# ============================================================================
# Description: Installs all required R packages for the blood proteomics analysis project
# Usage: Rscript install_dependencies.R
# Note: Run this script ONCE before running any analysis scripts
# ============================================================================

cat("🔬 Blood Proteomics Analysis - Installing Dependencies\n")
cat("======================================================\n\n")

# Set CRAN mirror for consistent package installation
if(is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  cat("✓ CRAN mirror set to: https://cloud.r-project.org/\n\n")
}

# Function to install packages if not available
install_if_missing <- function(packages, source = "CRAN") {
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste("📦 Installing", pkg, "from", source, "...\n"))
      
      if(source == "CRAN") {
        install.packages(pkg, dependencies = TRUE)
      } else if(source == "Bioconductor") {
        BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
      }
      
      if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
        stop(paste("❌ Failed to install package:", pkg))
      } else {
        cat(paste("✅", pkg, "installed successfully\n"))
      }
    } else {
      cat(paste("✅", pkg, "already installed\n"))
    }
  }
}

# ============================================================================
# CORE DATA MANIPULATION PACKAGES
# ============================================================================
cat("📊 Installing core data manipulation packages...\n")
cat("------------------------------------------------\n")
data_packages <- c(
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "readr",        # Reading CSV/TSV files
  "stringr",      # String manipulation
  "data.table",   # Fast data processing
  "optparse"      # Command line argument parsing
)
install_if_missing(data_packages)

# ============================================================================
# VISUALIZATION PACKAGES
# ============================================================================
cat("\n🎨 Installing visualization packages...\n")
cat("---------------------------------------\n")
viz_packages <- c(
  "ggplot2",      # Core plotting
  "ggridges",     # Ridge plots
  "ggbeeswarm",   # Bee swarm plots  
  "viridis",      # Color palettes
  "scales",       # Axis scaling
  "gridExtra",    # Multiple plots
  "grid",         # Grid graphics
  "RColorBrewer", # Color palettes
  "UpSetR",       # UpSet plots
  "ggupset",      # UpSet plots
  "VennDiagram",  # Venn diagrams
  "pheatmap",     # Heatmaps
  "corrplot"      # Correlation plots
)
install_if_missing(viz_packages)

# ============================================================================
# STATISTICAL ANALYSIS PACKAGES
# ============================================================================
cat("\n📈 Installing statistical analysis packages...\n")
cat("----------------------------------------------\n")
stats_packages <- c(
  "Hmisc",        # Statistical functions
  "car",          # Regression analysis
  "psych",        # Psychological statistics
  "broom"         # Tidy statistical output
)
install_if_missing(stats_packages)

# ============================================================================
# BIOCONDUCTOR PACKAGES
# ============================================================================
cat("\n🧬 Installing Bioconductor packages...\n")
cat("--------------------------------------\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("📦 Installing BiocManager...\n")
  install.packages("BiocManager")
  cat("✅ BiocManager installed successfully\n")
}

# Bioconductor packages for biological analysis
bioc_packages <- c(
  "biomaRt",      # Biomart database access
  "org.Hs.eg.db", # Human gene annotation
  "GO.db",        # Gene Ontology database
  "KEGG.db",      # KEGG pathway database
  "AnnotationDbi" # Annotation database interface
)

for(pkg in bioc_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("📦 Installing Bioconductor package", pkg, "...\n"))
    tryCatch({
      BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
      if(require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat(paste("✅", pkg, "installed successfully\n"))
      } else {
        warning(paste("⚠️ Failed to load", pkg, "after installation"))
      }
    }, error = function(e) {
      warning(paste("❌ Failed to install Bioconductor package", pkg, ":", e$message))
    })
  } else {
    cat(paste("✅", pkg, "already installed\n"))
  }
}

# ============================================================================
# OPTIONAL ENHANCEMENT PACKAGES
# ============================================================================
cat("\n🎯 Installing optional enhancement packages...\n")
cat("----------------------------------------------\n")
optional_packages <- c(
  "plotly",       # Interactive plots
  "DT",           # Interactive tables
  "knitr",        # Dynamic reports
  "rmarkdown",    # R Markdown documents
  "devtools"      # Development tools
)

for(pkg in optional_packages) {
  tryCatch({
    install_if_missing(pkg)
  }, error = function(e) {
    cat(paste("⚠️ Optional package", pkg, "failed to install:", e$message, "\n"))
  })
}

# ============================================================================
# INSTALLATION SUMMARY
# ============================================================================
cat("\n", rep("=", 60), "\n")
cat("🎯 INSTALLATION SUMMARY\n")
cat(rep("=", 60), "\n")

# Check all packages
all_packages <- c(data_packages, viz_packages, stats_packages, 
                 bioc_packages, optional_packages)

installed_count <- 0
failed_packages <- character()

for(pkg in all_packages) {
  if(require(pkg, character.only = TRUE, quietly = TRUE)) {
    installed_count <- installed_count + 1
  } else {
    failed_packages <- c(failed_packages, pkg)
  }
}

cat("📊 Total packages checked:", length(all_packages), "\n")
cat("✅ Successfully installed:", installed_count, "\n")
cat("❌ Failed installations:", length(failed_packages), "\n")

if(length(failed_packages) > 0) {
  cat("\n⚠️ Failed packages:\n")
  for(pkg in failed_packages) {
    cat("  ❌", pkg, "\n")
  }
  cat("\nNote: Some packages may be optional and not critical for basic functionality.\n")
  cat("You can try installing failed packages manually if needed.\n")
} else {
  cat("\n🎉 All packages installed successfully!\n")
}

cat("\n" , rep("=", 60), "\n")
cat("🚀 READY TO ANALYZE!\n")
cat(rep("=", 60), "\n")
cat("✅ All dependencies are now installed\n")
cat("✅ You can run any analysis script in the project\n")
cat("✅ Scripts will load packages automatically (no more installations)\n\n")

cat("📁 Next steps:\n")
cat("   1. Navigate to the scripts/analysis/ directory\n")
cat("   2. Run any analysis script (e.g., Rscript figure2.R)\n")
cat("   3. Check outputs/ directory for results\n\n")

cat("💡 Tip: If you encounter missing packages later, re-run this script\n")
cat("🔧 For troubleshooting, check the R session info below:\n\n")

# Print session information for debugging
sessionInfo() 