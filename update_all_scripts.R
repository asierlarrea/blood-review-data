#!/usr/bin/env Rscript

# Helper script to update all R analysis scripts to use the new dependency system
# This script identifies and updates any remaining scripts with old dependency code

message("Updating all R scripts to use new dependency management system...")

# Find all R scripts in the current directory
r_scripts <- list.files(pattern = "\\.R$", full.names = TRUE)
r_scripts <- r_scripts[!grepl("install_dependencies|load_packages|update_all_scripts", r_scripts)]

message(paste("Found", length(r_scripts), "R scripts to check"))

for(script in r_scripts) {
  message(paste("Checking:", basename(script)))
  
  # Read script content
  lines <- readLines(script, warn = FALSE)
  
  # Check if it has old dependency management
  has_install_func <- any(grepl("install_if_missing.*function", lines))
  has_install_call <- any(grepl("install_if_missing\\(", lines))
  has_cran_mirror <- any(grepl("options\\(repos.*CRAN", lines))
  
  if(has_install_func || has_install_call || has_cran_mirror) {
    message(paste("  → Needs updating:", basename(script)))
    
    # Find where dependency management starts and ends
    start_line <- which(grepl("# Set CRAN mirror|# Function to install", lines))[1]
    end_line <- which(grepl("install_if_missing\\(|install\\.packages\\(", lines))
    end_line <- max(end_line) + 1
    
    if(!is.na(start_line) && length(end_line) > 0) {
      # Extract package list
      pkg_line <- lines[grepl("required_packages.*<-.*c\\(", lines)]
      if(length(pkg_line) > 0) {
        # Extract packages from the line
        pkg_match <- regexpr("c\\([^)]+\\)", pkg_line)
        if(pkg_match > 0) {
          packages_str <- substr(pkg_line, pkg_match, pkg_match + attr(pkg_match, "match.length") - 1)
          
          # Create new loading code
          new_lines <- c(
            "# Load required packages (install via install_dependencies.R if needed)",
            "source(\"load_packages.R\")",
            paste("required_packages <-", packages_str),
            "load_packages(required_packages)"
          )
          
          # Replace old code with new code
          updated_lines <- c(
            lines[1:(start_line-1)],
            new_lines,
            lines[(end_line+1):length(lines)]
          )
          
          # Write updated script
          writeLines(updated_lines, script)
          message(paste("  ✓ Updated:", basename(script)))
        }
      }
    }
  } else {
    message(paste("  ✓ Already clean:", basename(script)))
  }
}

message("\nScript update completed!")
message("All scripts now use the centralized dependency management system.") 