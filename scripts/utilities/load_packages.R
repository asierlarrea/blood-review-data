# ============================================================================
# PACKAGE LOADING UTILITIES
# ============================================================================
# Description: Centralized package loading for blood proteomics analysis
# Note: All packages must be installed first using install_dependencies.R
# ============================================================================

# Function to load packages (assumes they are already installed)
load_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("❌ Package", pkg, "not found.", 
                "\n💡 Please run 'Rscript install_dependencies.R' first to install all dependencies."))
    }
  }
}

# Function to get the project root directory
get_project_root <- function() {
  # Start from current working directory and look for key directories
  current_dir <- getwd()
  
  # Check if we're already in the root (has data/, scripts/, outputs/ directories)
  if (all(dir.exists(c("data", "scripts", "outputs")))) {
    return(current_dir)
  }
  
  # If we're in a subdirectory, go up to find the root
  parent_dir <- dirname(current_dir)
  while (parent_dir != current_dir) {
    if (all(dir.exists(file.path(parent_dir, c("data", "scripts", "outputs"))))) {
      return(parent_dir)
    }
    current_dir <- parent_dir
    parent_dir <- dirname(current_dir)
  }
  
  # If not found, return current directory with warning
  warning("Project root not found. Using current directory.")
  return(getwd())
}

# Function to create output directories if they don't exist
ensure_output_dirs <- function() {
  root_dir <- get_project_root()
  dirs_to_create <- c(
    "outputs/plots",
    "outputs/tables", 
    "outputs/reports",
    "outputs/logs"
  )
  
  for (dir in dirs_to_create) {
    full_path <- file.path(root_dir, dir)
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = TRUE)
      message(paste("Created directory:", full_path))
    }
  }
}

# Function to get data file path
get_data_path <- function(filename, subdir = "raw") {
  root_dir <- get_project_root()
  return(file.path(root_dir, "data", subdir, filename))
}

# Function to get output path
get_output_path <- function(filename, subdir = "plots") {
  root_dir <- get_project_root()
  return(file.path(root_dir, "outputs", subdir, filename))
}
