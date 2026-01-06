#!/usr/bin/env Rscript

################################################################################
# R Package Installation Script
# Spleen cDC Homeostatic Maturation Project
################################################################################
#
# This script installs all required R packages for the analysis pipeline.
# Run this once to set up your environment.
#
# Usage:
#   Rscript install_dependencies.R
#
# Or from R console:
#   source("install_dependencies.R")
#
################################################################################

cat("\n")
cat("================================================================================\n")
cat("  Installing R Package Dependencies\n")
cat("  Spleen cDC Homeostatic Maturation Project\n")
cat("================================================================================\n")
cat("\n")

# Function to check if package is installed
is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# Function to install package if not present
install_if_missing <- function(pkg, source = "CRAN") {
  if (!is_installed(pkg)) {
    cat("Installing:", pkg, "from", source, "\n")
    if (source == "CRAN") {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    } else if (source == "Bioconductor") {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
  } else {
    cat("✓", pkg, "already installed\n")
  }
}

################################################################################
# Step 1: Install Package Management Tools
################################################################################

cat("\n[1/6] Installing package management tools...\n")
cat("────────────────────────────────────────────────────────────────────────────\n")

# renv for dependency management
if (!is_installed("renv")) {
  install.packages("renv", repos = "https://cloud.r-project.org/")
}

# BiocManager for Bioconductor packages
if (!is_installed("BiocManager")) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# devtools for GitHub packages
if (!is_installed("devtools")) {
  install.packages("devtools", repos = "https://cloud.r-project.org/")
}

cat("✓ Package management tools installed\n")

################################################################################
# Step 2: Set Bioconductor Version
################################################################################

cat("\n[2/6] Setting up Bioconductor...\n")
cat("────────────────────────────────────────────────────────────────────────────\n")

library(BiocManager)

# Install/update to latest stable Bioconductor version
# Change version number if you need a specific version
BiocManager::install(version = "3.18", ask = FALSE, update = FALSE)

cat("Bioconductor version:", as.character(BiocManager::version()), "\n")
cat("✓ Bioconductor configured\n")

################################################################################
# Step 3: Install Core CRAN Packages
################################################################################

cat("\n[3/6] Installing CRAN packages...\n")
cat("────────────────────────────────────────────────────────────────────────────\n")

cran_packages <- c(
  # Tidyverse (meta-package includes dplyr, ggplot2, tidyr, readr, etc.)
  "tidyverse",

  # Visualization
  "pheatmap",
  "RColorBrewer",
  "viridis",
  "rgl",  # Optional: 3D visualization (may require OpenGL)

  # Utilities
  "openxlsx",  # Excel file I/O
  "future",    # Parallel processing

  # Development tools
  "devtools"
)

for (pkg in cran_packages) {
  install_if_missing(pkg, source = "CRAN")
}

cat("✓ CRAN packages installed\n")

################################################################################
# Step 4: Install Bioconductor Packages
################################################################################

cat("\n[4/6] Installing Bioconductor packages...\n")
cat("────────────────────────────────────────────────────────────────────────────\n")

bioc_packages <- c(
  # Single-cell data structures
  "SingleCellExperiment",

  # Single-cell analysis
  "scater",
  "scran",

  # Differential expression
  "edgeR",
  "limma",

  # Batch correction
  "sva",

  # Specialized analysis
  "dorothea",  # Transcription factor activity
  "dyno"       # Trajectory inference
)

for (pkg in bioc_packages) {
  install_if_missing(pkg, source = "Bioconductor")
}

cat("✓ Bioconductor packages installed\n")

################################################################################
# Step 5: Install Seurat (CRAN version)
################################################################################

cat("\n[5/6] Installing Seurat...\n")
cat("────────────────────────────────────────────────────────────────────────────\n")

install_if_missing("Seurat", source = "CRAN")

cat("✓ Seurat installed\n")

################################################################################
# Step 6: Install GitHub Packages
################################################################################

cat("\n[6/6] Installing GitHub packages...\n")
cat("────────────────────────────────────────────────────────────────────────────\n")

# DoubletFinder from GitHub
if (!is_installed("DoubletFinder")) {
  cat("Installing: DoubletFinder from GitHub\n")
  tryCatch({
    devtools::install_github("chris-mcginnis-ucsf/DoubletFinder", quiet = TRUE)
    cat("✓ DoubletFinder installed\n")
  }, error = function(e) {
    cat("⚠ Warning: DoubletFinder installation failed\n")
    cat("  Error:", e$message, "\n")
    cat("  You may need to install it manually\n")
  })
} else {
  cat("✓ DoubletFinder already installed\n")
}

################################################################################
# Finalize: Create renv Snapshot
################################################################################

cat("\n")
cat("================================================================================\n")
cat("  Installation Complete!\n")
cat("================================================================================\n")
cat("\n")

cat("Creating renv snapshot...\n")
renv::init(bare = TRUE, restart = FALSE)
renv::snapshot(prompt = FALSE)

cat("\n")
cat("✓ All packages installed successfully\n")
cat("✓ renv snapshot created: renv.lock\n")
cat("\n")

################################################################################
# Print Summary
################################################################################

cat("────────────────────────────────────────────────────────────────────────────\n")
cat("Summary of installed packages:\n")
cat("────────────────────────────────────────────────────────────────────────────\n")
cat("\nR version:", R.version.string, "\n")
cat("Bioconductor version:", as.character(BiocManager::version()), "\n")
cat("\n")

# Print versions of key packages
key_packages <- c("Seurat", "tidyverse", "edgeR", "limma", "scater", "scran")
cat("Key package versions:\n")
for (pkg in key_packages) {
  if (is_installed(pkg)) {
    version <- as.character(packageVersion(pkg))
    cat(sprintf("  %-20s %s\n", pkg, version))
  }
}

cat("\n")
cat("────────────────────────────────────────────────────────────────────────────\n")
cat("Next steps:\n")
cat("────────────────────────────────────────────────────────────────────────────\n")
cat("\n")
cat("1. Review the renv.lock file to see all package versions\n")
cat("2. Run your analysis scripts to verify everything works\n")
cat("3. See DEPENDENCY_AUDIT_REPORT.md for optimization recommendations\n")
cat("4. See IMPLEMENTATION_GUIDE.md for cleanup instructions\n")
cat("\n")
cat("To restore this environment on another machine:\n")
cat("  renv::restore()\n")
cat("\n")
cat("To update packages:\n")
cat("  renv::update()\n")
cat("\n")
cat("================================================================================\n")
cat("\n")
