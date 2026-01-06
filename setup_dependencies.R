#!/usr/bin/env Rscript
# Setup script for initializing renv and installing dependencies
# for Spleen_cDC_Homeostatic_Maturation project

cat("========================================\n")
cat("Dependency Setup for cDC Analysis Project\n")
cat("========================================\n\n")

# Install renv if not present
if (!requireNamespace("renv", quietly = TRUE)) {
  cat("Installing renv package...\n")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

cat("\nInitializing renv for this project...\n")
renv::init()

cat("\nInstalling Bioconductor dependencies...\n")
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Set Bioconductor version (adjust as needed)
# BiocManager::install(version = "3.18")

# Install Bioconductor packages
bioc_packages <- c(
  "limma",
  "edgeR",
  "scater",
  "scran",
  "SingleCellExperiment",
  "dorothea"
)

cat("Installing Bioconductor packages:", paste(bioc_packages, collapse = ", "), "\n")
BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)

# Install CRAN packages
cran_packages <- c(
  "Seurat",
  "tidyverse",
  "gridExtra",
  "openxlsx",
  "RColorBrewer",
  "viridis",
  "pheatmap",
  "future",
  "dyno",
  "DoubletFinder"
)

cat("\nInstalling CRAN packages:", paste(cran_packages, collapse = ", "), "\n")
install.packages(cran_packages, repos = "https://cloud.r-project.org")

# Optional packages (commented out by default - uncomment if needed)
# install.packages(c("sqldf", "fields", "modes", "rgl", "gplots"))

cat("\nTaking snapshot of installed packages...\n")
renv::snapshot()

cat("\n========================================\n")
cat("Setup complete!\n")
cat("========================================\n")
cat("\nDependency lockfile created: renv.lock\n")
cat("To restore this environment on another machine, run:\n")
cat("  renv::restore()\n\n")

# Verify installation
cat("Verifying key packages...\n")
required_pkgs <- c("Seurat", "tidyverse", "limma", "edgeR", "scater")
for (pkg in required_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("  ✓", pkg, "\n")
  } else {
    cat("  ✗", pkg, "- FAILED\n")
  }
}

cat("\nRun renv::status() to check for any issues.\n")
