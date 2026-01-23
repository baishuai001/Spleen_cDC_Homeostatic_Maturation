################################################################################
## Script for analyzing ADT data using CITEViz
## CITEViz provides a FlowJo-like gating workflow for CITE-seq data
##
## Reference: CITEViz: interactively classify cell populations in CITE-Seq
##            via a flow cytometry-like gating workflow using R-Shiny
##            BMC Bioinformatics (2024)
##            https://github.com/maxsonBraunLab/CITEViz
##
## Input: Seurat object with processed ADT assay (seurat_obj_annotated.rds)
## Output: Interactive CITEViz Shiny app for ADT gating analysis
################################################################################

################################################################################
########################## INSTALLATION ########################################
################################################################################

# Install CITEViz from GitHub (run once)
# Note: The GitHub repo is named "CITE-Viz" (with hyphen) but package is "CITEViz"
if (!requireNamespace("CITEViz", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("maxsonBraunLab/CITE-Viz")
}

################################################################################
########################## LOAD PACKAGES #######################################
################################################################################

library(Seurat)
library(CITEViz)
library(dplyr)

################################################################################
########################## SETUP PATHS #########################################
################################################################################

# Specify your working directory and file paths
# Adjust these paths based on where your files are located
setwd("/home/user/Spleen_cDC_Homeostatic_Maturation/")

# Input: Seurat object from previous analysis
# This should be the output from 3.RNA-ADT_HPCscript.R and 4.script_CITEseq_SAM_WT_aggr.R
input_seurat <- "path/to/seurat_obj_annotated.rds"

# Example paths based on your project structure:
# Option 1: If using the full WT_aggr object
# input_seurat <- "SAM2and3_WT/results/Robjects/seuratObj_SAM2and3_WT_ADT.rds"

# Option 2: If using cDC1 subset
# input_seurat <- "SAM2and3_WT_subset1/results/Robjects/seuratObj_SAM2and3_WT_subset1_ADT.rds"

# Option 3: If using cDC2 subset
# input_seurat <- "SAM2and3_WT_subset2/results/Robjects/seuratObj_SAM2and3_WT_subset2_ADT.rds"

################################################################################
########################## LOAD SEURAT OBJECT ##################################
################################################################################

cat("Loading Seurat object...\n")
seuratObj <- readRDS(input_seurat)

# Check that the object has ADT data
if (!"ADT" %in% names(seuratObj@assays)) {
  stop("Error: Seurat object does not contain ADT assay.
       Please use an object generated from 3.RNA-ADT_HPCscript.R")
}

cat("Seurat object loaded successfully!\n")
cat("Number of cells:", ncol(seuratObj), "\n")
cat("Number of ADT markers:", nrow(seuratObj@assays$ADT), "\n")

################################################################################
########################## CHECK ADT NORMALIZATION #############################
################################################################################

# CITEViz requires normalized ADT data
# Check if ADT assay has normalized data
if (!"data" %in% slotNames(seuratObj@assays$ADT)) {
  cat("Warning: ADT assay may not be normalized.\n")
  cat("Performing CLR normalization...\n")
  seuratObj <- NormalizeData(seuratObj,
                             assay = "ADT",
                             normalization.method = "CLR",
                             verbose = TRUE)
}

################################################################################
########################## PREPARE OBJECT FOR CITEViz ##########################
################################################################################

# CITEViz works best with:
# 1. Normalized ADT data (CLR transformation) - already done in script 3
# 2. Dimensional reductions (UMAP/tSNE) - already done in script 3
# 3. Cell metadata including clusters - already done in script 4

# Verify dimensional reductions exist
cat("\nAvailable reductions:\n")
print(names(seuratObj@reductions))

# Verify ADT markers available
cat("\nAvailable ADT markers:\n")
print(rownames(seuratObj@assays$ADT))

# Display available metadata columns for grouping/coloring
cat("\nAvailable metadata columns:\n")
print(colnames(seuratObj@meta.data))

################################################################################
########################## LAUNCH CITEViz APP ##################################
################################################################################

cat("\n==============================================================\n")
cat("Launching CITEViz interactive app...\n")
cat("==============================================================\n")
cat("\nInstructions for using CITEViz:\n")
cat("1. The app will open in your web browser\n")
cat("2. Your Seurat object is already loaded\n")
cat("3. Use the 'Gating' tab to perform FlowJo-like gating:\n")
cat("   - Select ADT markers for X and Y axes\n")
cat("   - Draw gates to subset cell populations\n")
cat("   - Gates are immediately reflected in UMAP/tSNE plots\n")
cat("4. Use the 'QC' tab to view quality control metrics\n")
cat("5. Use the 'Feature Plot' tab to visualize marker expression\n")
cat("6. Save gated populations back to your Seurat object\n")
cat("\nKey features similar to FlowJo:\n")
cat("- Interactive 2D scatter plots of ADT markers\n")
cat("- Sequential gating (gate within gates)\n")
cat("- Real-time visualization of gated cells\n")
cat("- Export gated cell identities\n")
cat("==============================================================\n\n")

# Launch the CITEViz Shiny app with your Seurat object
# The app will open in your default web browser
run_app(seurat_object = seuratObj)

################################################################################
########################## SAVE GATED POPULATIONS ##############################
################################################################################

# After gating in CITEViz, you can save the results
# CITEViz allows you to export gated cell populations
# These can be saved as new metadata columns in your Seurat object

# Example workflow after gating (run after closing CITEViz app):
# 1. CITEViz creates a new metadata column with your gate names
# 2. Save the updated object:
# saveRDS(seuratObj, file = "path/to/seuratObj_with_CITEViz_gates.rds")

# Example: Extract specific gated populations for downstream analysis
# gated_cells <- WhichCells(seuratObj, expression = gate_name == "Your_Gate")
# subset_obj <- subset(seuratObj, cells = gated_cells)

################################################################################
########################## ALTERNATIVE: MANUAL GATING ##########################
################################################################################

# If you want to manually gate specific populations without the interactive app,
# you can use Seurat functions directly:

# Example 1: Gate CD11c+ CD8a+ cells (Resident cDC1s)
# DefaultAssay(seuratObj) <- "ADT"
# cd11c_high <- WhichCells(seuratObj, expression = `CD11c` > 2)
# cd8a_high <- WhichCells(seuratObj, expression = `CD8a` > 1.5)
# cDC1_gated <- intersect(cd11c_high, cd8a_high)
# seuratObj$manual_gate_cDC1 <- ifelse(colnames(seuratObj) %in% cDC1_gated,
#                                      "cDC1", "Other")

# Example 2: Gate CD11b+ CD172a+ cells (Resident cDC2s)
# cd11b_high <- WhichCells(seuratObj, expression = `CD11b` > 2)
# cd172a_high <- WhichCells(seuratObj, expression = `CD172a` > 1.5)
# cDC2_gated <- intersect(cd11b_high, cd172a_high)
# seuratObj$manual_gate_cDC2 <- ifelse(colnames(seuratObj) %in% cDC2_gated,
#                                      "cDC2", "Other")

# Visualize manual gates
# DimPlot(seuratObj, reduction = "ADT_umap", group.by = "manual_gate_cDC1")
# FeaturePlot(seuratObj, reduction = "ADT_umap",
#             features = c("adt_CD11c", "adt_CD8a"),
#             min.cutoff = "q05", max.cutoff = "q95")

################################################################################
########################## ADDITIONAL NOTES ####################################
################################################################################

# CITEViz Advantages:
# - Interactive gating similar to FlowJo
# - Immediate visualization in UMAP/tSNE space
# - Easy to refine cell population definitions
# - No need to manually determine threshold values
# - Export gates for reproducible analysis

# CITEViz Requirements:
# - Pre-processed Seurat object (completed in scripts 3 & 4)
# - Normalized ADT data (CLR transformation recommended)
# - Basic R knowledge
# - Familiarity with flow cytometry concepts

# Useful ADT markers for DC populations (based on your study):
# - CD11c: All DCs
# - CD8a: cDC1s
# - CD11b: cDC2s
# - CD172a (SIRPa): cDC2s
# - MHC-II (I-A/I-E): Mature DCs
# - CD86: Mature DCs
# - CCR7: Migratory DCs
# - ESAM: Resident cDC2s

# For questions or issues with CITEViz, see:
# GitHub: https://github.com/maxsonBraunLab/CITEViz
# Publication: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05762-1

################################################################################
