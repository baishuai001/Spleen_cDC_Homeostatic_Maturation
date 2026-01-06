# Spleen_cDC_Homeostatic_Maturation

## About

Depending on how an antigen is perceived, dendritic cells (DCs) mature in an immunogenic or tolerogenic manner, safeguarding the balance between immunity and tolerance. Whereas the pathways driving immunogenic maturation in response to infectious insults are well characterized, the signals driving tolerogenic maturation in homeostasis are still poorly understood. Here we demonstrate that engulfment of apoptotic cells triggers homeostatic maturation of conventional cDC1s in the spleen. This process can be modeled by engulfment of empty, non-adjuvanted lipid nanoparticles (LNPs), is marked by intracellular accumulation of cholesterol, and highly unique to type 1 DCs. Engulfment of apoptotic cells or cholesterol-rich LNPs leads to activation of the LXR pathway driving cellular cholesterol efflux and repression of immunogenic genes. In contrast, simultaneous engagement of TLR3 to mimic viral infection via administration of poly(I:C)-adjuvanted LNPs represses the LXR pathway, thus delaying cellular cholesterol efflux and inducing genes that promote T cell immunity.

These data demonstrate how DCs exploit the conserved cellular cholesterol efflux pathway to regulate induction of tolerance or immunity and reveal that administration of non-adjuvanted cholesterol-rich LNPs is a powerful platform for inducing tolerogenic DC maturation.


## Overview scripts

Here's an overview of the various R scripts used in processing the Bulk and CITE-Seq data in the manuscript Bosteels et al.:
- 1.script_analysisLimma_getDEgenes.R: Limma and edgeR workflow for running differential expression analysis on the Bulk RNA-seq data featured in FigS1 and FigS3
- 2.script_heatmaps_migvsres.R: Script for generating the Bulk RNA-seq heatmap between homeostatic mature and immature cDC1s shown in FigS1
- 3.RNA-ADT_HPCscript.R: Standard pipeline used on the High Performance Cluster for analyzing the RNA and ADT assays of the CITE-seq objects
- 4.script_CITEseq_SAM_WT_aggr.R: Script for downstream analysis of the WT_aggr object
- 5.script_CITEseq_SAM_WT_aggr_cDC1_subset.R: Script for downstream analysis of the cDC1 subset of the WT_aggr object
- 6.script_CITEseq_SAM_WT_aggr_cDC2_subset.R: Script for downstream analysis of the cDC2 subset of the WT_aggr object
- 7.script_CITEseq_SAM_WT_aggr_cDC1_subset_TI_ADT.R: Script for trajectory inference on the ADT assay of the cDC1 subset of the WT_aggr object
- 8.script_CITEseq_SAM_WT_aggr_cDC1_subset_TI_SCT.R: Script for trajectory inference on the SCT assay of the cDC1 subset of the WT_aggr object
- 9.script_CITEseq_SAM_WT_aggr_cDC1_subset_Dorothea.R: Script for TF analysis (DoRothEA) of the cDC1 subset of the WT_aggr object

## Dependencies

### Requirements
- R version: >= 4.0.0
- Bioconductor: >= 3.14

### Installation

#### Quick Setup (Recommended)
For reproducible dependency management with renv:
```r
# Run the setup script
source("setup_dependencies.R")
```

#### Manual Installation
If you prefer manual installation:
```r
# Install BiocManager for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("limma", "edgeR", "scater", "scran",
                       "SingleCellExperiment", "dorothea"))

# Install CRAN packages
install.packages(c("Seurat", "tidyverse", "gridExtra", "openxlsx",
                   "RColorBrewer", "viridis", "pheatmap", "future",
                   "dyno", "DoubletFinder"))
```

#### For Collaborators
If an `renv.lock` file is present:
```r
# Install renv
install.packages("renv")

# Restore exact package versions
renv::restore()
```

### Key Packages
- **Seurat**: Single-cell RNA-seq analysis framework
- **limma/edgeR**: Differential expression analysis for bulk RNA-seq
- **scater/scran**: Single-cell QC and normalization
- **dyno**: Trajectory inference
- **dorothea**: Transcription factor activity inference
- **tidyverse**: Data manipulation and visualization

### Dependency Documentation
See `DEPENDENCY_AUDIT.md` for a comprehensive analysis of all dependencies, including security recommendations and cleanup suggestions.

## Citation

Victor Bosteels et al., LXR signaling controls homeostatic dendritic cell maturation. Sci. Immunol. 8, eadd3955(2023). DOI: https://doi.org/10.1126/sciimmunol.add3955


