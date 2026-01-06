# Dependency Cleanup Implementation Guide

This guide provides step-by-step instructions for implementing the recommendations from the dependency audit.

---

## Quick Wins (Do These First)

### 1. Initialize renv for Dependency Management

**Why:** Ensures reproducibility and version control for all R packages

**Steps:**
```r
# In R console, from project root:
install.packages("renv")
renv::init()
renv::snapshot()
```

**Result:** Creates `renv.lock` file with all package versions

**Time:** 5-10 minutes

---

### 2. Remove Redundant tidyverse Package Calls

**Problem:** Loading tidyverse AND its component packages is redundant

**Current pattern in scripts:**
```r
library("tidyverse")
library("dplyr")      # REMOVE - included in tidyverse
library("ggplot2")    # REMOVE - included in tidyverse
library("tidyr")      # REMOVE - included in tidyverse
library("tibble")     # REMOVE - included in tidyverse
library("stringr")    # REMOVE - included in tidyverse
library("readr")      # REMOVE - included in tidyverse
```

**Recommended approach:**
```r
# Option 1: Use tidyverse (loads all at once)
library("tidyverse")

# Option 2: Load only what you need (more explicit)
library("dplyr")
library("ggplot2")
library("tidyr")
# etc.
```

**Scripts to update:**
- 3.RNA-ADT_HPCscript.R
- 4.script_CITEseq_SAM_WT_aggr.R
- All other scripts using tidyverse

**Time:** 15-20 minutes

---

### 3. Document R Environment

**Create file:** `R_ENVIRONMENT.md`

**Content:**
```markdown
# R Environment Information

## R Version
R version 4.x.x (update with your version)

## Bioconductor Version
Bioconductor version 3.18 (or your version)

## System Information
- Platform: [Linux/macOS/Windows]
- Date created: 2026-01-06

## Installation Instructions

### Step 1: Install R
Download R from: https://cran.r-project.org/

### Step 2: Install BiocManager
```r
install.packages("BiocManager")
BiocManager::install(version = "3.18")
```

### Step 3: Restore Package Environment
```r
install.packages("renv")
renv::restore()
```

### GitHub Packages
DoubletFinder requires manual installation:
```r
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
```

## Verification
```r
# Check R version
R.version.string

# Check Bioconductor version
BiocManager::version()

# Check installed packages
renv::status()
```
```

**Time:** 10 minutes

---

## Medium-Priority Updates

### 4. Replace sqldf with dplyr

**Problem:** sqldf is slower and heavier than dplyr

**Find usage:**
```bash
grep -n "sqldf" *.R
```

**Typical replacements:**

**Before (sqldf):**
```r
library(sqldf)
result <- sqldf("SELECT * FROM df1
                 INNER JOIN df2
                 ON df1.id = df2.id
                 WHERE df1.value > 10")
```

**After (dplyr):**
```r
library(dplyr)
result <- df1 %>%
  inner_join(df2, by = "id") %>%
  filter(value > 10)
```

**Time:** 30-60 minutes depending on complexity

---

### 5. Review DoubletFinder vs scDblFinder

**Current:** DoubletFinder (GitHub, last updated 2019)
**Alternative:** scDblFinder (Bioconductor, actively maintained)

**Test scDblFinder:**
```r
# Install
BiocManager::install("scDblFinder")

# Basic usage (comparable to DoubletFinder)
library(scDblFinder)
sce <- scDblFinder(sce)

# Results in sce$scDblFinder.class
table(sce$scDblFinder.class)
```

**Comparison:**
| Feature | DoubletFinder | scDblFinder |
|---------|---------------|-------------|
| Maintenance | Inactive | Active |
| Repository | GitHub | Bioconductor |
| Documentation | Good | Excellent |
| Performance | Good | Better |
| Integration | Seurat | SCE/Seurat |

**Recommendation:** Test scDblFinder on subset of data, compare results

**Time:** 2-3 hours for testing and validation

---

### 6. Modernize Plotting: gplots → pheatmap

**Problem:** gplots is less actively maintained

**Find usage:**
```bash
grep -n "gplots\|heatmap.2" *.R
```

**Current (gplots):**
```r
library("gplots")
heatmap.2(as.matrix(distsRL),
          Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=rev(hmcol),margin=c(13, 13))
```

**Alternative (pheatmap - already in project):**
```r
library("pheatmap")
pheatmap(as.matrix(distsRL),
         clustering_method = "ward.D",
         color = rev(hmcol),
         show_colnames = TRUE,
         show_rownames = TRUE)
```

**Files affected:**
- 1.script_analysisLimma_getDEgenes.R

**Time:** 30 minutes

---

### 7. Replace grid/gridExtra with patchwork

**Problem:** gridExtra is older, patchwork is more intuitive

**Install patchwork:**
```r
install.packages("patchwork")
```

**Before (gridExtra):**
```r
library(gridExtra)
grid.arrange(plot1, plot2, ncol = 2)
```

**After (patchwork):**
```r
library(patchwork)
plot1 + plot2  # Side by side
plot1 / plot2  # Stacked
(plot1 | plot2) / plot3  # Complex layouts
```

**Time:** 1 hour

---

## Lower-Priority Optimizations

### 8. Evaluate rgl Usage

**Current usage:** 3D PCA plots in script 1

**Issue:** Platform-dependent, OpenGL requirements

**Alternative 1: plotly (interactive)**
```r
library(plotly)
plot_ly(data = pca_df,
        x = ~PC1, y = ~PC2, z = ~PC3,
        color = ~group,
        type = "scatter3d",
        mode = "markers")
```

**Alternative 2: Static 2D projections**
```r
library(ggplot2)

# PC1 vs PC2
p1 <- ggplot(pca_df, aes(PC1, PC2, color = group)) + geom_point()

# PC1 vs PC3
p2 <- ggplot(pca_df, aes(PC1, PC3, color = group)) + geom_point()

# PC2 vs PC3
p3 <- ggplot(pca_df, aes(PC2, PC3, color = group)) + geom_point()

library(patchwork)
p1 + p2 + p3
```

**Recommendation:** Keep rgl for now unless installation issues arise

**Time:** 2-3 hours if replacing

---

### 9. Audit Rarely-Used Packages

#### Package: modes
**Usage check:**
```bash
grep -n "modes" *.R
```

**If minimal usage:** Replace with simple function
```r
# Instead of library(modes)
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
```

#### Package: fields
**Usage check:**
```bash
grep -n "fields" *.R | grep -v "^#"
```

**Determine:** Which specific function is being used?
- If only using one function, consider lighter alternatives

**Time:** 1-2 hours for thorough audit

---

## Cleanup Script Templates

### Template 1: Standard Package Loading Block

**Create:** `R/setup_packages.R`

```r
# Standard package loading for analysis scripts
# Source this at the beginning of each script

# Core Bioconductor packages
suppressPackageStartupMessages({
  library("SingleCellExperiment")
  library("scater")
  library("scran")
  library("edgeR")
  library("limma")
  library("sva")
})

# Single-cell analysis
suppressPackageStartupMessages({
  library("Seurat")
  library("DoubletFinder")
})

# Data manipulation (tidyverse)
suppressPackageStartupMessages({
  library("tidyverse")  # Includes dplyr, ggplot2, tidyr, readr, etc.
})

# Visualization
suppressPackageStartupMessages({
  library("pheatmap")
  library("RColorBrewer")
  library("viridis")
})

# Utilities
suppressPackageStartupMessages({
  library("openxlsx")
  library("future")
})

# Specialized packages (load as needed)
# library("dorothea")  # TF analysis
# library("dyno")      # Trajectory inference
# library("rgl")       # 3D visualization
```

**Usage in scripts:**
```r
source("R/setup_packages.R")
# Now proceed with analysis
```

---

### Template 2: Package Installation Script

**Create:** `install_dependencies.R`

```r
#!/usr/bin/env Rscript

# Package installation script for Spleen cDC Maturation project
# Run this once to set up your R environment

cat("Installing R package dependencies...\n\n")

# Install renv for dependency management
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Set Bioconductor version
BiocManager::install(version = "3.18", ask = FALSE)

# Core CRAN packages
cran_packages <- c(
  "tidyverse",
  "openxlsx",
  "pheatmap",
  "RColorBrewer",
  "viridis",
  "rgl",
  "future",
  "devtools"
)

cat("Installing CRAN packages...\n")
install.packages(cran_packages)

# Bioconductor packages
bioc_packages <- c(
  "SingleCellExperiment",
  "scater",
  "scran",
  "edgeR",
  "limma",
  "sva",
  "dorothea",
  "dyno"
)

cat("\nInstalling Bioconductor packages...\n")
BiocManager::install(bioc_packages)

# Seurat (from CRAN)
cat("\nInstalling Seurat...\n")
install.packages("Seurat")

# GitHub packages
cat("\nInstalling GitHub packages...\n")
devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")

# Initialize renv
cat("\nInitializing renv...\n")
renv::init()

cat("\n✓ Installation complete!\n")
cat("Run renv::snapshot() to save the current state\n")
```

**Usage:**
```bash
Rscript install_dependencies.R
```

---

## Verification Checklist

After making changes, verify everything still works:

- [ ] All scripts run without errors
- [ ] Package loading is successful
- [ ] No deprecated function warnings
- [ ] Results match previous analysis (use saved .rds files for comparison)
- [ ] renv.lock is up to date
- [ ] Documentation is complete

**Verification command:**
```r
# Check for missing packages
renv::status()

# Test script execution
source("1.script_analysisLimma_getDEgenes.R")
```

---

## Rollback Plan

If issues arise:

1. **Restore from renv snapshot:**
```r
renv::restore()
```

2. **Revert specific package:**
```r
renv::install("package@version")
```

3. **Git revert (if using version control):**
```bash
git checkout HEAD~1 -- script.R
```

---

## Timeline Estimate

| Phase | Tasks | Time Estimate |
|-------|-------|---------------|
| **Phase 1** | renv + redundancy cleanup | 1-2 hours |
| **Phase 2** | sqldf, gplots, gridExtra | 3-4 hours |
| **Phase 3** | DoubletFinder testing | 2-3 hours |
| **Phase 4** | Deep cleanup (modes, fields) | 2-3 hours |
| **Phase 5** | Testing & verification | 2-3 hours |
| **Total** | | 10-15 hours |

---

## Questions or Issues?

If you encounter problems:

1. Check package documentation: `?package_name`
2. Check vignettes: `browseVignettes("package_name")`
3. Search issues: GitHub repos for each package
4. Bioconductor support: https://support.bioconductor.org/

---

## Success Metrics

After implementation, you should have:

✓ `renv.lock` file with all dependencies
✓ Cleaner, more maintainable code
✓ Faster package loading times
✓ Better reproducibility
✓ Modern, actively-maintained packages
✓ Complete documentation

