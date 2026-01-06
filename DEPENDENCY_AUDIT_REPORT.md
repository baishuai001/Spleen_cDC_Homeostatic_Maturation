# Dependency Audit Report
## Spleen_cDC_Homeostatic_Maturation Project

**Date:** 2026-01-06
**Total Packages Identified:** 29

---

## Executive Summary

This R bioinformatics project uses 29 packages for analyzing Bulk RNA-Seq and CITE-Seq data. The analysis reveals several opportunities for optimization:

- **Redundancy Issues:** Multiple overlapping packages (tidyverse ecosystem)
- **Potential Security Concerns:** Legacy packages requiring updates
- **Bloat Concerns:** Unnecessary package inclusions
- **Missing Best Practices:** No formal dependency management (renv/packrat)

---

## Complete Package Inventory

### Core Bioconductor Packages (9)
1. **edgeR** - Differential expression analysis for RNA-seq
2. **limma** - Linear models for microarray and RNA-seq data
3. **scater** - Single-cell analysis tools
4. **scran** - Single-cell RNA-seq normalization
5. **SingleCellExperiment** - S4 class for single-cell data
6. **sva** - Surrogate variable analysis (batch correction)
7. **dorothea** - Transcription factor activity inference
8. **dyno** - Trajectory inference framework
9. **DoubletFinder** - Doublet detection in single-cell data (GitHub)

### Data Manipulation & Tidyverse (7)
10. **tidyverse** - Meta-package (dplyr, ggplot2, tidyr, readr, tibble, stringr, etc.)
11. **dplyr** - Data manipulation
12. **tidyr** - Data tidying
13. **tibble** - Modern data frames
14. **stringr** - String manipulation
15. **readr** - Fast file reading
16. **sqldf** - SQL queries on data frames

### Visualization Packages (6)
17. **ggplot2** - Grammar of graphics plotting
18. **pheatmap** - Pretty heatmaps
19. **gplots** - Various plotting functions
20. **RColorBrewer** - Color palettes
21. **viridis** - Color scales
22. **rgl** - 3D visualization

### Single-Cell Analysis (1)
23. **Seurat** - Comprehensive single-cell analysis

### Utility Packages (6)
24. **openxlsx** - Excel file I/O
25. **future** - Parallel processing
26. **fields** - Spatial data tools and utilities
27. **modes** - Mode estimation
28. **grid** - Low-level graphics
29. **gridExtra** - Grid graphics extensions

---

## Critical Issues & Recommendations

### 🚨 HIGH PRIORITY: Package Redundancy

**Issue:** Significant redundancy in tidyverse packages

**Current State:**
```r
library("tidyverse")  # Meta-package includes 8 core packages
library("dplyr")      # REDUNDANT - already in tidyverse
library("ggplot2")    # REDUNDANT - already in tidyverse
library("tidyr")      # REDUNDANT - already in tidyverse
library("tibble")     # REDUNDANT - already in tidyverse
library("stringr")    # REDUNDANT - already in tidyverse
library("readr")      # REDUNDANT - already in tidyverse
```

**Recommendation:**
- **Remove individual tidyverse package calls** - Just keep `library("tidyverse")`
- **OR** load only specific packages needed and remove tidyverse
- **Impact:** Reduces redundancy, cleaner code, no functional change

**Estimated Cleanup:**
- Remove 6 redundant library() calls across scripts

---

### 🔧 MEDIUM PRIORITY: Outdated/Deprecated Packages

#### 1. **DoubletFinder** (GitHub Installation)
- **Issue:** Installed from GitHub, not CRAN/Bioconductor
- **Risk:** No official maintenance guarantee, version pinning issues
- **Recommendation:**
  - Consider alternatives like `scDblFinder` (Bioconductor, actively maintained)
  - If keeping DoubletFinder, use renv to pin GitHub commit/version

#### 2. **gplots**
- **Status:** Less actively maintained
- **Modern Alternative:** `pheatmap`, `ComplexHeatmap` (already using pheatmap)
- **Recommendation:** Migrate remaining gplots usage to pheatmap or ggplot2 alternatives

#### 3. **rgl** (3D Visualization)
- **Issue:** Platform-dependent, OpenGL requirements, can cause installation issues
- **Usage:** Only in script 1 for 3D PCA plots
- **Recommendation:**
  - Consider `plotly` for interactive 3D plots (web-based, more portable)
  - Or use static 2D PCA projections (more reproducible)

#### 4. **sqldf**
- **Issue:** Heavy dependency (requires RSQLite, DBI), slower than native dplyr
- **Usage:** Only in one script (file 4)
- **Recommendation:** Replace with dplyr equivalents for better performance

---

### 🛡️ SECURITY & MAINTENANCE

#### Package Update Status

**Known Concerns:**
1. **No dependency management file** (renv.lock, DESCRIPTION)
   - **Risk:** Irreproducible analysis, version conflicts
   - **Recommendation:** Implement `renv` immediately

2. **Bioconductor Version Tracking**
   - **Issue:** No Bioconductor version specified
   - **Risk:** Package compatibility issues across Bioc versions
   - **Recommendation:** Document Bioconductor version (3.18 current as of Jan 2025)

**Security Recommendations:**
```r
# Create renv snapshot
renv::init()
renv::snapshot()

# Document R and Bioconductor versions
R.version
BiocManager::version()
```

---

### 🗑️ BLOAT ANALYSIS

#### Rarely Used Packages (Consider Removing)

1. **modes** - Mode estimation
   - **Usage:** Appears in only 1 script
   - **Alternative:** Simple custom function or built-in alternatives
   - **Recommendation:** Remove and implement simple mode calculation if needed

2. **fields** - Spatial data tools
   - **Usage:** Not clear from scripts what specific function is needed
   - **Issue:** Heavy package with many dependencies
   - **Recommendation:** Replace with specific function from lighter package

3. **grid/gridExtra**
   - **Usage:** Grid graphics arrangement
   - **Alternative:** `patchwork` (modern, intuitive), or `cowplot`
   - **Recommendation:** Consider migration to patchwork

#### Package Size Impact
Large packages that significantly increase installation time/size:
- **Seurat** (~50MB) - ESSENTIAL, cannot remove
- **tidyverse** (~200MB total) - ESSENTIAL for analysis
- **rgl** (~15MB + OpenGL) - OPTIONAL, consider removing
- **fields** (~20MB) - REVIEW usage, possibly removable

---

## Optimization Recommendations

### Phase 1: Immediate Actions (High Impact, Low Risk)

1. **Remove redundant tidyverse package calls**
   ```r
   # BEFORE: 7 separate library calls
   # AFTER: 1 library call
   library("tidyverse")
   ```
   - Files to update: All 9 R scripts
   - Impact: Cleaner code, no functional change

2. **Implement renv for dependency management**
   ```r
   renv::init()
   renv::snapshot()
   ```
   - Creates: `renv.lock` file
   - Benefit: Reproducible environment

3. **Create package documentation**
   - File: `DEPENDENCIES.md`
   - Content: R version, Bioconductor version, package versions
   - Include installation instructions

### Phase 2: Modernization (Medium Impact, Medium Risk)

4. **Replace deprecated packages**
   - `gplots` → `pheatmap` or `ComplexHeatmap`
   - `sqldf` → `dplyr` operations
   - `grid/gridExtra` → `patchwork`

5. **Evaluate DoubletFinder alternatives**
   - Test `scDblFinder` as alternative
   - If keeping DoubletFinder, pin to specific commit

6. **Review rgl usage**
   - Consider `plotly` for 3D plots
   - Or use 2D projections for better reproducibility

### Phase 3: Deep Cleanup (Low Priority, Higher Risk)

7. **Remove low-utility packages**
   - Audit `modes` usage
   - Audit `fields` usage
   - Replace or remove as appropriate

8. **Consolidate visualization approaches**
   - Standardize on ggplot2-based ecosystem
   - Remove legacy plotting packages

---

## Dependency Management Best Practices

### Recommended Setup

Create `.Rprofile` with:
```r
source("renv/activate.R")
options(repos = BiocManager::repositories())
```

Create `DESCRIPTION` file:
```
Package: SpleencDCMaturation
Type: Project
Version: 1.0.0
Depends:
    R (>= 4.2.0)
Imports:
    Seurat (>= 4.3.0),
    tidyverse (>= 2.0.0),
    SingleCellExperiment,
    scater,
    scran,
    edgeR,
    limma,
    dorothea,
    dyno,
    sva,
    openxlsx,
    pheatmap,
    RColorBrewer,
    viridis,
    future
Suggests:
    rgl,
    DoubletFinder
Remotes:
    chris-mcginnis-ucsf/DoubletFinder
```

---

## Security Vulnerability Check

**Method:** No automated R package vulnerability scanner equivalent to npm audit exists, but:

1. **CRAN Package Status:** Check package CRAN status
2. **Bioconductor Support:** Verify packages are actively maintained
3. **GitHub Issues:** Review open security issues

**Current Status (as of Jan 2026):**
- No known critical CVEs in core packages
- All Bioconductor packages appear actively maintained
- DoubletFinder last updated 2019 (⚠️ potential maintenance concern)

**Recommendations:**
- Monitor: https://github.com/ropenscilabs/r-security
- Regular updates: `BiocManager::install()` to update all packages
- Subscribe to R security mailing lists

---

## Implementation Priority Matrix

| Priority | Action | Impact | Effort | Files Affected |
|----------|--------|--------|--------|----------------|
| 🔴 HIGH | Remove tidyverse redundancy | High | Low | All 9 scripts |
| 🔴 HIGH | Implement renv | High | Low | Project root |
| 🟡 MEDIUM | Replace sqldf | Medium | Low | 1 script |
| 🟡 MEDIUM | Document versions | Medium | Low | New file |
| 🟡 MEDIUM | Evaluate DoubletFinder | Medium | Medium | 1 script |
| 🟢 LOW | Replace gplots | Low | Medium | 1-2 scripts |
| 🟢 LOW | Modernize grid/gridExtra | Low | Medium | Multiple |
| 🟢 LOW | Review rgl usage | Low | Medium | 1 script |

---

## Conclusion

The project has a solid foundation with appropriate bioinformatics packages, but suffers from:

1. **Redundancy** - Multiple tidyverse package loads
2. **No formal dependency management** - Missing renv/packrat
3. **Legacy packages** - Some packages have modern alternatives
4. **Potential bloat** - Some rarely-used heavy packages

**Recommended First Steps:**
1. Remove redundant tidyverse library calls (30 min effort)
2. Initialize renv for reproducibility (15 min effort)
3. Document R and package versions (15 min effort)
4. Review and test DoubletFinder alternatives (2-3 hours)

**Expected Benefits:**
- Improved code clarity
- Better reproducibility
- Reduced installation complexity
- Modernized toolchain
- Maintained functionality

---

## Appendix: Package Categories

### Must Keep (Core Functionality)
- Seurat, SingleCellExperiment, scater, scran
- edgeR, limma, sva
- tidyverse (or components)
- dorothea, dyno
- openxlsx (for Excel export)

### Should Review
- DoubletFinder (consider scDblFinder)
- sqldf (replace with dplyr)
- gplots (replace with pheatmap)
- rgl (consider plotly or 2D alternatives)

### Can Potentially Remove
- modes (simple custom function)
- fields (check actual usage)
- grid/gridExtra (use patchwork instead)

### Keep for Visualization
- ggplot2, pheatmap, RColorBrewer, viridis

### Keep for Utilities
- future (parallel processing)
