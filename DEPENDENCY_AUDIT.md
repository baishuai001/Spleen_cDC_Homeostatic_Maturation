# Dependency Audit Report

**Project:** Spleen_cDC_Homeostatic_Maturation
**Date:** 2026-01-06
**Analysis Type:** R Package Dependencies

## Executive Summary

This R bioinformatics project uses 29 distinct packages across 9 scripts for analyzing CITE-seq and bulk RNA-seq data. The analysis identified several issues:

- **Redundant package loading** (same packages loaded multiple times)
- **Overlapping dependencies** (tidyverse conflicts with individual packages)
- **Unnecessary packages** (base R packages loaded explicitly)
- **No dependency management** (missing renv/packrat configuration)
- **Security considerations** (no version pinning or vulnerability tracking)

## Current Dependencies

### Core Analysis Packages
| Package | Purpose | Used In | Notes |
|---------|---------|---------|-------|
| Seurat | Single-cell RNA-seq analysis | Scripts 3-9 | Core package, actively maintained |
| limma | Differential expression | Script 1 | Bioconductor, well-maintained |
| edgeR | Differential expression | Script 1 | Bioconductor, well-maintained |
| dyno | Trajectory inference | Scripts 7-8 | May have stability concerns |
| dorothea | TF activity inference | Script 9 | Part of decoupleR ecosystem |

### Data Processing Packages
| Package | Purpose | Used In | Notes |
|---------|---------|---------|-------|
| tidyverse | Data manipulation suite | Scripts 2-3, 5-8 | Meta-package with 8+ sub-packages |
| dplyr | Data manipulation | All scripts | **Redundant** - included in tidyverse |
| tidyr | Data tidying | Script 9 | **Redundant** - included in tidyverse |
| readr | Data reading | Script 3 | **Redundant** - included in tidyverse |
| stringr | String manipulation | Script 3 | **Redundant** - included in tidyverse |
| tibble | Modern data frames | Script 9 | **Redundant** - included in tidyverse |

### Single-Cell Specific Packages
| Package | Purpose | Used In | Notes |
|---------|---------|---------|-------|
| scater | Single-cell QC | Scripts 3-6, 9 | Well-maintained |
| scran | Single-cell normalization | Script 3 | Well-maintained |
| SingleCellExperiment | Data structure | Script 3 | Bioconductor standard |
| DoubletFinder | Doublet detection | Script 3 | Specialized tool |

### Visualization Packages
| Package | Purpose | Used In | Notes |
|---------|---------|---------|-------|
| ggplot2 | Plotting | Scripts 1-2 | **Redundant** - included in tidyverse |
| pheatmap | Heatmaps | Scripts 2, 9 | Good alternative to base heatmap |
| gridExtra | Grid layouts | Scripts 4-6, 9 | Useful for multi-panel plots |
| grid | Grid graphics | Script 2 | **Unnecessary** - base R package |
| RColorBrewer | Color palettes | Scripts 1-2, 5, 7-8 | Standard for color schemes |
| viridis | Color palettes | Scripts 5-6 | Modern accessible palettes |
| gplots | Advanced plots | Script 1 | Overlaps with ggplot2 |
| rgl | 3D visualization | Script 1 | Heavy dependency, rarely used |

### Utility Packages
| Package | Purpose | Used In | Notes |
|---------|---------|---------|-------|
| openxlsx | Excel I/O | Scripts 1, 4-6, 9 | **Loaded 25+ times** across scripts |
| sqldf | SQL queries on data frames | Script 2 | Can be replaced with dplyr |
| future | Parallel processing | Script 5 | Good for performance |
| sva | Batch effect correction | Script 1 | Specialized, appropriate |
| fields | Spatial statistics | Script 3 | Unclear if necessary |
| modes | Mode calculation | Script 3 | Single function, could inline |

## Critical Issues

### 1. Redundant Package Loading (HIGH PRIORITY)

**Problem:** Packages are loaded multiple times within and across scripts.

**Examples:**
```r
# Script 5.script_CITEseq_SAM_WT_aggr_cDC1_subset.R
library('dplyr')        # Line 4
library('openxlsx')     # Lines 255, 285, 317, 347, 386, 635, 662, 687, 719, 743, 757
library('tidyverse')    # Line 759
library('dplyr')        # Line 758 (again!)
library('RColorBrewer') # Lines 569, 905, 919, 1016, 1051
```

**Impact:**
- Cluttered code
- Potential version conflicts
- Namespace pollution
- Slower script initialization

**Recommendation:**
- Load each package **once** at the script beginning
- Remove all duplicate library() calls within scripts

### 2. Tidyverse Conflicts (HIGH PRIORITY)

**Problem:** Individual tidyverse packages loaded alongside the tidyverse meta-package.

The `tidyverse` package includes:
- ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats

**Examples:**
```r
# Script 3.RNA-ADT_HPCscript.R
library("tidyverse")  # Line 26
library("dplyr")      # Line 34 - REDUNDANT
library("stringr")    # Line 35 - REDUNDANT
library("readr")      # Line 36 - REDUNDANT

# Script 5 loads tidyverse AND dplyr, tibble
# Script 2 loads tidyverse AND dplyr, ggplot2
```

**Impact:**
- Unnecessary memory usage
- Confusing dependency tree
- Mask warnings on package load

**Recommendation:**
- **Option A (Recommended):** Use tidyverse for all scripts needing multiple tidy packages
- **Option B:** Load only specific packages needed (e.g., just dplyr if that's all you need)
- Remove individual package loads when tidyverse is present

### 3. Unnecessary Base R Packages (MEDIUM PRIORITY)

**Problem:** Loading packages that are part of base R installation.

**Example:**
```r
# Script 2.script_heatmaps_migvsres.R, line 119
library('grid')  # grid is part of base R, auto-loaded with many graphics packages
```

**Recommendation:**
- Remove `library('grid')` - it loads automatically when needed

### 4. No Dependency Management System (HIGH PRIORITY)

**Problem:** No `renv.lock`, `DESCRIPTION`, or package version tracking.

**Impact:**
- Cannot reproduce analysis with exact package versions
- No protection against breaking changes in package updates
- Difficult to share/collaborate
- Cannot audit for security vulnerabilities in specific versions

**Recommendation:**
Implement `renv` for dependency management:

```r
# Initialize renv in project root
renv::init()

# This will:
# 1. Scan all R scripts for dependencies
# 2. Create renv.lock with exact versions
# 3. Create project-local library
# 4. Make analysis reproducible
```

### 5. Potentially Unnecessary Packages (MEDIUM PRIORITY)

**Packages to review:**

1. **sqldf** (Script 2)
   - Used for SQL queries on data frames
   - Can be replaced with dplyr (already loaded)
   - Adds Java dependency overhead
   - **Recommendation:** Migrate to dplyr::filter(), dplyr::select(), etc.

2. **fields** (Script 3)
   - Spatial statistics package
   - Large package with many dependencies
   - **Action needed:** Verify if actually used in analysis
   - If only using 1-2 functions, consider alternatives

3. **modes** (Script 3)
   - Provides mode calculation
   - Very small functionality
   - **Recommendation:** Implement custom mode function if rarely used

4. **rgl** (Script 1, line 308)
   - 3D visualization, OpenGL-based
   - Heavy system dependencies
   - **Action needed:** Check if 3D plots are actually generated
   - **Alternative:** plotly for interactive 3D (if needed)

5. **gplots** (Script 1)
   - Overlaps with ggplot2 functionality
   - Consider using ggplot2 consistently

## Security & Maintenance Concerns

### Packages with Known Update Patterns

1. **Seurat** - Major version updates can break code
   - Currently no version pinning
   - Recommendation: Pin to major version (e.g., Seurat ^4.0 or ^5.0)

2. **DoubletFinder** - Less actively maintained
   - Last CRAN update may be outdated
   - Consider alternatives like scDblFinder (Bioconductor)

3. **dyno** - Development has slowed
   - Part of dynverse ecosystem
   - Monitor for updates or consider alternatives

### Bioconductor vs CRAN Mixing

Several Bioconductor packages are mixed with CRAN packages:
- limma, edgeR, scater, scran, SingleCellExperiment (Bioconductor)
- Seurat, tidyverse, etc. (CRAN)

**Recommendation:**
- Document Bioconductor version alongside R version
- Use BiocManager::install() for all Bioconductor packages
- Test compatibility when updating either CRAN or Bioconductor packages

## Performance Considerations

### Heavy Dependencies
1. **tidyverse** - Loads 8+ packages even if only using 1-2
   - **Recommendation:** For production, load only needed packages
   - For analysis scripts, tidyverse is acceptable

2. **rgl** - Requires OpenGL system libraries
   - **Recommendation:** Remove if not essential

### Parallelization
- **future** package is loaded (Script 5, line 693)
- Good for performance with large Seurat objects
- **Recommendation:** Ensure plan(multicore) or plan(multisession) is set appropriately

## Recommendations Summary

### Immediate Actions (High Priority)

1. **Remove duplicate library() calls**
   - Consolidate all package loads to script beginning
   - Remove 25+ redundant openxlsx loads
   - Remove duplicate dplyr, RColorBrewer loads

2. **Resolve tidyverse conflicts**
   - Remove individual package loads when tidyverse is present
   - Choose consistent strategy: all tidyverse OR individual packages

3. **Implement renv dependency management**
   ```r
   # Run once in project root
   install.packages("renv")
   renv::init()
   renv::snapshot()
   ```

4. **Remove unnecessary packages**
   - Remove `library('grid')` from script 2

### Medium Priority Actions

5. **Audit rarely-used packages**
   - Check if fields, modes, rgl, gplots are actually needed
   - Replace sqldf with dplyr equivalents

6. **Create dependency documentation**
   - Document minimum R version required
   - Document Bioconductor version
   - List system dependencies (if any)

7. **Version pinning for critical packages**
   - Pin Seurat to specific major version
   - Pin tidyverse version
   - Document in DESCRIPTION or renv.lock

### Long-term Recommendations

8. **Regular dependency audits**
   - Run `renv::status()` regularly
   - Check for package updates quarterly
   - Test updates in isolated environment before deploying

9. **Security scanning**
   - Monitor Bioconductor/CRAN security advisories
   - Use `oysteR` package to scan for vulnerabilities:
   ```r
   install.packages("oysteR")
   oysteR::audit_installed_r_pkgs()
   ```

10. **Consider lightweight alternatives**
    - data.table instead of dplyr for very large datasets
    - Simpler alternatives to heavy visualization packages

## Implementation Plan

### Phase 1: Cleanup (Estimated effort: 2-4 hours)

1. Create backup branch
2. Refactor each script to load packages only once at top
3. Remove tidyverse/individual package conflicts
4. Remove unnecessary packages
5. Test each script still runs

### Phase 2: Dependency Management (Estimated effort: 1-2 hours)

1. Install and initialize renv
2. Generate renv.lock file
3. Commit renv.lock and .Rprofile to version control
4. Update README with renv instructions

### Phase 3: Optimization (Estimated effort: 4-8 hours)

1. Audit rarely-used packages
2. Replace sqldf with dplyr
3. Test alternatives for heavy dependencies
4. Document system requirements

## Package Replacement Suggestions

| Current Package | Replacement | Reason |
|----------------|-------------|--------|
| sqldf | dplyr | Already loaded, no Java dependency |
| modes | Custom function | Too small to justify dependency |
| rgl | plotly or remove | Lighter, more modern |
| gplots | ggplot2 | Consistency |
| DoubletFinder | scDblFinder | More actively maintained |

## Example: Clean Package Loading Pattern

### Before (Script 5):
```r
library('Seurat')       # Line 3
library('dplyr')        # Line 4
# ... 700 lines ...
library(openxlsx)       # Line 757
library(dplyr)          # Line 758 - duplicate!
library(tidyverse)      # Line 759 - conflicts with dplyr!
# ... more code ...
library(openxlsx)       # Line 757 - duplicate again!
```

### After:
```r
# Load packages once at beginning
library('Seurat')
library('tidyverse')    # Includes dplyr, tidyr, ggplot2, etc.
library('gridExtra')
library('scater')
library('openxlsx')
library('RColorBrewer')
library('viridis')
library('future')

# Set future plan for parallel processing
plan(multicore, workers = 4)

# Rest of script...
```

## Testing Checklist

After implementing changes:
- [ ] All scripts run without errors
- [ ] No "masked object" warnings for critical functions
- [ ] Outputs match previous results (statistical reproducibility)
- [ ] renv::status() shows all dependencies captured
- [ ] renv.lock file committed to git
- [ ] README updated with dependency installation instructions
- [ ] Collaborators can reproduce environment with renv::restore()

## Resources

- [renv documentation](https://rstudio.github.io/renv/)
- [Bioconductor installation guide](https://bioconductor.org/install/)
- [tidyverse style guide](https://style.tidyverse.org/)
- [R package security with oysteR](https://github.com/sonatype-nexus-community/oysteR)

## Conclusion

This project has a robust set of dependencies for bioinformatics analysis, but suffers from:
1. **Redundant package loading** causing code bloat
2. **No version control** for dependencies
3. **tidyverse conflicts** with individual packages
4. **Some potentially unnecessary packages**

Implementing the recommended changes will:
- Improve code clarity and maintainability
- Enable reproducible research
- Reduce security risks through version tracking
- Potentially improve performance by removing unnecessary packages

**Priority:** Implement renv and remove duplicate library() calls as first steps.
