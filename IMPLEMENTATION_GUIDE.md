# Dependency Cleanup Implementation Guide

This guide provides specific commands and code changes to implement the recommendations from DEPENDENCY_AUDIT.md.

## Quick Start: Initialize renv

Run these commands in R console from the project root:

```r
# Install renv if not already installed
install.packages("renv")

# Initialize renv for this project
renv::init()

# Take a snapshot of current dependencies
renv::snapshot()
```

This creates:
- `renv.lock` - Lockfile with exact package versions
- `renv/` - Project library directory
- `.Rprofile` - Auto-loads renv on project open

## Script-by-Script Cleanup

### Script 1: 1.script_analysisLimma_getDEgenes.R

**Current package loads (scattered throughout):**
```r
library("edgeR")        # Line 5
library("limma")        # Line 6
library("ggplot2")      # Line 7
library("dplyr")        # Line 8
library('sva')          # Line 160
library("rgl")          # Line 308
library("RColorBrewer") # Line 350
library("gplots")       # Line 351
library('openxlsx')     # Line 462
```

**Recommended cleanup:**
Replace lines 5-8 with this consolidated block, and remove all other library() calls:

```r
# Load all required packages at script start
library("edgeR")
library("limma")
library("ggplot2")
library("dplyr")
library("sva")
library("RColorBrewer")
library("openxlsx")
# library("rgl")     # Uncomment only if 3D plots are needed
# library("gplots")  # Consider migrating to ggplot2
```

**Lines to delete:**
- Line 160: `library('sva')`
- Line 308: `library("rgl")` (or move to top with comment)
- Line 350: `library("RColorBrewer")`
- Line 351: `library("gplots")`
- Line 462: `library('openxlsx')`

---

### Script 2: 2.script_heatmaps_migvsres.R

**Current package loads:**
```r
library(dplyr)          # Line 3
library(sqldf)          # Line 4
library(ggplot2)        # Line 5
library('pheatmap')     # Line 118
library('grid')         # Line 119 - UNNECESSARY (base R)
library(tidyverse)      # Line 149 - CONFLICTS with lines 3, 5
library(RColorBrewer)   # Line 196
```

**Recommended cleanup:**
Replace lines 3-5 with:

```r
# Load all required packages at script start
library(tidyverse)      # Includes dplyr, ggplot2, and more
library(pheatmap)
library(RColorBrewer)
# library(sqldf)        # TODO: Migrate SQL queries to dplyr
```

**Lines to delete:**
- Line 3: `library(dplyr)` - included in tidyverse
- Line 4: `library(sqldf)` - recommend migrating to dplyr
- Line 5: `library(ggplot2)` - included in tidyverse
- Line 118: `library('pheatmap')` - moved to top
- Line 119: `library('grid')` - not needed
- Line 149: `library(tidyverse)` - moved to top
- Line 196: `library(RColorBrewer)` - moved to top

---

### Script 3: 3.RNA-ADT_HPCscript.R

**Current package loads:**
```r
library("tidyverse")          # Line 26
library("Seurat")             # Line 27
library("SingleCellExperiment") # Line 28
library("scater")             # Line 29
library("scran")              # Line 30
library('DoubletFinder')      # Line 31
library('fields')             # Line 32 - VERIFY IF NEEDED
library('modes')              # Line 33 - CONSIDER REMOVING
library("dplyr")              # Line 34 - REDUNDANT (in tidyverse)
library("stringr")            # Line 35 - REDUNDANT (in tidyverse)
library("readr")              # Line 36 - REDUNDANT (in tidyverse)
```

**Recommended cleanup:**
Replace lines 26-36 with:

```r
# Load all required packages at script start
library("tidyverse")           # Includes dplyr, stringr, readr, etc.
library("Seurat")
library("SingleCellExperiment")
library("scater")
library("scran")
library('DoubletFinder')
# library('fields')            # TODO: Verify if actually used in analysis
# library('modes')             # TODO: Replace with custom mode function if needed
```

**Lines to delete:**
- Line 32: `library('fields')` - verify necessity
- Line 33: `library('modes')` - consider custom function
- Line 34: `library("dplyr")` - redundant
- Line 35: `library("stringr")` - redundant
- Line 36: `library("readr")` - redundant

---

### Script 4: 4.script_CITEseq_SAM_WT_aggr.R

**Current package loads:**
```r
library('Seurat')     # Line 3
library('dplyr')      # Line 4
library('gridExtra')  # Line 5
library('scater')     # Line 6
library('openxlsx')   # Lines 177, 207, 239, 269, 299
```

**Recommended cleanup:**
Replace lines 3-6 with:

```r
# Load all required packages at script start
library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('openxlsx')
```

**Lines to delete:**
- Line 177: `library('openxlsx')`
- Line 207: `library('openxlsx')`
- Line 239: `library('openxlsx')`
- Line 269: `library('openxlsx')`
- Line 299: `library('openxlsx')`

---

### Script 5: 5.script_CITEseq_SAM_WT_aggr_cDC1_subset.R (MOST ISSUES)

**Current package loads:**
```r
library('Seurat')          # Line 3
library('dplyr')           # Line 4
library('gridExtra')       # Line 5
library('scater')          # Line 6
library('openxlsx')        # Lines 255, 285, 317, 347, 386, 635, 662, 687, 719, 743, 757
library("RColorBrewer")    # Lines 569, 905, 919, 1016, 1051
library(future)            # Line 693
library(openxlsx)          # Line 757 (again!)
library(dplyr)             # Line 758 - DUPLICATE!
library(tidyverse)         # Line 759 - CONFLICTS!
library(viridis)           # Line 850
library("RColorBrewer")    # Lines 905, 919, 1016, 1051 - DUPLICATES!
```

**Recommended cleanup:**
Replace lines 3-6 with:

```r
# Load all required packages at script start
library('Seurat')
library('tidyverse')       # Includes dplyr and many others
library('gridExtra')
library('scater')
library('openxlsx')
library('RColorBrewer')
library('viridis')
library('future')

# Configure parallel processing
plan(multicore, workers = 4)  # Adjust workers as needed
```

**Lines to delete (ALL of these):**
- Line 4: `library('dplyr')` - included in tidyverse
- Lines 255, 285, 317, 347, 386, 635, 662, 687, 719, 743: All `library('openxlsx')`
- Lines 569, 905, 919, 1016, 1051: All `library("RColorBrewer")`
- Line 693: `library(future)` - moved to top
- Lines 757-759: `library(openxlsx)`, `library(dplyr)`, `library(tidyverse)`
- Line 850: `library(viridis)` - moved to top

**This removes 25+ redundant library calls from this script alone!**

---

### Script 6: 6.script_CITEseq_SAM_WT_aggr_cDC2_subset.R

**Current package loads:**
```r
library('Seurat')     # Line 3
library('dplyr')      # Line 4
library('gridExtra')  # Line 5
library('scater')     # Line 6
library('openxlsx')   # Lines 171, 201, 233, 263, 313, 398
library(viridis)      # Line 359
library(dplyr)        # Line 399 - DUPLICATE!
library(tidyverse)    # Line 400 - CONFLICTS!
```

**Recommended cleanup:**
Replace lines 3-6 with:

```r
# Load all required packages at script start
library('Seurat')
library('tidyverse')       # Includes dplyr
library('gridExtra')
library('scater')
library('openxlsx')
library('viridis')
```

**Lines to delete:**
- Line 4: `library('dplyr')` - included in tidyverse
- Lines 171, 201, 233, 263, 313, 398: All `library('openxlsx')`
- Line 359: `library(viridis)`
- Lines 399-400: `library(dplyr)`, `library(tidyverse)`

---

### Script 7: 7.script_CITEseq_SAM_WT_aggr_cDC1_subset_TI_ADT.R

**Current package loads:**
```r
library(dyno)          # Line 3
library(tidyverse)     # Line 4
library(Seurat)        # Line 5
library(RColorBrewer)  # Line 138
```

**Recommended cleanup:**
Replace lines 3-5 with:

```r
# Load all required packages at script start
library(dyno)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
```

**Lines to delete:**
- Line 138: `library(RColorBrewer)` - moved to top

---

### Script 8: 8.script_CITEseq_SAM_WT_aggr_cDC1_subset_TI_SCT.R

**Current package loads:**
```r
library(dyno)          # Line 3
library(tidyverse)     # Line 4
library(Seurat)        # Line 5
library(RColorBrewer)  # Line 137
```

**Recommended cleanup:**
Replace lines 3-5 with:

```r
# Load all required packages at script start
library(dyno)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
```

**Lines to delete:**
- Line 137: `library(RColorBrewer)` - moved to top

---

### Script 9: 9.script_CITEseq_SAM_WT_aggr_cDC1_subset_Dorothea.R

**Current package loads:**
```r
library('Seurat')     # Line 3
library('dplyr')      # Line 4
library('gridExtra')  # Line 5
library('scater')     # Line 6
library('dorothea')   # Line 7
library('tidyr')      # Line 8
library('pheatmap')   # Line 9
library('tibble')     # Line 10
library('openxlsx')   # Lines 129, 151
```

**Recommended cleanup:**
Replace lines 3-10 with:

```r
# Load all required packages at script start
library('Seurat')
library('tidyverse')       # Includes dplyr, tidyr, tibble
library('gridExtra')
library('scater')
library('dorothea')
library('pheatmap')
library('openxlsx')
```

**Lines to delete:**
- Line 4: `library('dplyr')` - included in tidyverse
- Line 8: `library('tidyr')` - included in tidyverse
- Line 10: `library('tibble')` - included in tidyverse
- Lines 129, 151: `library('openxlsx')`

---

## Summary of Changes

| Script | Redundant Library Calls Removed | Key Issue |
|--------|--------------------------------|-----------|
| Script 1 | 5 calls | Scattered throughout file |
| Script 2 | 6 calls | grid + tidyverse conflicts |
| Script 3 | 3 calls | tidyverse sub-packages loaded individually |
| Script 4 | 5 calls | openxlsx loaded 5 times |
| Script 5 | **25+ calls** | openxlsx loaded 11 times, multiple duplicates |
| Script 6 | 8 calls | tidyverse + dplyr conflict |
| Script 7 | 1 call | Minor cleanup |
| Script 8 | 1 call | Minor cleanup |
| Script 9 | 5 calls | tidyverse sub-packages |

**Total: ~60 redundant library() calls removed**

---

## After Cleanup: Create renv Lockfile

Once all scripts are cleaned up:

```r
# In R console, from project root:

# Initialize renv (if not done already)
renv::init()

# Check status
renv::status()

# Create snapshot
renv::snapshot()

# Verify lockfile created
file.exists("renv.lock")
```

---

## Git Commit Strategy

```bash
# Stage the changes
git add DEPENDENCY_AUDIT.md IMPLEMENTATION_GUIDE.md

# Commit documentation
git commit -m "Add dependency audit and implementation guide

- Identified 60+ redundant library() calls across 9 scripts
- Found tidyverse conflicts with individual package loads
- Recommended renv for reproducible dependency management
- Provided script-by-script cleanup instructions"

# After implementing cleanup in scripts:
git add *.R renv.lock .Rprofile renv/

git commit -m "Refactor: Clean up package dependencies

- Consolidate all library() calls to script beginnings
- Remove 60+ redundant library() calls
- Resolve tidyverse conflicts (dplyr, ggplot2, tidyr, etc.)
- Remove unnecessary base R package loads (grid)
- Initialize renv for reproducible dependency management

Scripts modified:
- All 9 R scripts: Move library() calls to top, remove duplicates
- Initialize renv with renv.lock for version tracking"
```

---

## Validation Checklist

After implementing changes:

- [ ] Each script loads packages only once at the beginning
- [ ] No tidyverse + individual package conflicts
- [ ] No duplicate library() calls within any script
- [ ] `library('grid')` removed from Script 2
- [ ] All scripts run without errors
- [ ] Output files match previous results (reproducibility check)
- [ ] `renv.lock` file created and committed
- [ ] `.Rprofile` committed
- [ ] `renv/settings.dcf` committed
- [ ] `renv/activate.R` committed (usually)
- [ ] README updated with renv usage instructions

---

## Future Maintenance

### Update Dependencies
```r
# Check for outdated packages
renv::status()

# Update specific package
renv::update("Seurat")

# Update all packages
renv::update()

# Snapshot changes
renv::snapshot()
```

### Audit for Vulnerabilities
```r
# Install security auditing package
install.packages("oysteR")

# Audit installed packages
oysteR::audit_installed_r_pkgs()
```

### Share with Collaborators
```r
# Collaborators clone repo and run:
renv::restore()

# This installs exact package versions from renv.lock
```

---

## Questions & Troubleshooting

### Q: Script fails after removing library() call
**A:** The package may be used later in the script. Move library() to the top section instead of deleting entirely.

### Q: "Object masked" warnings after cleanup
**A:** This is normal when packages have overlapping function names. Usually safe to ignore, but verify functions work as expected.

### Q: renv::init() hangs or fails
**A:**
- Ensure all scripts are syntax-valid R code
- Check internet connection for package downloads
- Try `renv::init(bare = TRUE)` then manually add packages

### Q: Large renv/ directory
**A:**
- Add `renv/library/` to `.gitignore` (usually done automatically)
- Only commit `renv.lock`, `.Rprofile`, and `renv/settings.dcf`
- The actual packages are restored via `renv::restore()`, not committed

---

## Alternative: Lightweight Cleanup (No renv)

If you want to clean up without implementing renv:

1. Apply the script-by-script changes above
2. Skip the renv initialization
3. Document required packages in README.md
4. Manually track working package versions in comments

**Example README section:**
```markdown
## Dependencies

Required R version: >= 4.0.0
Required Bioconductor version: >= 3.14

### CRAN Packages
- Seurat (tested with v4.3.0)
- tidyverse (tested with v2.0.0)
- [... etc ...]

### Bioconductor Packages
- limma
- edgeR
- scater
- scran

### Installation
Install CRAN packages:
install.packages(c("Seurat", "tidyverse", ...))

Install Bioconductor packages:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "scater", "scran"))
```

This is less robust than renv but better than nothing.
