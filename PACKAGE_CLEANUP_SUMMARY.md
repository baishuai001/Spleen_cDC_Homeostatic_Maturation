# Package Cleanup Summary - Quick Reference

## Immediate Actions (High Priority)

### 1. Remove Unused Package: sqldf

**File:** `2.script_heatmaps_migvsres.R`
**Line:** 4
**Issue:** Package is loaded but never used in the script
**Action:** Delete line 4: `library(sqldf)`
**Impact:** Removes unnecessary dependency
**Risk:** None (not used)

---

### 2. Fix Tidyverse Redundancy

**Problem:** Loading tidyverse AND individual component packages

#### Script: 3.RNA-ADT_HPCscript.R
**Current (Lines 26-36):**
```r
library("tidyverse")  # Line 26
library("dplyr")      # Line 34 - REDUNDANT
library("stringr")    # Line 35 - REDUNDANT
library("readr")      # Line 36 - REDUNDANT
```
**Action:** Remove lines 34-36 (dplyr, stringr, readr already in tidyverse)

#### Script: 4.script_CITEseq_SAM_WT_aggr.R
Check for similar pattern

#### Script: 2.script_heatmaps_migvsres.R
**Current (Lines 3-6):**
```r
library(dplyr)        # Line 3
library(sqldf)        # Line 4 - REMOVE (not used)
library(ggplot2)      # Line 5
```
**Later:**
```r
library(tidyverse)    # Line 149
```
**Action:**
- Remove line 4 (sqldf)
- Option A: Keep tidyverse (line 149), remove dplyr and ggplot2 (lines 3, 5)
- Option B: Remove tidyverse (line 149), keep dplyr and ggplot2

---

## Package Replacement Recommendations

### DoubletFinder → scDblFinder

**Current:** GitHub installation, last updated 2019
**Alternative:** Bioconductor package, actively maintained

**Migration path:**
```r
# Old (DoubletFinder)
library(DoubletFinder)
seuratObj <- runDoubletFinder(seuratObj, dimsToUse, minPCT = 1, maxPCT = 15)

# New (scDblFinder)
library(scDblFinder)
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(seuratObj)
sce <- scDblFinder(sce)
seuratObj$doublet_class <- sce$scDblFinder.class
```

**Benefit:** Better maintained, faster, integrated with Bioconductor

---

## Complete Package List with Status

| Package | Status | Action | Priority |
|---------|--------|--------|----------|
| **Seurat** | ✅ Keep | Essential | - |
| **tidyverse** | ✅ Keep | Consolidate usage | HIGH |
| **dplyr** | ⚠️ Redundant | Remove when tidyverse loaded | HIGH |
| **ggplot2** | ⚠️ Redundant | Remove when tidyverse loaded | HIGH |
| **tidyr** | ⚠️ Redundant | Remove when tidyverse loaded | HIGH |
| **tibble** | ⚠️ Redundant | Remove when tidyverse loaded | HIGH |
| **stringr** | ⚠️ Redundant | Remove when tidyverse loaded | HIGH |
| **readr** | ⚠️ Redundant | Remove when tidyverse loaded | HIGH |
| **sqldf** | ❌ Remove | Not used | HIGH |
| **DoubletFinder** | 🔄 Replace | Use scDblFinder | MEDIUM |
| **gplots** | 🔄 Replace | Use pheatmap (already loaded) | MEDIUM |
| **grid** | ✅ Keep | Needed for pheatmap | - |
| **gridExtra** | 🔄 Consider | Evaluate vs patchwork | LOW |
| **rgl** | ⚠️ Review | Platform-dependent | LOW |
| **modes** | ⚠️ Review | Check usage | LOW |
| **fields** | ⚠️ Review | Check usage | LOW |
| **edgeR** | ✅ Keep | Essential | - |
| **limma** | ✅ Keep | Essential | - |
| **scater** | ✅ Keep | Essential | - |
| **scran** | ✅ Keep | Essential | - |
| **SingleCellExperiment** | ✅ Keep | Essential | - |
| **sva** | ✅ Keep | Batch correction | - |
| **dorothea** | ✅ Keep | TF analysis | - |
| **dyno** | ✅ Keep | Trajectory inference | - |
| **pheatmap** | ✅ Keep | Visualization | - |
| **RColorBrewer** | ✅ Keep | Color palettes | - |
| **viridis** | ✅ Keep | Color scales | - |
| **openxlsx** | ✅ Keep | Excel I/O | - |
| **future** | ✅ Keep | Parallel processing | - |

**Legend:**
- ✅ Keep: Essential package
- ⚠️ Redundant: Loaded multiple times or included in meta-package
- ❌ Remove: Not used
- 🔄 Replace: Better alternative exists

---

## File-by-File Cleanup Checklist

### ✅ 1.script_analysisLimma_getDEgenes.R
**Changes:**
```diff
- library("RColorBrewer")  # Line 350 - Already loaded earlier
- library("gplots")        # Line 351 - Replace with pheatmap
```

### ✅ 2.script_heatmaps_migvsres.R
**Changes:**
```diff
  library(dplyr)
- library(sqldf)           # Line 4 - NOT USED, REMOVE
  library(ggplot2)
```

### ✅ 3.RNA-ADT_HPCscript.R
**Changes:**
```diff
  library("tidyverse")
  library("Seurat")
  library("SingleCellExperiment")
  library("scater")
  library("scran")
  library('DoubletFinder')
  library('fields')
  library('modes')
- library("dplyr")         # REDUNDANT with tidyverse
- library("stringr")       # REDUNDANT with tidyverse
- library("readr")         # REDUNDANT with tidyverse
```

### 📋 Check all other scripts for similar patterns

---

## Quick Start: Minimal Changes

**Step 1:** Create backup
```bash
mkdir backup
cp *.R backup/
```

**Step 2:** Remove sqldf (safe, not used)
```bash
sed -i '/^library(sqldf)/d' 2.script_heatmaps_migvsres.R
```

**Step 3:** Fix tidyverse redundancy in script 3
Edit `3.RNA-ADT_HPCscript.R`:
- Remove individual package loads that are in tidyverse

**Step 4:** Test scripts
```r
source("1.script_analysisLimma_getDEgenes.R")
# etc.
```

---

## Expected Outcomes

### Before Cleanup:
- 29 package library() calls
- Multiple redundant loads
- 1 unused package (sqldf)
- Unclear dependencies

### After Cleanup:
- ~22 unique packages
- No redundant loads
- All packages actively used
- Clear dependency structure
- renv.lock for reproducibility

### Benefits:
- ⚡ Faster script loading
- 📦 Cleaner code
- 🔒 Better reproducibility
- 📝 Clearer documentation
- 🛡️ Easier maintenance

---

## Testing Checklist

After making changes, verify:

```r
# Test each script loads without errors
test_scripts <- c(
  "1.script_analysisLimma_getDEgenes.R",
  "2.script_heatmaps_migvsres.R",
  "3.RNA-ADT_HPCscript.R",
  "4.script_CITEseq_SAM_WT_aggr.R",
  "5.script_CITEseq_SAM_WT_aggr_cDC1_subset.R",
  "6.script_CITEseq_SAM_WT_aggr_cDC2_subset.R",
  "7.script_CITEseq_SAM_WT_aggr_cDC1_subset_TI_ADT.R",
  "8.script_CITEseq_SAM_WT_aggr_cDC1_subset_TI_SCT.R",
  "9.script_CITEseq_SAM_WT_aggr_cDC1_subset_Dorothea.R"
)

for (script in test_scripts) {
  cat("Testing:", script, "\n")
  tryCatch({
    source(script)
    cat("✓ Success\n\n")
  }, error = function(e) {
    cat("✗ Error:", e$message, "\n\n")
  })
}
```

---

## Rollback Instructions

If problems occur:

```bash
# Restore from backup
cp backup/*.R .

# Or use git
git checkout -- *.R
```

---

## Next Steps

1. ✅ Review this summary
2. ✅ Make high-priority changes (redundancy + sqldf)
3. ✅ Test all scripts
4. ✅ Initialize renv: `renv::init()`
5. ✅ Document changes in git commit
6. ⏭️ Plan DoubletFinder evaluation
7. ⏭️ Consider modernizing plotting packages

---

## Questions?

- Which packages am I actually using? → Check DEPENDENCY_AUDIT_REPORT.md
- How do I implement changes? → Check IMPLEMENTATION_GUIDE.md
- What about security? → No known CVEs, but update regularly
- Should I update all packages? → Yes, but test after updating

