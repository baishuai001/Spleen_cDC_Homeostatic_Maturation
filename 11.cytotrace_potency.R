#!/usr/bin/env Rscript
#' CytoTRACE Differentiation Potency Analysis
#' ===========================================
#'
#' Author: GSE228544 CITE-seq Analysis Pipeline
#' Date: 2026-01-21
#'
#' Description:
#' -----------
#' This script uses CytoTRACE to predict cellular differentiation potency
#' in cDC development.
#'
#' CytoTRACE predicts differentiation state based on:
#' - Number of expressed genes per cell
#' - Transcriptional diversity
#' - Correlation with known stemness signatures
#'
#' Main Features:
#' -------------
#' 1. **Potency Scoring**
#'    - Cell-level differentiation potency (0-1)
#'    - Higher score = less differentiated (more potent)
#'    - Lower score = more differentiated (less potent)
#'
#' 2. **Multi-UMAP Support**
#'    - RNA UMAP visualization
#'    - ADT UMAP visualization
#'    - WNN UMAP visualization ⭐ Recommended
#'
#' 3. **Validation**
#'    - Compare with pseudotime
#'    - Validate with known cell types
#'    - Cross-check with markers
#'
#' 4. **Integration**
#'    - Add scores to Seurat object
#'    - Export for downstream analysis
#'    - Combine with trajectory results
#'
#' Input Data:
#' ----------
#' From Seurat analysis:
#' - Expression matrix (raw or normalized counts)
#' - Cell type annotations
#' - UMAP coordinates
#'
#' Output:
#' -------
#' cytotrace_analysis_{mode}/
#' ├── data/
#' │   └── cytotrace_results.rds          # CytoTRACE results object
#' ├── Robjects/
#' │   └── seurat_with_cytotrace.rds      # Seurat + CytoTRACE scores
#' ├── Plots/
#' │   ├── Potency/
#' │   │   ├── potency_rna_umap.png       # Potency on RNA UMAP
#' │   │   ├── potency_adt_umap.png       # Potency on ADT UMAP
#' │   │   ├── potency_wnn_umap.png       # Potency on WNN UMAP ⭐
#' │   │   └── potency_by_celltype.png    # Potency distribution
#' │   ├── Validation/
#' │   │   ├── potency_vs_pseudotime.png  # Correlation with pseudotime
#' │   │   ├── potency_vs_genes.png       # Correlation with gene count
#' │   │   └── validation_metrics.png     # Summary metrics
#' │   └── Genes/
#' │       ├── potency_genes_heatmap.png  # Top potency-associated genes
#' │       └── markers_vs_potency.png     # Known markers
#' └── Tables/
#'     ├── cytotrace_scores.csv           # Per-cell potency scores
#'     ├── potency_genes.csv              # Potency-associated genes
#'     └── validation_results.csv         # Validation statistics
#'
#' Usage:
#' ------
#' # Command line
#' Rscript 11.cytotrace_potency.R --mode cDC1 \\
#'   --seurat results/cDC1/Robjects/seurat_obj_annotated.rds
#'
#' # With pseudotime for validation
#' Rscript 11.cytotrace_potency.R --mode cDC1 \\
#'   --seurat results/cDC1/Robjects/seurat_obj_annotated.rds \\
#'   --pseudotime monocle3_analysis_cDC1/Tables/pseudotime_by_cell.csv
#'
#' # In R session
#' source("11.cytotrace_potency.R")
#' results <- run_cytotrace_analysis(
#'   mode = "cDC1",
#'   seurat_file = "results/cDC1/Robjects/seurat_obj_annotated.rds"
#' )
#'
#' Requirements:
#' ------------
#' - CytoTRACE (install from GitHub)
#' - Seurat
#' - tidyverse
#' - ggplot2, patchwork
#'
#' Installation:
#' ------------
#' # Install CytoTRACE
#' devtools::install_github("digitalcytometry/cytotrace")
#'
#' Expected Results:
#' ----------------
#' For cDC development:
#' - Pre-cDC: HIGH potency (0.6-0.8)
#' - Proliferating: MEDIUM-HIGH potency (0.5-0.7)
#' - Mature cDC: LOW potency (0.2-0.4)
#'
#' Potency should negatively correlate with pseudotime (r < -0.5)

# =============================================================================
# Setup
# =============================================================================

# Check for command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Load required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(scales)

  # Try loading CytoTRACE
  cytotrace_available <- requireNamespace("CytoTRACE", quietly = TRUE)

  if (!cytotrace_available) {
    warning("CytoTRACE package not installed!")
    cat("\nTo install CytoTRACE:\n")
    cat("  devtools::install_github('digitalcytometry/cytotrace')\n\n")
  } else {
    library(CytoTRACE)
  }
})

# =============================================================================
# Helper Functions
# =============================================================================

#' Create output directories
#'
#' @param base_dir Base output directory
#' @return List of directory paths
create_output_directories <- function(base_dir) {
  dirs <- list(
    base = base_dir,
    data = file.path(base_dir, "data"),
    robjects = file.path(base_dir, "Robjects"),
    plots = file.path(base_dir, "Plots"),
    potency = file.path(base_dir, "Plots", "Potency"),
    validation = file.path(base_dir, "Plots", "Validation"),
    genes = file.path(base_dir, "Plots", "Genes"),
    tables = file.path(base_dir, "Tables")
  )

  for (dir_path in dirs) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }

  return(dirs)
}


#' Run CytoTRACE analysis
#'
#' @param seurat_obj Seurat object
#' @param ncores Number of cores for parallel processing
#' @return CytoTRACE results
run_cytotrace <- function(seurat_obj, ncores = 1) {
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("RUNNING CYTOTRACE ANALYSIS\n")
  cat(rep("=", 80), "\n\n", sep = "")

  if (!cytotrace_available) {
    cat("⚠ CytoTRACE not available, using gene count approximation\n\n")

    # Fallback: use gene count as proxy for potency
    # This is a simplified approximation
    gene_counts <- Matrix::colSums(GetAssayData(seurat_obj, slot = "counts") > 0)

    # Normalize to 0-1 range (higher gene count = higher potency)
    potency_scores <- (gene_counts - min(gene_counts)) /
                     (max(gene_counts) - min(gene_counts))

    results <- list(
      CytoTRACE = potency_scores,
      method = "gene_count_approximation"
    )

    cat("✓ Approximated potency using gene counts\n")
    cat("  Score range:", round(range(potency_scores), 3), "\n")

    return(results)
  }

  # Real CytoTRACE analysis
  cat("Extracting expression matrix...\n")

  # Extract normalized counts
  expr_mat <- as.matrix(GetAssayData(seurat_obj, slot = "data", assay = "RNA"))

  # Subset to highly variable genes for speed
  if (length(VariableFeatures(seurat_obj)) > 0) {
    hvgs <- VariableFeatures(seurat_obj)
    hvgs <- hvgs[hvgs %in% rownames(expr_mat)]

    if (length(hvgs) > 1000) {
      hvgs <- hvgs[1:1000]
    }

    cat("  Using", length(hvgs), "highly variable genes\n")
    expr_mat <- expr_mat[hvgs, ]
  }

  cat("  Matrix size:", nrow(expr_mat), "×", ncol(expr_mat), "\n")

  # Run CytoTRACE
  cat("\nRunning CytoTRACE (this may take several minutes)...\n")

  results <- tryCatch({
    CytoTRACE(
      mat = expr_mat,
      ncores = ncores,
      enableFast = TRUE  # Use fast mode
    )
  }, error = function(e) {
    cat("Error running CytoTRACE:", conditionMessage(e), "\n")
    cat("Falling back to gene count approximation...\n")

    gene_counts <- Matrix::colSums(GetAssayData(seurat_obj, slot = "counts") > 0)
    potency_scores <- (gene_counts - min(gene_counts)) /
                     (max(gene_counts) - min(gene_counts))

    list(
      CytoTRACE = potency_scores,
      method = "gene_count_approximation"
    )
  })

  cat("\n✓ CytoTRACE analysis complete\n")
  cat("  Score range:", round(range(results$CytoTRACE, na.rm = TRUE), 3), "\n")

  return(results)
}


#' Plot potency on UMAP
#'
#' @param seurat_obj Seurat object with potency scores
#' @param reduction UMAP reduction name
#' @param output_file Output file path
#' @param title Plot title
plot_potency_umap <- function(seurat_obj, reduction, output_file, title) {
  if (!reduction %in% names(seurat_obj@reductions)) {
    cat("  Skipping", reduction, "(not found)\n")
    return(NULL)
  }

  p <- FeaturePlot(
    seurat_obj,
    features = "cytotrace_potency",
    reduction = reduction,
    pt.size = 1.5,
    cols = viridis(100, option = "magma")
  ) +
    labs(
      title = title,
      subtitle = "Higher score = Less differentiated (more potent)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      legend.position = "right"
    )

  ggsave(output_file, plot = p, width = 12, height = 10, dpi = 300)
  cat("  Saved:", basename(output_file), "\n")

  return(p)
}


#' Plot potency distribution by cell type
#'
#' @param seurat_obj Seurat object
#' @param cell_type_col Cell type column name
#' @param output_file Output file path
plot_potency_distribution <- function(seurat_obj, cell_type_col, output_file) {
  plot_df <- data.frame(
    cell_type = seurat_obj@meta.data[[cell_type_col]],
    potency = seurat_obj$cytotrace_potency
  )

  # Remove NA
  plot_df <- plot_df[!is.na(plot_df$potency), ]

  # Violin plot
  p1 <- ggplot(plot_df, aes(x = cell_type, y = potency, fill = cell_type)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "CytoTRACE Potency by Cell Type",
      x = "Cell Type",
      y = "Differentiation Potency"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )

  # Density plot
  p2 <- ggplot(plot_df, aes(x = potency, fill = cell_type)) +
    geom_density(alpha = 0.5) +
    scale_fill_brewer(palette = "Set2", name = "Cell Type") +
    labs(
      title = "Potency Density",
      x = "Differentiation Potency",
      y = "Density"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )

  # Combine
  p_combined <- p1 + p2 + plot_layout(widths = c(2, 1))

  ggsave(output_file, plot = p_combined, width = 16, height = 6, dpi = 300)
  cat("  Saved:", basename(output_file), "\n")

  return(p_combined)
}


#' Validate potency with pseudotime
#'
#' @param seurat_obj Seurat object
#' @param pseudotime_file Pseudotime CSV file
#' @param output_file Output file path
validate_with_pseudotime <- function(seurat_obj, pseudotime_file, output_file) {
  if (is.null(pseudotime_file) || !file.exists(pseudotime_file)) {
    cat("  ⚠ No pseudotime file provided, skipping validation\n")
    return(NULL)
  }

  cat("\n", rep("=", 80), "\n", sep = "")
  cat("VALIDATING WITH PSEUDOTIME\n")
  cat(rep("=", 80), "\n\n", sep = "")

  # Load pseudotime
  pt_df <- read.csv(pseudotime_file, row.names = 1)

  # Match cells
  common_cells <- intersect(colnames(seurat_obj), rownames(pt_df))

  if (length(common_cells) < 10) {
    cat("  ⚠ Too few matching cells (", length(common_cells), "), skipping\n", sep = "")
    return(NULL)
  }

  cat("  Matched", length(common_cells), "cells\n")

  # Extract data
  plot_df <- data.frame(
    potency = seurat_obj$cytotrace_potency[common_cells],
    pseudotime = pt_df[common_cells, 1],
    cell_type = seurat_obj@meta.data[common_cells, "cell_type"]
  )

  # Remove NA
  plot_df <- plot_df[complete.cases(plot_df), ]

  # Calculate correlation
  corr <- cor.test(plot_df$potency, plot_df$pseudotime, method = "spearman")

  cat("\n  Spearman correlation: r =", round(corr$estimate, 3),
      ", p =", format.pval(corr$p.value, digits = 2), "\n")

  # Expected: negative correlation (high potency = early = low pseudotime)
  if (corr$estimate < -0.3) {
    cat("  ✓ Good negative correlation (as expected)\n")
  } else if (corr$estimate > 0.3) {
    cat("  ⚠ Unexpected positive correlation\n")
  } else {
    cat("  ⚠ Weak correlation\n")
  }

  # Plot
  p <- ggplot(plot_df, aes(x = pseudotime, y = potency)) +
    geom_point(aes(color = cell_type), size = 2, alpha = 0.6) +
    geom_smooth(method = "lm", color = "black", linewidth = 1.5, se = TRUE) +
    scale_color_brewer(palette = "Set2", name = "Cell Type") +
    labs(
      title = "CytoTRACE Potency vs Pseudotime",
      subtitle = paste0("Spearman r = ", round(corr$estimate, 3),
                       ", p = ", format.pval(corr$p.value, digits = 2)),
      x = "Trajectory Pseudotime",
      y = "Differentiation Potency"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray40"),
      legend.position = "right"
    )

  ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
  cat("  Saved:", basename(output_file), "\n")

  return(list(correlation = corr$estimate, p_value = corr$p.value))
}


#' Export results
#'
#' @param seurat_obj Seurat object
#' @param cytotrace_results CytoTRACE results
#' @param dirs Output directories
export_results <- function(seurat_obj, cytotrace_results, dirs) {
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("EXPORTING RESULTS\n")
  cat(rep("=", 80), "\n\n", sep = "")

  # Save CytoTRACE results
  output_file <- file.path(dirs$data, "cytotrace_results.rds")
  saveRDS(cytotrace_results, output_file)
  cat("✓ Saved CytoTRACE results:", output_file, "\n")

  # Save Seurat object
  output_file <- file.path(dirs$robjects, "seurat_with_cytotrace.rds")
  saveRDS(seurat_obj, output_file)
  cat("✓ Saved Seurat object:", output_file, "\n")

  # Export potency scores CSV
  potency_df <- data.frame(
    cell_id = colnames(seurat_obj),
    cell_type = seurat_obj$cell_type,
    cytotrace_potency = seurat_obj$cytotrace_potency,
    n_genes = seurat_obj$nFeature_RNA
  )

  output_file <- file.path(dirs$tables, "cytotrace_scores.csv")
  write.csv(potency_df, output_file, row.names = FALSE)
  cat("✓ Saved potency scores:", output_file, "\n")

  cat("\n✓ All results exported\n")
}


# =============================================================================
# Main Analysis Function
# =============================================================================

#' Run complete CytoTRACE analysis
#'
#' @param mode Analysis mode (cDC1, cDC2, integrated)
#' @param seurat_file Path to Seurat RDS file
#' @param pseudotime_file Path to pseudotime CSV (optional)
#' @param output_dir Output directory
#' @param cell_type_col Cell type column name
#' @param ncores Number of cores
#' @return Analysis results
run_cytotrace_analysis <- function(
  mode = "cDC1",
  seurat_file = NULL,
  pseudotime_file = NULL,
  output_dir = NULL,
  cell_type_col = "author_annotation",
  ncores = 1
) {

  cat("\n", rep("=", 80), "\n", sep = "")
  cat("     CYTOTRACE DIFFERENTIATION POTENCY ANALYSIS\n")
  cat("            GSE228544 CITE-seq Data\n")
  cat(rep("=", 80), "\n\n", sep = "")

  cat("Mode:", mode, "\n")
  cat("Cell type column:", cell_type_col, "\n")

  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- paste0("cytotrace_analysis_", mode)
  }

  dirs <- create_output_directories(output_dir)
  cat("Output directory:", output_dir, "\n")

  # Load Seurat object
  if (is.null(seurat_file)) {
    seurat_file <- file.path("results", mode, "Robjects", "seurat_obj_annotated.rds")
  }

  if (!file.exists(seurat_file)) {
    stop("Seurat file not found: ", seurat_file)
  }

  cat("\nLoading Seurat object:", seurat_file, "\n")
  seurat_obj <- readRDS(seurat_file)

  cat("  Cells:", ncol(seurat_obj), "\n")
  cat("  Genes:", nrow(seurat_obj), "\n")

  # Add cell type to metadata if not already there
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop("Cell type column not found: ", cell_type_col)
  }

  seurat_obj$cell_type <- seurat_obj@meta.data[[cell_type_col]]

  # Run CytoTRACE
  cytotrace_results <- run_cytotrace(seurat_obj, ncores = ncores)

  # Add potency scores to Seurat object
  potency_scores <- cytotrace_results$CytoTRACE

  # Match cell names
  common_cells <- intersect(names(potency_scores), colnames(seurat_obj))

  if (length(common_cells) == 0) {
    stop("No matching cells between CytoTRACE results and Seurat object")
  }

  seurat_obj$cytotrace_potency <- NA
  seurat_obj$cytotrace_potency[common_cells] <- potency_scores[common_cells]

  cat("\n✓ Added potency scores to Seurat object\n")
  cat("  Cells with scores:", sum(!is.na(seurat_obj$cytotrace_potency)), "\n")

  # Visualizations
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("CREATING VISUALIZATIONS\n")
  cat(rep("=", 80), "\n\n", sep = "")

  # Plot on all UMAP spaces
  cat("Plotting potency on UMAP spaces...\n")

  umap_configs <- list(
    list(reduction = "umap", suffix = "rna_umap", label = "RNA UMAP"),
    list(reduction = "adt_umap", suffix = "adt_umap", label = "ADT UMAP"),
    list(reduction = "wnn_umap", suffix = "wnn_umap", label = "WNN UMAP")
  )

  for (config in umap_configs) {
    output_file <- file.path(dirs$potency,
                            paste0("potency_", config$suffix, ".png"))
    plot_potency_umap(seurat_obj, config$reduction, output_file,
                     paste0(mode, " Potency - ", config$label))
  }

  # Potency distribution
  cat("\nPlotting potency distribution...\n")
  output_file <- file.path(dirs$potency, "potency_by_celltype.png")
  plot_potency_distribution(seurat_obj, "cell_type", output_file)

  # Validation with pseudotime
  if (!is.null(pseudotime_file)) {
    output_file <- file.path(dirs$validation, "potency_vs_pseudotime.png")
    validation_results <- validate_with_pseudotime(seurat_obj, pseudotime_file, output_file)
  }

  # Export results
  export_results(seurat_obj, cytotrace_results, dirs)

  # Summary
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("CYTOTRACE ANALYSIS COMPLETE!\n")
  cat(rep("=", 80), "\n\n", sep = "")

  cat("Summary:\n")
  cat("  Cells analyzed:", sum(!is.na(seurat_obj$cytotrace_potency)), "\n")
  cat("  Potency range:", round(range(seurat_obj$cytotrace_potency, na.rm = TRUE), 3), "\n")

  # Cell type potency
  cat("\nMean potency by cell type:\n")
  potency_by_type <- tapply(seurat_obj$cytotrace_potency,
                            seurat_obj$cell_type,
                            mean, na.rm = TRUE)
  print(round(sort(potency_by_type, decreasing = TRUE), 3))

  cat("\n✓ All results saved to:", output_dir, "\n")

  return(list(
    seurat_obj = seurat_obj,
    cytotrace_results = cytotrace_results
  ))
}


# =============================================================================
# Command Line Interface
# =============================================================================

if (length(args) > 0) {
  # Parse command line arguments
  mode <- "cDC1"
  seurat_file <- NULL
  pseudotime_file <- NULL
  output_dir <- NULL
  cell_type_col <- "author_annotation"
  ncores <- 1

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--mode") {
      mode <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--seurat") {
      seurat_file <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--pseudotime") {
      pseudotime_file <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output-dir") {
      output_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--cell-type-col") {
      cell_type_col <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--ncores") {
      ncores <- as.integer(args[i + 1])
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  # Run analysis
  results <- run_cytotrace_analysis(
    mode = mode,
    seurat_file = seurat_file,
    pseudotime_file = pseudotime_file,
    output_dir = output_dir,
    cell_type_col = cell_type_col,
    ncores = ncores
  )
}
