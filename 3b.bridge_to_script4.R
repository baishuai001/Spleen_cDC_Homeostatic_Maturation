#!/usr/bin/env Rscript

################################################################################
# 衔接脚本：将 CITE-seq_Complete_Analysis_v3.Rmd 输出转换为脚本4兼容格式
################################################################################
#
# 功能说明：
# 1. 加载 v3.Rmd 生成的 pre_wnn.rds 对象（包含 ADT 和 RNA 独立分析）
# 2. 重命名降维结果以匹配脚本4的命名规范
# 3. 调整列名以符合脚本4的期望
# 4. 保存为脚本4可直接使用的格式
#
# 使用方法：
# Rscript 3b.bridge_to_script4.R [analysis_mode]
# 例如：Rscript 3b.bridge_to_script4.R cDC1
#
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
})

# ==================== 参数设置 ====================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  analysis_mode <- args[1]
} else {
  analysis_mode <- "cDC1"  # 默认模式
}

# 验证模式
valid_modes <- c("cDC1", "cDC2", "integrated", "full")
if (!analysis_mode %in% valid_modes) {
  stop("❌ 无效的分析模式: ", analysis_mode,
       "\n   可选: ", paste(valid_modes, collapse = ", "))
}

cat("========================================\n")
cat("  衔接脚本：v3.Rmd → 脚本4\n")
cat("  分析模式:", analysis_mode, "\n")
cat("========================================\n\n")

# ==================== 1. 加载对象 ====================
output_dir <- "results/"
project_name <- paste0("GSE228544_", analysis_mode)

# 优先使用 pre_wnn.rds（更"干净"，只包含 ADT 和 RNA）
input_file <- paste0(output_dir, "Robjects/", project_name, "_pre_wnn.rds")

if (!file.exists(input_file)) {
  cat("⚠️  Pre-WNN 对象不存在，尝试使用最终对象...\n")
  input_file <- paste0(output_dir, "Robjects/", project_name, "_final.rds")
}

if (!file.exists(input_file)) {
  stop("❌ 找不到输入文件: ", input_file,
       "\n   请先运行 CITE-seq_Complete_Analysis_v3.Rmd")
}

cat("📂 加载对象:", input_file, "\n")
seuratObj <- readRDS(input_file)

cat("✓ 对象加载完成\n")
cat("  细胞数:", ncol(seuratObj), "\n")
cat("  基因数:", nrow(seuratObj), "\n")
cat("  ADT 数:", nrow(seuratObj[["ADT"]]), "\n\n")

# ==================== 2. 检查现有降维 ====================
cat("📊 检查现有降维结果...\n")
reductions_available <- names(seuratObj@reductions)
cat("  现有降维:", paste(reductions_available, collapse = ", "), "\n\n")

# ==================== 3. 重命名降维结果（匹配脚本4命名） ====================
cat("🔄 重命名降维结果以匹配脚本4期望...\n\n")

# --- RNA UMAP: umap → SCT_umap ---
if ("umap" %in% reductions_available) {
  cat("  [1/4] RNA UMAP: umap → SCT_umap\n")

  # 复制降维对象
  seuratObj[["SCT_umap"]] <- seuratObj[["umap"]]

  # 修正列名
  old_colnames <- colnames(seuratObj[["SCT_umap"]]@cell.embeddings)
  cat("    当前列名:", paste(old_colnames, collapse = ", "), "\n")

  colnames(seuratObj[["SCT_umap"]]@cell.embeddings) <- c("sctUMAP_1", "sctUMAP_2")
  cat("    新列名: sctUMAP_1, sctUMAP_2\n")
  cat("    ✓ 完成\n\n")
} else {
  cat("  [1/4] ⚠️  未找到 RNA UMAP (umap)\n\n")
}

# --- RNA tSNE: tsne → SCT_tsne ---
if ("tsne" %in% reductions_available) {
  cat("  [2/4] RNA tSNE: tsne → SCT_tsne\n")

  # 复制降维对象
  seuratObj[["SCT_tsne"]] <- seuratObj[["tsne"]]

  # 修正列名
  old_colnames <- colnames(seuratObj[["SCT_tsne"]]@cell.embeddings)
  cat("    当前列名:", paste(old_colnames, collapse = ", "), "\n")

  colnames(seuratObj[["SCT_tsne"]]@cell.embeddings) <- c("sctTSNE_1", "sctTSNE_2")
  cat("    新列名: sctTSNE_1, sctTSNE_2\n")
  cat("    ✓ 完成\n\n")
} else {
  cat("  [2/4] ⚠️  未找到 RNA tSNE (tsne)\n\n")
}

# --- ADT UMAP: adt_umap → ADT_umap ---
if ("adt_umap" %in% reductions_available) {
  cat("  [3/4] ADT UMAP: adt_umap → ADT_umap\n")

  # 复制降维对象
  seuratObj[["ADT_umap"]] <- seuratObj[["adt_umap"]]

  # 修正列名
  old_colnames <- colnames(seuratObj[["ADT_umap"]]@cell.embeddings)
  cat("    当前列名:", paste(old_colnames, collapse = ", "), "\n")

  colnames(seuratObj[["ADT_umap"]]@cell.embeddings) <- c("adtUMAP_1", "adtUMAP_2")
  cat("    新列名: adtUMAP_1, adtUMAP_2\n")
  cat("    ✓ 完成\n\n")
} else {
  cat("  [3/4] ⚠️  未找到 ADT UMAP (adt_umap)\n\n")
}

# --- ADT tSNE: adt_tsne → ADT_tsne ---
if ("adt_tsne" %in% reductions_available) {
  cat("  [4/4] ADT tSNE: adt_tsne → ADT_tsne\n")

  # 复制降维对象
  seuratObj[["ADT_tsne"]] <- seuratObj[["adt_tsne"]]

  # 修正列名
  old_colnames <- colnames(seuratObj[["ADT_tsne"]]@cell.embeddings)
  cat("    当前列名:", paste(old_colnames, collapse = ", "), "\n")

  colnames(seuratObj[["ADT_tsne"]]@cell.embeddings) <- c("adtTSNE_1", "adtTSNE_2")
  cat("    新列名: adtTSNE_1, adtTSNE_2\n")
  cat("    ✓ 完成\n\n")
} else {
  cat("  [4/4] ⚠️  未找到 ADT tSNE (adt_tsne)\n\n")
}

# ==================== 4. 验证所有必需组件 ====================
cat("========================================\n")
cat("✅ 验证脚本4所需组件\n")
cat("========================================\n\n")

# 检查降维
cat("【降维结果】\n")
required_reductions <- c("SCT_umap", "SCT_tsne", "ADT_umap", "ADT_tsne")
all_reductions_ok <- TRUE

for (red in required_reductions) {
  if (red %in% names(seuratObj@reductions)) {
    # 检查列名
    colnames_current <- colnames(seuratObj[[red]]@cell.embeddings)
    cat("  ✓", red, ":", paste(colnames_current, collapse = ", "), "\n")
  } else {
    cat("  ✗", red, "(缺失)\n")
    all_reductions_ok <- FALSE
  }
}

# 检查 metadata
cat("\n【Meta.data 列】\n")
required_metadata <- c("SCT_clusters", "ADT_clusters", "nCount_RNA", "subsets_Mito_percent")
all_metadata_ok <- TRUE

for (meta in required_metadata) {
  if (meta %in% colnames(seuratObj@meta.data)) {
    if (grepl("clusters", meta)) {
      n_clusters <- length(unique(seuratObj@meta.data[[meta]]))
      cat("  ✓", meta, "(", n_clusters, "clusters )\n")
    } else if (grepl("percent|nCount", meta)) {
      range_val <- range(seuratObj@meta.data[[meta]], na.rm = TRUE)
      cat("  ✓", meta, "( 范围:", round(range_val[1], 2), "-", round(range_val[2], 2), ")\n")
    } else {
      cat("  ✓", meta, "\n")
    }
  } else {
    cat("  ✗", meta, "(缺失)\n")
    all_metadata_ok <- FALSE
  }
}

# 检查 Assays
cat("\n【Assays】\n")
required_assays <- c("RNA", "SCT", "ADT")
all_assays_ok <- TRUE

for (assay in required_assays) {
  if (assay %in% names(seuratObj@assays)) {
    n_features <- nrow(seuratObj[[assay]])
    cat("  ✓", assay, "(", n_features, "features )\n")
  } else {
    cat("  ✗", assay, "(缺失)\n")
    all_assays_ok <- FALSE
  }
}

cat("\n")

# ==================== 5. 创建输出目录结构 ====================
cat("📁 创建脚本4兼容的目录结构...\n")

output_folder <- paste0("SAM2and3_", analysis_mode, "/")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, "results/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, "results/Robjects/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, "results/QC/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, "results/Annotation/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_folder, "results/Marker_lists/"), showWarnings = FALSE, recursive = TRUE)

cat("  ✓ 目录结构已创建:", output_folder, "\n\n")

# ==================== 6. 保存兼容对象 ====================
cat("💾 保存脚本4兼容对象...\n")

output_file <- paste0(output_folder, "results/Robjects/seuratObj_", analysis_mode, "_for_script4.rds")
saveRDS(seuratObj, output_file)

file_size_mb <- round(file.size(output_file) / 1024 / 1024, 2)
cat("  ✓ 对象已保存:", output_file, "\n")
cat("  文件大小:", file_size_mb, "MB\n\n")

# ==================== 7. 生成摘要报告 ====================
cat("========================================\n")
cat("📊 衔接完成摘要\n")
cat("========================================\n\n")

cat("【输入】\n")
cat("  文件:", basename(input_file), "\n\n")

cat("【输出】\n")
cat("  目录:", output_folder, "\n")
cat("  文件:", basename(output_file), "\n\n")

cat("【数据概览】\n")
cat("  分析模式:", analysis_mode, "\n")
cat("  细胞数:", ncol(seuratObj), "\n")
cat("  RNA 特征:", nrow(seuratObj[["RNA"]]), "\n")
cat("  ADT 特征:", nrow(seuratObj[["ADT"]]), "\n")
cat("  SCT clusters:", length(unique(seuratObj$SCT_clusters)), "\n")
cat("  ADT clusters:", length(unique(seuratObj$ADT_clusters)), "\n\n")

cat("【降维结果】\n")
cat("  ✓ SCT_umap ( sctUMAP_1, sctUMAP_2 )\n")
cat("  ✓ SCT_tsne ( sctTSNE_1, sctTSNE_2 )\n")
cat("  ✓ ADT_umap ( adtUMAP_1, adtUMAP_2 )\n")
cat("  ✓ ADT_tsne ( adtTSNE_1, adtTSNE_2 )\n\n")

# 最终状态检查
if (all_reductions_ok && all_metadata_ok && all_assays_ok) {
  cat("========================================\n")
  cat("  ✅ 衔接成功！\n")
  cat("  现在可以运行脚本4了\n")
  cat("========================================\n\n")

  cat("【下一步】\n")
  cat("1. 修改脚本4的开头：\n")
  cat("   analysis_mode <- \"", analysis_mode, "\"\n", sep = "")
  cat("   sampleName <- \"SAM2and3_", analysis_mode, "\"\n", sep = "")
  cat("   seuratObj <- readRDS(\"", output_folder, "results/Robjects/seuratObj_",
      analysis_mode, "_for_script4.rds\")\n\n", sep = "")
  cat("2. 运行脚本4：\n")
  cat("   Rscript 4.script_CITEseq_SAM_WT_aggr.R\n\n")

} else {
  cat("========================================\n")
  cat("  ⚠️  衔接完成，但存在缺失组件\n")
  cat("  请检查上述验证结果\n")
  cat("========================================\n\n")
}
