#!/usr/bin/env Rscript
# 修复 Seurat 对象以兼容 CITEViz

library(Seurat)

cat("=== 诊断和修复 Seurat 对象 ===\n\n")

# 加载原始对象
input_file <- "results/cDC1/Robjects/seurat_obj_annotated.rds"
cat("加载:", input_file, "\n")
seurat_obj <- readRDS(input_file)

cat("\n原始对象信息:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  Assays:", paste(Assays(seurat_obj), collapse=", "), "\n")
cat("  降维:", paste(names(seurat_obj@reductions), collapse=", "), "\n")

# 问题 1: 检查 ADT assay
cat("\n检查 ADT assay...\n")
if (!"ADT" %in% Assays(seurat_obj)) {
  stop("错误：未找到 ADT assay")
}

# 检查 ADT 数据层
DefaultAssay(seurat_obj) <- "ADT"
adt_assay <- seurat_obj@assays$ADT

cat("  ADT 特征数:", nrow(adt_assay), "\n")
cat("  ADT 名称示例:\n")
print(head(rownames(adt_assay), 10))

# 问题 2: 确保 ADT 有正确的归一化数据
cat("\n检查 ADT 归一化...\n")

# 对于 Seurat v5，可能需要 JoinLayers
if (packageVersion("Seurat") >= "5.0.0") {
  cat("  检测到 Seurat v5，合并数据层...\n")
  seurat_obj <- JoinLayers(seurat_obj, assay = "ADT")
}

# 确保有 data slot
if (is.null(slot(adt_assay, "data")) || length(slot(adt_assay, "data")) == 0) {
  cat("  执行 CLR 归一化...\n")
  seurat_obj <- NormalizeData(seurat_obj,
                               normalization.method = 'CLR',
                               margin = 2,
                               assay = "ADT")
} else {
  cat("  ✓ ADT 数据已归一化\n")
}

# 问题 3: 确保降维名称兼容
cat("\n检查降维名称...\n")
reductions <- names(seurat_obj@reductions)
cat("  现有降维:", paste(reductions, collapse=", "), "\n")

# CITEViz 可能期望特定的名称格式
# 如果有 wnn.umap，也创建一个 wnn_umap 的副本
if ("wnn.umap" %in% reductions && !"wnn_umap" %in% reductions) {
  cat("  创建 wnn_umap 副本（兼容性）...\n")
  seurat_obj@reductions$wnn_umap <- seurat_obj@reductions$wnn.umap
}

if ("umap" %in% reductions && !"rna_umap" %in% reductions) {
  cat("  创建 rna_umap 副本（兼容性）...\n")
  seurat_obj@reductions$rna_umap <- seurat_obj@reductions$umap
}

# 问题 4: 简化对象（移除不必要的数据）
cat("\n简化对象...\n")

# 只保留必要的 assays
keep_assays <- c("RNA", "ADT")
for (assay_name in Assays(seurat_obj)) {
  if (!assay_name %in% keep_assays) {
    cat("  移除 assay:", assay_name, "\n")
    seurat_obj[[assay_name]] <- NULL
  }
}

# 问题 5: 确保元数据简洁
cat("\n检查元数据...\n")
meta_cols <- colnames(seurat_obj@meta.data)
cat("  元数据列数:", length(meta_cols), "\n")

# 保存修复后的对象
output_file <- "seurat_obj_citeviz_compatible.rds"
cat("\n保存兼容对象:", output_file, "\n")
saveRDS(seurat_obj, output_file)

# 最终验证
cat("\n=== 修复后的对象 ===\n")
cat("细胞数:", ncol(seurat_obj), "\n")
cat("Assays:", paste(Assays(seurat_obj), collapse=", "), "\n")
cat("降维:", paste(names(seurat_obj@reductions), collapse=", "), "\n")
cat("ADT 数量:", nrow(seurat_obj@assays$ADT), "\n")

cat("\n✓ 修复完成！\n")
cat("现在可以用这个文件启动 CITEViz:\n")
cat("  Rscript run_citeviz.R（需要修改脚本中的文件路径）\n")
cat("  或在 CITEViz 网页上上传:", output_file, "\n")
