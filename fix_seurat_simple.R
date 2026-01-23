#!/usr/bin/env Rscript
# 简化版 Seurat 对象修复脚本（避免 Seurat v5 复杂性）

library(Seurat)

cat("=== 简化版 Seurat 对象修复 ===\n\n")

# 加载原始对象
input_file <- "results/cDC1/Robjects/seurat_obj_annotated.rds"
cat("加载:", input_file, "\n")
seurat_obj <- readRDS(input_file)

cat("\n原始对象:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  Assays:", paste(Assays(seurat_obj), collapse=", "), "\n")

# 方案：创建一个全新的简化 Seurat 对象（兼容性最好）
cat("\n创建简化的 Seurat 对象...\n")

# 1. 提取 RNA 数据
cat("  提取 RNA 数据...\n")
DefaultAssay(seurat_obj) <- "RNA"
rna_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

# 2. 提取 ADT 数据
cat("  提取 ADT 数据...\n")
DefaultAssay(seurat_obj) <- "ADT"
adt_counts <- GetAssayData(seurat_obj, assay = "ADT", slot = "counts")

cat("  RNA 特征:", nrow(rna_counts), "\n")
cat("  ADT 特征:", nrow(adt_counts), "\n")

# 3. 创建新的 Seurat 对象（只包含必要成分）
cat("\n创建新对象...\n")
seurat_new <- CreateSeuratObject(
  counts = rna_counts,
  meta.data = seurat_obj@meta.data
)

# 4. 添加 ADT assay
cat("  添加 ADT assay...\n")
seurat_new[["ADT"]] <- CreateAssayObject(counts = adt_counts)

# 5. 归一化
cat("  归一化 RNA...\n")
DefaultAssay(seurat_new) <- "RNA"
seurat_new <- NormalizeData(seurat_new, verbose = FALSE)

cat("  归一化 ADT (CLR)...\n")
DefaultAssay(seurat_new) <- "ADT"
seurat_new <- NormalizeData(seurat_new,
                             normalization.method = 'CLR',
                             margin = 2,
                             verbose = FALSE)

# 6. 复制降维结果
cat("  复制降维结果...\n")
for (reduction_name in names(seurat_obj@reductions)) {
  seurat_new@reductions[[reduction_name]] <- seurat_obj@reductions[[reduction_name]]
  cat("    ✓", reduction_name, "\n")
}

# 7. 创建兼容的降维名称
if ("wnn.umap" %in% names(seurat_new@reductions)) {
  seurat_new@reductions$wnn_umap <- seurat_new@reductions$wnn.umap
  cat("    ✓ wnn_umap (别名)\n")
}

if ("umap" %in% names(seurat_new@reductions)) {
  seurat_new@reductions$rna_umap <- seurat_new@reductions$umap
  cat("    ✓ rna_umap (别名)\n")
}

# 8. 保存
output_file <- "seurat_obj_citeviz_ready.rds"
cat("\n保存:", output_file, "\n")
saveRDS(seurat_new, output_file)

# 验证
cat("\n=== 新对象信息 ===\n")
cat("细胞数:", ncol(seurat_new), "\n")
cat("Assays:", paste(Assays(seurat_new), collapse=", "), "\n")
cat("降维:", paste(names(seurat_new@reductions), collapse=", "), "\n")
cat("ADT 数量:", nrow(seurat_new@assays$ADT), "\n")

# 测试数据访问
cat("\n测试数据访问...\n")
tryCatch({
  DefaultAssay(seurat_new) <- "ADT"
  test_data <- GetAssayData(seurat_new, assay = "ADT", slot = "data")
  cat("  ✓ ADT 数据可访问\n")
  cat("  ✓ 维度:", dim(test_data), "\n")
}, error = function(e) {
  cat("  ✗ 错误:", e$message, "\n")
})

cat("\n✓ 完成！\n")
cat("新文件:", output_file, "\n")
cat("此文件应该可以直接在 CITEViz 中使用\n")
