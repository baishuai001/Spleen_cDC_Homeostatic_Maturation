#!/usr/bin/env Rscript
# CITEViz 快速启动脚本 - 下采样版本（解决卡顿问题）

library(Seurat)
library(CITEViz)

cat("=================================================\n")
cat("CITEViz 快速启动脚本（下采样版）\n")
cat("=================================================\n\n")

# 参数设置
DATA_PATH <- "results/cDC1/Robjects/seurat_obj_annotated.rds"
SAMPLE_SIZE <- 5000  # 采样细胞数（可调整：3000-8000）
SEED <- 42

cat("加载 Seurat 对象...\n")
seurat_obj <- readRDS(DATA_PATH)

cat("原始数据:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数:", nrow(seurat_obj), "\n")

# 智能采样：保留每种细胞类型的代表性
cat("\n进行智能下采样...\n")
cat("  目标细胞数:", SAMPLE_SIZE, "\n")

if ("author_annotation" %in% colnames(seurat_obj@meta.data)) {
  # 按细胞类型分层采样
  cell_types <- seurat_obj$author_annotation
  cells_per_type <- table(cell_types)

  cat("  采样策略: 按细胞类型比例分层采样\n")

  # 计算每种细胞类型应采样的数量
  proportions <- cells_per_type / sum(cells_per_type)
  sample_per_type <- round(proportions * SAMPLE_SIZE)

  # 确保至少每种类型有一些细胞
  sample_per_type[sample_per_type < 50] <- min(50, cells_per_type[sample_per_type < 50])

  # 采样
  set.seed(SEED)
  sampled_cells <- c()
  for (ct in names(sample_per_type)) {
    type_cells <- colnames(seurat_obj)[cell_types == ct]
    n_sample <- min(sample_per_type[ct], length(type_cells))
    sampled_cells <- c(sampled_cells, sample(type_cells, n_sample))
  }

  seurat_obj_subset <- seurat_obj[, sampled_cells]

} else {
  # 简单随机采样
  cat("  采样策略: 随机采样\n")
  set.seed(SEED)

  if (ncol(seurat_obj) > SAMPLE_SIZE) {
    sampled_cells <- sample(colnames(seurat_obj), SAMPLE_SIZE)
    seurat_obj_subset <- seurat_obj[, sampled_cells]
  } else {
    seurat_obj_subset <- seurat_obj
  }
}

cat("\n采样后数据:\n")
cat("  细胞数:", ncol(seurat_obj_subset), "\n")
cat("  减少:", ncol(seurat_obj) - ncol(seurat_obj_subset), "个细胞",
    "(", round((1 - ncol(seurat_obj_subset)/ncol(seurat_obj))*100, 1), "%)\n")

# 显示细胞类型分布
if ("author_annotation" %in% colnames(seurat_obj_subset@meta.data)) {
  cat("\n采样后细胞类型分布:\n")
  print(table(seurat_obj_subset$author_annotation))
}

# 可选：进一步简化 ADT 数据（只保留高变 ADT）
if ("ADT" %in% Assays(seurat_obj_subset)) {
  n_adt <- nrow(seurat_obj_subset@assays$ADT)
  if (n_adt > 50) {
    cat("\n✓ ADT 数据已包含，共", n_adt, "个标记\n")
    cat("  提示：如果仍然卡顿，可以只保留关键 ADT 标记\n")
  }
}

# 启动 CITEViz
cat("\n启动 CITEViz（优化版）...\n")
cat("提示：\n")
cat("  - 采样后的数据应该更流畅\n")
cat("  - 如果仍然卡顿，减小 SAMPLE_SIZE\n")
cat("  - 建议：3000-5000 细胞最流畅\n\n")

# 启动应用
if (exists("run_app", where = "package:CITEViz")) {
  run_app(seurat_object = seurat_obj_subset)
} else if (exists("launchApp", where = "package:CITEViz")) {
  launchApp(seurat_obj_subset)
} else {
  cat("可用函数:\n")
  print(ls("package:CITEViz"))
  stop("未找到启动函数")
}
