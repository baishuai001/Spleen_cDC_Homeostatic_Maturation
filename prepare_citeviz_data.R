#!/usr/bin/env Rscript
# 准备 CITEViz 数据文件（下采样版本）

library(Seurat)

cat("=================================================\n")
cat("准备 CITEViz 数据文件\n")
cat("=================================================\n\n")

# 参数设置
ORIGINAL_FILE <- "results/cDC1/Robjects/seurat_obj_annotated.rds"
COMPATIBLE_FILE <- "results/integrated/Robjects/seurat_obj_citeviz_compatible.rds"  # 可能已存在的兼容文件
OUTPUT_FAST <- "citeviz_fast_5k.rds"      # 快速版（5000细胞）
OUTPUT_MEDIUM <- "citeviz_medium_8k.rds"  # 中等版（8000细胞）
OUTPUT_FULL <- "citeviz_full_compatible.rds"  # 完整版（所有细胞）

# 智能加载：优先使用已有的兼容文件
if (file.exists(COMPATIBLE_FILE)) {
  cat("✓ 发现已有兼容文件:", COMPATIBLE_FILE, "\n")
  cat("  跳过兼容性处理，直接使用\n\n")
  seurat_obj <- readRDS(COMPATIBLE_FILE)
} else if (file.exists(OUTPUT_FULL)) {
  cat("✓ 发现本地完整兼容文件:", OUTPUT_FULL, "\n")
  cat("  跳过兼容性处理，直接使用\n\n")
  seurat_obj <- readRDS(OUTPUT_FULL)
} else {
  cat("未找到兼容文件，从原始数据创建...\n")
  cat("  加载:", ORIGINAL_FILE, "\n")
  seurat_obj <- readRDS(ORIGINAL_FILE)

  # 稍后会创建兼容对象
  NEEDS_COMPATIBILITY <- TRUE
}

cat("\n数据信息:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数:", nrow(seurat_obj@assays$RNA), "\n")
cat("  ADT数:", nrow(seurat_obj@assays$ADT), "\n\n")

# 创建兼容对象的函数
create_compatible <- function(seurat_obj) {
  cat("创建兼容对象...\n")

  # 提取数据
  DefaultAssay(seurat_obj) <- "RNA"
  rna_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

  DefaultAssay(seurat_obj) <- "ADT"
  adt_counts <- GetAssayData(seurat_obj, assay = "ADT", slot = "counts")

  # 创建新对象
  seurat_new <- CreateSeuratObject(
    counts = rna_counts,
    meta.data = seurat_obj@meta.data
  )

  seurat_new[["ADT"]] <- CreateAssayObject(counts = adt_counts)

  # 归一化
  DefaultAssay(seurat_new) <- "RNA"
  seurat_new <- NormalizeData(seurat_new, verbose = FALSE)

  DefaultAssay(seurat_new) <- "ADT"
  seurat_new <- NormalizeData(seurat_new,
                               normalization.method = 'CLR',
                               margin = 2,
                               verbose = FALSE)

  # 复制降维
  for (reduction_name in names(seurat_obj@reductions)) {
    seurat_new@reductions[[reduction_name]] <- seurat_obj@reductions[[reduction_name]]
  }

  # 创建兼容的降维名称
  if ("wnn.umap" %in% names(seurat_new@reductions)) {
    seurat_new@reductions$wnn_umap <- seurat_new@reductions$wnn.umap
  }
  if ("umap" %in% names(seurat_new@reductions)) {
    seurat_new@reductions$rna_umap <- seurat_new@reductions$umap
  }

  return(seurat_new)
}

# 智能采样函数
smart_sample <- function(seurat_obj, n_cells) {
  if (ncol(seurat_obj) <= n_cells) {
    return(seurat_obj)
  }

  cat("  采样到", n_cells, "个细胞...\n")

  if ("author_annotation" %in% colnames(seurat_obj@meta.data)) {
    # 按细胞类型比例采样
    cell_types <- seurat_obj$author_annotation
    cells_per_type <- table(cell_types)

    proportions <- cells_per_type / sum(cells_per_type)
    sample_per_type <- round(proportions * n_cells)
    sample_per_type[sample_per_type < 30] <- min(30, cells_per_type[sample_per_type < 30])

    set.seed(42)
    sampled_cells <- c()
    for (ct in names(sample_per_type)) {
      type_cells <- colnames(seurat_obj)[cell_types == ct]
      n_sample <- min(sample_per_type[ct], length(type_cells))
      sampled_cells <- c(sampled_cells, sample(type_cells, n_sample))
    }

    seurat_obj <- seurat_obj[, sampled_cells]
  } else {
    # 随机采样
    set.seed(42)
    sampled_cells <- sample(colnames(seurat_obj), n_cells)
    seurat_obj <- seurat_obj[, sampled_cells]
  }

  cat("  采样后细胞数:", ncol(seurat_obj), "\n")
  return(seurat_obj)
}

# 如果需要，先创建完整的兼容对象
if (exists("NEEDS_COMPATIBILITY") && NEEDS_COMPATIBILITY) {
  cat("\n创建完整兼容对象（首次运行）...\n")
  seurat_obj <- create_compatible(seurat_obj)
  # 保存完整版
  saveRDS(seurat_obj, OUTPUT_FULL)
  cat("✓ 保存完整兼容对象:", OUTPUT_FULL, "\n")
}

# 生成文件 1：快速版（推荐）
cat("\n生成文件 1/3: 快速版（5000 细胞）\n")
cat("用途：快速探索，响应最快，推荐首次使用\n")
seurat_fast <- smart_sample(seurat_obj, 5000)
# 不需要再次 create_compatible，因为 seurat_obj 已经是兼容的
saveRDS(seurat_fast, OUTPUT_FAST)
cat("✓ 保存:", OUTPUT_FAST, "\n")
cat("  文件大小:", round(file.size(OUTPUT_FAST) / 1024^2, 1), "MB\n")

# 生成文件 2：中等版
cat("\n生成文件 2/3: 中等版（8000 细胞）\n")
cat("用途：平衡性能和细节，适合常规分析\n")
seurat_medium <- smart_sample(seurat_obj, 8000)
saveRDS(seurat_medium, OUTPUT_MEDIUM)
cat("✓ 保存:", OUTPUT_MEDIUM, "\n")
cat("  文件大小:", round(file.size(OUTPUT_MEDIUM) / 1024^2, 1), "MB\n")

# 生成文件 3：完整版（如果还没有）
if (!exists("NEEDS_COMPATIBILITY") || !NEEDS_COMPATIBILITY) {
  cat("\n生成文件 3/3: 完整版（复制现有兼容对象）\n")
  cat("用途：最终分析，包含所有数据（可能较慢）\n")
  # 只是复制，不需要重新创建
  file.copy(
    if (file.exists(COMPATIBLE_FILE)) COMPATIBLE_FILE else OUTPUT_FULL,
    OUTPUT_FULL,
    overwrite = TRUE
  )
  cat("✓ 保存:", OUTPUT_FULL, "\n")
  cat("  文件大小:", round(file.size(OUTPUT_FULL) / 1024^2, 1), "MB\n")
}

# 总结
cat("\n=================================================\n")
cat("数据准备完成！\n")
cat("=================================================\n\n")

cat("生成的文件:\n")
cat("1. ", OUTPUT_FAST, " (快速版，推荐首次使用)\n")
cat("2. ", OUTPUT_MEDIUM, " (中等版，平衡性能)\n")
cat("3. ", OUTPUT_FULL, " (完整版，最详细)\n\n")

cat("使用方法:\n")
cat("1. 在 R 中启动 CITEViz:\n")
cat("   library(CITEViz)\n")
cat("   run_app()  # 或 launchApp()\n\n")

cat("2. 在 CITEViz 网页界面点击 'Browse...'\n")
cat("3. 上传其中一个文件（推荐从快速版开始）\n\n")

cat("提示:\n")
cat("- 快速版: 响应最快，适合快速浏览和门控\n")
cat("- 中等版: 兼顾性能和细节\n")
cat("- 完整版: 包含所有细胞，但可能卡顿\n")
