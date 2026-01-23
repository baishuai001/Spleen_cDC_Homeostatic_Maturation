#!/usr/bin/env Rscript
# CITEViz 一键启动脚本（自动处理兼容性）

library(Seurat)
library(CITEViz)

cat("=================================================\n")
cat("CITEViz 一键启动\n")
cat("=================================================\n\n")

# 参数设置
ORIGINAL_FILE <- "results/cDC1/Robjects/seurat_obj_annotated.rds"
COMPATIBLE_FILE <- "seurat_obj_citeviz_ready.rds"
USE_FAST_MODE <- TRUE  # 使用快速模式（下采样）
SAMPLE_SIZE <- 5000    # 快速模式的细胞数

# 步骤 1：检查是否已有兼容文件
if (!file.exists(COMPATIBLE_FILE)) {
  cat("步骤 1/3: 创建兼容的 Seurat 对象...\n")
  cat("（首次运行需要几分钟，之后会复用此文件）\n\n")

  # 加载原始对象
  cat("  加载原始数据:", ORIGINAL_FILE, "\n")
  seurat_obj <- readRDS(ORIGINAL_FILE)

  # 创建简化对象
  cat("  提取 RNA 数据...\n")
  DefaultAssay(seurat_obj) <- "RNA"
  rna_counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

  cat("  提取 ADT 数据...\n")
  DefaultAssay(seurat_obj) <- "ADT"
  adt_counts <- GetAssayData(seurat_obj, assay = "ADT", slot = "counts")

  cat("  创建新 Seurat 对象...\n")
  seurat_new <- CreateSeuratObject(
    counts = rna_counts,
    meta.data = seurat_obj@meta.data
  )

  cat("  添加 ADT assay...\n")
  seurat_new[["ADT"]] <- CreateAssayObject(counts = adt_counts)

  cat("  归一化数据...\n")
  DefaultAssay(seurat_new) <- "RNA"
  seurat_new <- NormalizeData(seurat_new, verbose = FALSE)

  DefaultAssay(seurat_new) <- "ADT"
  seurat_new <- NormalizeData(seurat_new,
                               normalization.method = 'CLR',
                               margin = 2,
                               verbose = FALSE)

  cat("  复制降维结果...\n")
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

  cat("  保存兼容文件:", COMPATIBLE_FILE, "\n")
  saveRDS(seurat_new, COMPATIBLE_FILE)

  seurat_obj <- seurat_new
  cat("  ✓ 兼容对象创建完成\n\n")
} else {
  cat("步骤 1/3: 加载兼容的 Seurat 对象...\n")
  cat("  文件:", COMPATIBLE_FILE, "\n")
  seurat_obj <- readRDS(COMPATIBLE_FILE)
  cat("  ✓ 加载成功\n\n")
}

# 步骤 2：下采样（如果启用快速模式）
if (USE_FAST_MODE && ncol(seurat_obj) > SAMPLE_SIZE) {
  cat("步骤 2/3: 下采样以提高响应速度...\n")
  cat("  原始细胞数:", ncol(seurat_obj), "\n")
  cat("  目标细胞数:", SAMPLE_SIZE, "\n")

  # 智能采样：保留每种细胞类型的代表性
  if ("author_annotation" %in% colnames(seurat_obj@meta.data)) {
    cell_types <- seurat_obj$author_annotation
    cells_per_type <- table(cell_types)

    proportions <- cells_per_type / sum(cells_per_type)
    sample_per_type <- round(proportions * SAMPLE_SIZE)
    sample_per_type[sample_per_type < 50] <- min(50, cells_per_type[sample_per_type < 50])

    set.seed(42)
    sampled_cells <- c()
    for (ct in names(sample_per_type)) {
      type_cells <- colnames(seurat_obj)[cell_types == ct]
      n_sample <- min(sample_per_type[ct], length(type_cells))
      sampled_cells <- c(sampled_cells, sample(type_cells, n_sample))
    }

    seurat_obj <- seurat_obj[, sampled_cells]
  } else {
    set.seed(42)
    sampled_cells <- sample(colnames(seurat_obj), SAMPLE_SIZE)
    seurat_obj <- seurat_obj[, sampled_cells]
  }

  cat("  采样后细胞数:", ncol(seurat_obj), "\n")
  cat("  ✓ 下采样完成\n\n")
} else {
  cat("步骤 2/3: 跳过下采样\n")
  cat("  使用完整数据集\n\n")
}

# 步骤 3：启动 CITEViz
cat("步骤 3/3: 启动 CITEViz...\n")
cat("=================================================\n")
cat("数据信息:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数:", nrow(seurat_obj@assays$RNA), "\n")
cat("  ADT 数:", nrow(seurat_obj@assays$ADT), "\n")
cat("  降维:", paste(names(seurat_obj@reductions), collapse=", "), "\n")
cat("=================================================\n\n")

cat("提示:\n")
cat("  - CITEViz 将在浏览器中打开\n")
cat("  - 使用类似 FlowJo 的界面分析 ADT 数据\n")
cat("  - 按 Ctrl+C 停止应用\n\n")

# 尝试不同的启动函数
if (exists("run_app", where = "package:CITEViz")) {
  cat("启动 CITEViz...\n\n")
  run_app(seurat_object = seurat_obj)
} else if (exists("launchApp", where = "package:CITEViz")) {
  cat("启动 CITEViz...\n\n")
  launchApp(seurat_obj)
} else if (exists("CITEViz", where = "package:CITEViz")) {
  cat("启动 CITEViz...\n\n")
  CITEViz(seurat_obj)
} else {
  cat("错误：未找到 CITEViz 启动函数\n")
  cat("可用函数:\n")
  print(ls("package:CITEViz"))
  stop("请检查 CITEViz 包是否正确安装")
}
