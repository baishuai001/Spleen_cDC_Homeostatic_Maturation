#!/usr/bin/env Rscript
# ============================================================================
# Seurat ADT 分析 - CITEViz 替代方案
# ============================================================================
# 功能：使用 Seurat 原生函数进行 ADT（抗体）数据分析
# 优点：无需 Shiny，100% Seurat v5 兼容，生成出版级图表
# 输出：PDF 格式的静态图表
# ============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(pheatmap)

cat("\n", strrep("=", 80), "\n")
cat("Seurat ADT 分析（CITEViz 替代方案）\n")
cat(strrep("=", 80), "\n\n")

# ============================================================================
# 1. 自动检测输入文件
# ============================================================================
cat("检测 Seurat 对象文件...\n")

possible_input_files <- c(
  "results/integrated/Robjects/seurat_obj_annotated.rds",
  "results/integrated/Robjects/seurat_obj_for_annotation.rds",
  "results/cDC1/Robjects/seurat_obj_annotated.rds",
  "results/integrated/seurat_obj_annotated.rds",
  "seurat_obj_annotated.rds",
  "seurat_obj.rds"
)

INPUT_FILE <- NULL
for (file_path in possible_input_files) {
  if (file.exists(file_path)) {
    INPUT_FILE <- file_path
    cat("  ✓ 找到文件:", file_path, "\n")
    break
  }
}

if (is.null(INPUT_FILE)) {
  cat("\n✗ 错误：未找到 Seurat 对象文件\n\n")
  cat("查找的位置:\n")
  for (file_path in possible_input_files) {
    cat("  -", file_path, "\n")
  }
  cat("\n")
  cat("解决方法:\n")
  cat("  1. 运行下游分析生成 Seurat 对象:\n")
  cat("     rmarkdown::render('4.downstream_analysis.Rmd')\n\n")
  cat("  2. 或运行文件检测脚本:\n")
  cat("     Rscript find_seurat_object.R\n\n")
  cat("  3. 或手动指定文件路径:\n")
  cat("     修改本脚本第 28 行的 possible_input_files 列表\n\n")
  stop("缺少必需的输入文件")
}

cat("\n加载 Seurat 对象...\n")
seurat_obj <- readRDS(INPUT_FILE)
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数:", nrow(seurat_obj), "\n")
cat("  Assay:", paste(names(seurat_obj@assays), collapse=", "), "\n\n")

# 检查是否有 ADT 数据
if (!"ADT" %in% names(seurat_obj@assays)) {
  cat("✗ 错误：此 Seurat 对象没有 ADT 数据\n")
  cat("  可用的 Assay:", paste(names(seurat_obj@assays), collapse=", "), "\n\n")
  stop("缺少 ADT 数据")
}

# 输出目录
output_dir <- "results/integrated/adt_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("输出目录:", output_dir, "\n")
cat(strrep("=", 80), "\n\n")

# 1. ADT 双参数散点图（类似 FlowJo）
cat("1. 生成 ADT 双参数散点图...\n")

DefaultAssay(seurat_obj) <- "ADT"

# 常见的 cDC 标记组合
adt_pairs <- list(
  c("adt-CD11c", "adt-CD317"),    # cDC1 vs pDC
  c("adt-CD11c", "adt-Sirpa"),    # cDC1 vs cDC2
  c("adt-CD80", "adt-CD86"),      # 成熟标记
  c("adt-MHC-II", "adt-CD11c")    # 呈递能力
)

for (pair in adt_pairs) {
  if (all(pair %in% rownames(seurat_obj@assays$ADT))) {
    p <- FeatureScatter(seurat_obj,
                        feature1 = pair[1],
                        feature2 = pair[2],
                        group.by = "author_annotation",
                        pt.size = 0.5) +
      ggtitle(paste(pair[1], "vs", pair[2]))

    filename <- paste0(output_dir, "/scatter_",
                       gsub("adt-", "", pair[1]), "_vs_",
                       gsub("adt-", "", pair[2]), ".png")
    ggsave(filename, p, width = 10, height = 8, dpi = 150)
    cat("  保存:", filename, "\n")
  }
}

# 2. ADT 在 UMAP 上的表达
cat("\n2. 生成 ADT-UMAP 叠加图...\n")

DefaultAssay(seurat_obj) <- "ADT"
key_adts <- head(rownames(seurat_obj@assays$ADT), 12)  # 前12个

# 检测可用的 UMAP
available_umaps <- names(seurat_obj@reductions)
cat("  检测到的降维:", paste(available_umaps, collapse=", "), "\n")

# 选择 UMAP（优先顺序）
umap_name <- NULL
for (name in c("wnn.umap", "wnn_umap", "umap")) {
  if (name %in% available_umaps) {
    umap_name <- name
    break
  }
}

if (is.null(umap_name)) {
  cat("  警告：未找到 UMAP，跳过此步骤\n")
} else {
  cat("  使用 UMAP:", umap_name, "\n")

  plots <- lapply(key_adts, function(adt) {
    tryCatch({
      FeaturePlot(seurat_obj, features = adt,
                  reduction = umap_name,
                  pt.size = 0.5) +
        scale_color_viridis_c(option = "magma")
    }, error = function(e) {
      cat("    跳过", adt, ":", e$message, "\n")
      NULL
    })
  })

  # 过滤掉 NULL
  plots <- plots[!sapply(plots, is.null)]

  if (length(plots) > 0) {
    p_combined <- wrap_plots(plots, ncol = 4)
    ggsave(paste0(output_dir, "/adt_umap_overlay.png"),
           p_combined, width = 16, height = 12, dpi = 150)
    cat("  保存:", paste0(output_dir, "/adt_umap_overlay.png"), "\n")
  }
}

# 3. ADT 热图（各细胞类型的平均表达）
cat("\n3. 生成 ADT 热图...\n")

library(pheatmap)

# 计算平均表达
DefaultAssay(seurat_obj) <- "ADT"
avg_adt <- AverageExpression(seurat_obj,
                              group.by = "author_annotation",
                              assays = "ADT")$ADT

# 检查并处理 NA/Inf 值
if (any(is.na(avg_adt)) || any(is.infinite(avg_adt))) {
  cat("  检测到 NA/Inf 值，进行清理...\n")
  avg_adt[is.na(avg_adt)] <- 0
  avg_adt[is.infinite(avg_adt)] <- 0
}

# 移除全为 0 的行
avg_adt <- avg_adt[rowSums(abs(avg_adt)) > 0, , drop = FALSE]

# 绘制热图
if (nrow(avg_adt) > 0 && ncol(avg_adt) > 0) {
  png(paste0(output_dir, "/adt_heatmap_by_celltype.png"),
      width = 10, height = 8, units = "in", res = 150)
  pheatmap(avg_adt,
           scale = "row",
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "ADT Expression by Cell Type")
  dev.off()
  cat("  保存:", paste0(output_dir, "/adt_heatmap_by_celltype.png"), "\n")
} else {
  cat("  警告：数据不足，跳过热图生成\n")
}

# 4. 门控分析示例（基于阈值）
cat("\n4. 执行门控分析...\n")

# 示例：CD11c+ CD317- 定义 cDC1
if (all(c("adt-CD11c", "adt-CD317") %in% rownames(seurat_obj@assays$ADT))) {
  # Seurat v5: 使用 layer 而不是 slot
  cd11c <- GetAssayData(seurat_obj, assay = "ADT", layer = "data")["adt-CD11c", ]
  cd317 <- GetAssayData(seurat_obj, assay = "ADT", layer = "data")["adt-CD317", ]

  # 定义门控
  seurat_obj$gated_cDC1 <- ifelse(cd11c > 2 & cd317 < 1, "cDC1_gated", "Other")

  # 可视化门控结果
  p1 <- FeatureScatter(seurat_obj,
                       feature1 = "adt-CD11c",
                       feature2 = "adt-CD317",
                       group.by = "gated_cDC1",
                       pt.size = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
    ggtitle("Manual Gating: CD11c+ CD317-")

  # 使用之前检测到的 UMAP 名称
  if (!is.null(umap_name)) {
    p2 <- DimPlot(seurat_obj, group.by = "gated_cDC1",
                  reduction = umap_name) +
      ggtitle("Gated cDC1 on UMAP")

    p_gating <- p1 + p2
    ggsave(paste0(output_dir, "/manual_gating_example.png"),
           p_gating, width = 14, height = 6, dpi = 150)
    cat("  保存:", paste0(output_dir, "/manual_gating_example.png"), "\n")
  } else {
    # 只保存散点图
    ggsave(paste0(output_dir, "/manual_gating_example.png"),
           p1, width = 10, height = 8, dpi = 150)
    cat("  保存:", paste0(output_dir, "/manual_gating_example.png"), "\n")
  }
}

# 5. ADT 相关性矩阵
cat("\n5. 计算 ADT 相关性...\n")

# Seurat v5: 使用 layer 而不是 slot
adt_mat <- GetAssayData(seurat_obj, assay = "ADT", layer = "data")

# 转换为矩阵并计算相关性
adt_mat_dense <- as.matrix(adt_mat)

# 处理 NA 值
if (any(is.na(adt_mat_dense))) {
  cat("  检测到 NA 值，进行清理...\n")
  adt_mat_dense[is.na(adt_mat_dense)] <- 0
}

adt_cor <- cor(t(adt_mat_dense), method = "spearman", use = "pairwise.complete.obs")

# 检查相关性矩阵
if (!any(is.na(adt_cor)) && !any(is.infinite(adt_cor))) {
  png(paste0(output_dir, "/adt_correlation.png"),
      width = 10, height = 10, units = "in", res = 150)
  pheatmap(adt_cor,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "ADT Correlation Matrix",
           fontsize_row = 8,
           fontsize_col = 8)
  dev.off()
  cat("  保存:", paste0(output_dir, "/adt_correlation.png"), "\n")
} else {
  cat("  警告：相关性矩阵包含 NA/Inf，跳过此步骤\n")
}

cat("\n✓ 分析完成！\n")
cat("所有结果保存在:", output_dir, "\n")
cat("\n提示：这些静态图分析速度快，不卡顿\n")
cat("如需交互式分析，使用 run_citeviz_fast.R（下采样版本）\n")
