#!/usr/bin/env Rscript
# 使用 Seurat 进行 ADT 分析 - 更快的替代方案

library(Seurat)
library(ggplot2)
library(patchwork)

# 加载数据
seurat_obj <- readRDS("results/cDC1/Robjects/seurat_obj_annotated.rds")

# 输出目录
output_dir <- "results/cDC1/ADT_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("开始 ADT 分析（Seurat 版本 - 无需 Shiny）\n")
cat("输出目录:", output_dir, "\n\n")

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

plots <- lapply(key_adts, function(adt) {
  FeaturePlot(seurat_obj, features = adt,
              reduction = "wnn.umap",
              pt.size = 0.5) +
    scale_color_viridis_c(option = "magma")
})

p_combined <- wrap_plots(plots, ncol = 4)
ggsave(paste0(output_dir, "/adt_umap_overlay.png"),
       p_combined, width = 16, height = 12, dpi = 150)
cat("  保存:", paste0(output_dir, "/adt_umap_overlay.png\n"))

# 3. ADT 热图（各细胞类型的平均表达）
cat("\n3. 生成 ADT 热图...\n")

library(pheatmap)

# 计算平均表达
DefaultAssay(seurat_obj) <- "ADT"
avg_adt <- AverageExpression(seurat_obj,
                              group.by = "author_annotation",
                              assays = "ADT")$ADT

# 绘制热图
png(paste0(output_dir, "/adt_heatmap_by_celltype.png"),
    width = 10, height = 8, units = "in", res = 150)
pheatmap(avg_adt,
         scale = "row",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "ADT Expression by Cell Type")
dev.off()
cat("  保存:", paste0(output_dir, "/adt_heatmap_by_celltype.png\n"))

# 4. 门控分析示例（基于阈值）
cat("\n4. 执行门控分析...\n")

# 示例：CD11c+ CD317- 定义 cDC1
if (all(c("adt-CD11c", "adt-CD317") %in% rownames(seurat_obj@assays$ADT))) {
  cd11c <- GetAssayData(seurat_obj, assay = "ADT", slot = "data")["adt-CD11c", ]
  cd317 <- GetAssayData(seurat_obj, assay = "ADT", slot = "data")["adt-CD317", ]

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

  p2 <- DimPlot(seurat_obj, group.by = "gated_cDC1",
                reduction = "wnn.umap") +
    ggtitle("Gated cDC1 on UMAP")

  p_gating <- p1 + p2
  ggsave(paste0(output_dir, "/manual_gating_example.png"),
         p_gating, width = 14, height = 6, dpi = 150)
  cat("  保存:", paste0(output_dir, "/manual_gating_example.png\n"))
}

# 5. ADT 相关性矩阵
cat("\n5. 计算 ADT 相关性...\n")

adt_mat <- GetAssayData(seurat_obj, assay = "ADT", slot = "data")
adt_cor <- cor(t(as.matrix(adt_mat)), method = "spearman")

png(paste0(output_dir, "/adt_correlation.png"),
    width = 10, height = 10, units = "in", res = 150)
pheatmap(adt_cor,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "ADT Correlation Matrix",
         fontsize_row = 8,
         fontsize_col = 8)
dev.off()
cat("  保存:", paste0(output_dir, "/adt_correlation.png\n"))

cat("\n✓ 分析完成！\n")
cat("所有结果保存在:", output_dir, "\n")
cat("\n提示：这些静态图分析速度快，不卡顿\n")
cat("如需交互式分析，使用 run_citeviz_fast.R（下采样版本）\n")
