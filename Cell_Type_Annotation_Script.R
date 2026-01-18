################################################################################
#                    细胞类型标注交互式脚本                                      #
#                    Cell Type Annotation Script                                #
################################################################################
#
# 使用说明：
# 1. 在RStudio中打开此脚本
# 2. 按顺序运行每个代码块
# 3. 根据提示查看结果并输入标注
#
# 作者: baishuai
# 日期: 2024
################################################################################

# ============================================================================ #
#                            第0步：设置工作环境                                #
# ============================================================================ #

# 清理环境
rm(list = ls())
gc()

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(RColorBrewer)
  library(writexl)
})

cat("\n", strrep("=", 60), "\n")
cat("           细胞类型标注工作流程 v1.0                    \n")
cat(strrep("=", 60), "\n\n")

# ============================================================================ #
#                            第1步：加载您的Seurat对象                          #
# ============================================================================ #

# ★★★ 请修改这里：填入您的RDS文件路径 ★★★
rds_file <- "your_seurat_object.rds"  # <-- 改成您的文件路径

cat("正在加载Seurat对象:", rds_file, "\n")

# 检查文件是否存在
if (!file.exists(rds_file)) {
  stop("❌ 找不到文件: ", rds_file, "\n",
       "   请检查文件路径是否正确\n",
       "   例如: 'data/GSE228544/seurat_clustered.rds'")
}

# 加载对象
seurat_obj <- readRDS(rds_file)
cat("✓ Seurat对象加载成功！\n\n")

# 显示基本信息
cat("═══════════════════════════════════════════════════════════\n")
cat("数据概览:\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("  细胞总数:", ncol(seurat_obj), "\n")
cat("  基因总数:", nrow(seurat_obj), "\n")
cat("  Assays:  ", paste(Assays(seurat_obj), collapse = ", "), "\n")
cat("  Reductions:", paste(Reductions(seurat_obj), collapse = ", "), "\n")
cat("\n")

# ============================================================================ #
#                            第2步：查看现有聚类结果                            #
# ============================================================================ #

cat("═══════════════════════════════════════════════════════════\n")
cat("现有聚类信息:\n")
cat("═══════════════════════════════════════════════════════════\n")

# 尝试找到聚类列
possible_cluster_cols <- c("seurat_clusters", "SCT_clusters", "RNA_clusters",
                            "ADT_clusters", "wsnn_res.0.8", "wsnn_res.0.6")

available_clusters <- intersect(possible_cluster_cols, colnames(seurat_obj@meta.data))

if (length(available_clusters) > 0) {
  cat("找到的聚类列:\n")
  for (col in available_clusters) {
    n_clusters <- length(unique(seurat_obj@meta.data[[col]]))
    cat("  •", col, ":", n_clusters, "个clusters\n")
  }

  # 使用第一个找到的聚类
  cluster_col <- available_clusters[1]
  cat("\n→ 将使用:", cluster_col, "\n\n")
} else {
  # 显示所有meta.data列让用户选择
  cat("未找到标准聚类列，现有meta.data列:\n")
  print(colnames(seurat_obj@meta.data))
  cat("\n★ 请手动设置 cluster_col 变量 ★\n")
  cluster_col <- "seurat_clusters"  # 默认值
}

# 显示每个cluster的细胞数
cat("每个Cluster的细胞数:\n")
print(table(seurat_obj@meta.data[[cluster_col]]))
cat("\n")

# ============================================================================ #
#                            第3步：可视化聚类结果                              #
# ============================================================================ #

cat("═══════════════════════════════════════════════════════════\n")
cat("生成可视化图片...\n")
cat("═══════════════════════════════════════════════════════════\n")

# 设置颜色
n_clusters <- length(unique(seurat_obj@meta.data[[cluster_col]]))
if (n_clusters <= 12) {
  cluster_colors <- brewer.pal(max(3, n_clusters), "Set3")
} else {
  cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_clusters)
}

# 绘制UMAP
p_umap <- DimPlot(seurat_obj,
                  group.by = cluster_col,
                  reduction = "umap",
                  label = TRUE,
                  label.size = 5,
                  repel = TRUE,
                  cols = cluster_colors) +
  labs(title = paste("UMAP -", cluster_col)) +
  theme(legend.position = "right")

print(p_umap)

cat("✓ UMAP图已显示\n\n")

# ============================================================================ #
#                            第4步：生成Marker基因列表                          #
# ============================================================================ #

cat("═══════════════════════════════════════════════════════════\n")
cat("生成Marker基因列表（这可能需要几分钟）...\n")
cat("═══════════════════════════════════════════════════════════\n")

# 设置Idents
Idents(seurat_obj) <- cluster_col

# 查找所有cluster的markers
all_markers <- FindAllMarkers(seurat_obj,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.5,
                               verbose = TRUE)

cat("✓ Marker基因查找完成！\n\n")

# 显示每个cluster的top markers
cat("═══════════════════════════════════════════════════════════\n")
cat("每个Cluster的Top 5 Marker基因:\n")
cat("═══════════════════════════════════════════════════════════\n")

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

# 逐个cluster显示
for (cl in sort(unique(top_markers$cluster))) {
  cat("\n", strrep("-", 50), "\n")
  cat("Cluster", cl, ":\n")
  cl_markers <- top_markers %>% filter(cluster == cl)
  for (i in 1:nrow(cl_markers)) {
    cat(sprintf("  %d. %s (FC=%.2f, pct.1=%.0f%%, pct.2=%.0f%%)\n",
                i,
                cl_markers$gene[i],
                cl_markers$avg_log2FC[i],
                cl_markers$pct.1[i] * 100,
                cl_markers$pct.2[i] * 100))
  }
}
cat("\n")

# 保存markers到Excel
markers_file <- "cluster_markers.xlsx"
write_xlsx(all_markers, markers_file)
cat("✓ 完整Marker列表已保存到:", markers_file, "\n\n")

# ============================================================================ #
#                            第5步：可视化关键Marker基因                        #
# ============================================================================ #

cat("═══════════════════════════════════════════════════════════\n")
cat("可视化关键DC Marker基因...\n")
cat("═══════════════════════════════════════════════════════════\n")

# cDC相关的关键markers (会自动过滤不存在的基因)
dc_markers <- c(
  # cDC1 markers
  "Xcr1", "Clec9a", "Itgae", "Irf8", "Batf3", "Id2", "Cadm1",
  # cDC2 markers
  "Sirpa", "Cd24a", "Notch2", "Irf4", "Klf4", "Cd209a", "Esam",
  # 成熟/迁移 markers
  "Ccr7", "Fscn1", "Relb", "Cd40", "Cd86", "Cd80",
  # 增殖 markers
  "Mki67", "Top2a", "Pcna",
  # 炎症 markers
  "Il12b", "Tnf", "Il1b",
  # 通用DC markers
  "Itgax", "H2-Ab1"
)

# 过滤存在的基因
existing_markers <- dc_markers[dc_markers %in% rownames(seurat_obj)]
cat("找到", length(existing_markers), "个DC相关marker基因\n\n")

if (length(existing_markers) > 0) {
  # DotPlot
  p_dot <- DotPlot(seurat_obj,
                   features = existing_markers,
                   group.by = cluster_col,
                   cols = c("lightgrey", "red")) +
    RotatedAxis() +
    labs(title = "DC Marker Genes Expression")

  print(p_dot)

  # FeaturePlot (显示前6个)
  markers_to_show <- head(existing_markers, 6)
  p_feature <- FeaturePlot(seurat_obj,
                           features = markers_to_show,
                           reduction = "umap",
                           cols = c("lightgrey", "red"),
                           ncol = 3)
  print(p_feature)
}

# ============================================================================ #
#                            第6步：创建标注 (核心步骤)                          #
# ============================================================================ #

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║           现在开始进行细胞类型标注！                       ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("请根据上面的Marker基因信息，为每个cluster命名。\n")
cat("\n")
cat("常见DC亚型参考:\n")
cat("  • Resident cDC1     : Xcr1+, Clec9a+, Itgae+\n")
cat("  • Migratory cDC1    : Xcr1+, Ccr7+, Fscn1+\n")
cat("  • Resident cDC2     : Sirpa+, Cd24a+, Notch2+\n")
cat("  • Migratory cDC2    : Sirpa+, Ccr7+, Fscn1+\n")
cat("  • CD209a+ cDC2      : Sirpa+, Cd209a+, Esam-\n")
cat("  • Proliferating DC  : Mki67+, Top2a+\n")
cat("  • Inflammatory DC   : Il12b+, Tnf+, Il1b+\n")
cat("  • pre-cDC1/pre-cDC2 : 较低的成熟markers\n")
cat("\n")

# ★★★ 在这里填入您的标注 ★★★
# 根据您的cluster数量和marker基因，修改下面的标注

cluster_annotations <- c(
  "0" = "待标注 - Cluster 0",
  "1" = "待标注 - Cluster 1",
  "2" = "待标注 - Cluster 2",
  "3" = "待标注 - Cluster 3",
  "4" = "待标注 - Cluster 4",
  "5" = "待标注 - Cluster 5",
  "6" = "待标注 - Cluster 6",
  "7" = "待标注 - Cluster 7",
  "8" = "待标注 - Cluster 8",
  "9" = "待标注 - Cluster 9",
  "10" = "待标注 - Cluster 10",
  "11" = "待标注 - Cluster 11",
  "12" = "待标注 - Cluster 12"
  # 如果有更多cluster，继续添加...
)

# ============================================================================ #
#                            第7步：应用标注                                    #
# ============================================================================ #

cat("═══════════════════════════════════════════════════════════\n")
cat("应用细胞类型标注...\n")
cat("═══════════════════════════════════════════════════════════\n")

# 获取当前clusters
current_clusters <- as.character(seurat_obj@meta.data[[cluster_col]])

# 映射到细胞类型
seurat_obj$cell_type <- sapply(current_clusters, function(x) {
  if (x %in% names(cluster_annotations)) {
    return(cluster_annotations[[x]])
  } else {
    return(paste0("Unknown_", x))
  }
})

# 显示标注结果
cat("\n标注结果统计:\n")
print(table(seurat_obj$cell_type))
cat("\n")

# ============================================================================ #
#                            第8步：可视化标注结果                              #
# ============================================================================ #

cat("═══════════════════════════════════════════════════════════\n")
cat("生成标注后的可视化...\n")
cat("═══════════════════════════════════════════════════════════\n")

# 标注前后对比图
p1 <- DimPlot(seurat_obj,
              group.by = cluster_col,
              reduction = "umap",
              label = TRUE,
              label.size = 4) +
  labs(title = "原始Clusters") +
  NoLegend()

p2 <- DimPlot(seurat_obj,
              group.by = "cell_type",
              reduction = "umap",
              label = TRUE,
              label.size = 3,
              repel = TRUE) +
  labs(title = "细胞类型标注") +
  theme(legend.position = "right")

# 组合图
p_compare <- p1 | p2
print(p_compare)

# 保存对比图
ggsave("annotation_comparison.png", p_compare, width = 18, height = 8, dpi = 300)
cat("✓ 对比图已保存: annotation_comparison.png\n")

# ============================================================================ #
#                            第9步：保存标注后的对象                            #
# ============================================================================ #

cat("\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("保存标注后的Seurat对象...\n")
cat("═══════════════════════════════════════════════════════════\n")

# 保存对象
output_file <- "seurat_annotated.rds"
saveRDS(seurat_obj, output_file)
cat("✓ 已保存到:", output_file, "\n\n")

# ============================================================================ #
#                            第10步：生成标注报告                               #
# ============================================================================ #

cat("═══════════════════════════════════════════════════════════\n")
cat("生成标注报告...\n")
cat("═══════════════════════════════════════════════════════════\n")

# 创建统计表
annotation_stats <- as.data.frame(table(seurat_obj$cell_type))
colnames(annotation_stats) <- c("Cell_Type", "Count")
annotation_stats$Percentage <- round(annotation_stats$Count / sum(annotation_stats$Count) * 100, 2)
annotation_stats <- annotation_stats[order(-annotation_stats$Count), ]

# 打印表格
cat("\n细胞类型统计:\n")
print(annotation_stats, row.names = FALSE)

# 保存到Excel
write_xlsx(annotation_stats, "cell_type_statistics.xlsx")
cat("\n✓ 统计表已保存: cell_type_statistics.xlsx\n")

# ============================================================================ #
#                            完成！                                            #
# ============================================================================ #

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║                    标注流程完成！                          ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("生成的文件:\n")
cat("  1. cluster_markers.xlsx      - 所有cluster的marker基因\n")
cat("  2. annotation_comparison.png - 标注前后对比图\n")
cat("  3. seurat_annotated.rds      - 标注后的Seurat对象\n")
cat("  4. cell_type_statistics.xlsx - 细胞类型统计表\n")
cat("\n")
cat("下一步建议:\n")
cat("  • 检查标注结果是否合理\n")
cat("  • 如需修改，编辑 cluster_annotations 向量后重新运行第7-10步\n")
cat("  • 进行下游分析（差异表达、富集分析等）\n")
cat("\n")
