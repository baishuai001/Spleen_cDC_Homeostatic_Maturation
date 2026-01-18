################################################################################
#                    细胞类型标注辅助脚本                                        #
#                    基于文献标准的自动化标注建议                                  #
################################################################################
#
# 使用说明：
# 1. 确保已运行完 4.downstream_analysis.Rmd
# 2. 修改下方的 analysis_mode 和 output_dir
# 3. 运行脚本获取标注建议
#
# 作者: baishuai
# 参考: GSE228544 原文献标注标准
################################################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
})

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║          细胞类型标注辅助工具 v1.0                         ║\n")
cat("║          基于 GSE228544 原文献标准                         ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# ============================================================================ #
#                          配置区域（请修改）                                    #
# ============================================================================ #

# ★★★ 选择分析模式 ★★★
analysis_mode <- "integrated"  # "cDC1", "cDC2", 或 "integrated"

# ★★★ 设置路径 ★★★
output_dir <- paste0("downstream_results/", analysis_mode)

# ============================================================================ #
#                          标志基因定义（文献标准）                               #
# ============================================================================ #

# cDC1 vs cDC2 谱系特异性标志
lineage_markers <- list(
  cDC1 = c("Xcr1", "Clec9a", "Itgae", "Irf8", "Batf3", "Id2", "Cadm1", "Tlr3"),
  cDC2 = c("Sirpa", "Cd24a", "Notch2", "Irf4", "Klf4", "Esam", "Cd209a", "Clec4a4")
)

# 发育阶段标志
stage_markers <- list(
  # 前体阶段
  precursor = c("Kit", "Flt3", "Ly6d"),

  # 增殖阶段
  proliferating = c("Mki67", "Top2a", "Pcna", "Ccnb1", "Cdk1", "Mcm5"),

  # 未成熟阶段（低成熟markers）
  immature = c("H2-Aa", "H2-Ab1", "Cd74"),  # MHCII基础表达

  # 成熟阶段
  mature = c("Ccr7", "Fscn1", "Relb", "Cd40", "Cd86", "Cd80", "Il12b", "Socs2"),

  # LXR相关（文献核心通路）
  lxr_pathway = c("Nr1h3", "Nr1h2", "Abca1", "Abcg1", "Srebf1", "Fasn")
)

# ============================================================================ #
#                          加载数据                                             #
# ============================================================================ #

cat("════════════════════════════════════════════════════════════\n")
cat("【第1步】加载数据\n")
cat("════════════════════════════════════════════════════════════\n\n")

# 加载Seurat对象
seurat_file <- file.path(output_dir, "Robjects", "seurat_obj_for_annotation.rds")
if (!file.exists(seurat_file)) {
  stop("找不到Seurat对象文件: ", seurat_file,
       "\n请确认已运行 4.downstream_analysis.Rmd")
}

seurat_obj <- readRDS(seurat_file)
cat("Seurat对象加载成功\n")
cat("  细胞数:", ncol(seurat_obj), "\n")

# 加载marker文件
markers_file <- file.path(output_dir, "Tables", "RNA_markers_for_RNA_clusters.xlsx")
if (!file.exists(markers_file)) {
  stop("找不到marker文件: ", markers_file)
}

# 读取所有sheet
sheet_names <- getSheetNames(markers_file)
all_markers <- lapply(sheet_names, function(sheet) {
  read.xlsx(markers_file, sheet = sheet)
})
names(all_markers) <- sheet_names

cat("  Marker文件加载成功，共", length(sheet_names), "个clusters\n\n")

# ============================================================================ #
#                          分析每个Cluster                                      #
# ============================================================================ #

cat("════════════════════════════════════════════════════════════\n")
cat("【第2步】分析每个Cluster的特征\n")
cat("════════════════════════════════════════════════════════════\n\n")

# 分析函数
analyze_cluster <- function(cluster_markers, cluster_id) {

  result <- list(
    cluster = cluster_id,
    cell_count = sum(seurat_obj$SCT_clusters == cluster_id),
    lineage = NA,
    lineage_score = 0,
    stage = NA,
    stage_evidence = "",
    suggested_annotation = NA,
    confidence = "Low",
    top_markers = ""
  )

  if (nrow(cluster_markers) == 0) {
    return(result)
  }

  # 获取top markers
  top_genes <- head(cluster_markers$gene, 20)
  result$top_markers <- paste(head(cluster_markers$gene, 5), collapse = ", ")

  # ==================== 判断谱系 ====================
  cDC1_markers_found <- intersect(top_genes, lineage_markers$cDC1)
  cDC2_markers_found <- intersect(top_genes, lineage_markers$cDC2)

  cDC1_score <- length(cDC1_markers_found)
  cDC2_score <- length(cDC2_markers_found)

  if (cDC1_score > cDC2_score && cDC1_score >= 2) {
    result$lineage <- "cDC1"
    result$lineage_score <- cDC1_score
  } else if (cDC2_score > cDC1_score && cDC2_score >= 2) {
    result$lineage <- "cDC2"
    result$lineage_score <- cDC2_score
  } else if (cDC1_score == 1 && cDC2_score == 0) {
    result$lineage <- "cDC1 (可能)"
    result$lineage_score <- cDC1_score
  } else if (cDC2_score == 1 && cDC1_score == 0) {
    result$lineage <- "cDC2 (可能)"
    result$lineage_score <- cDC2_score
  } else {
    result$lineage <- "不确定"
    result$lineage_score <- 0
  }

  # ==================== 判断发育阶段 ====================
  prolif_found <- intersect(top_genes, stage_markers$proliferating)
  mature_found <- intersect(top_genes, stage_markers$mature)
  precursor_found <- intersect(top_genes, stage_markers$precursor)

  evidence_parts <- c()

  # 判断逻辑
  if (length(prolif_found) >= 2) {
    result$stage <- "Proliferating"
    evidence_parts <- c(evidence_parts, paste("增殖:", paste(prolif_found, collapse = "+")))
  } else if (length(mature_found) >= 3) {
    result$stage <- "Late Mature"
    evidence_parts <- c(evidence_parts, paste("成熟:", paste(mature_found, collapse = "+")))
  } else if (length(mature_found) >= 1) {
    result$stage <- "Early Mature"
    evidence_parts <- c(evidence_parts, paste("成熟:", paste(mature_found, collapse = "+")))
  } else if (length(precursor_found) >= 1) {
    result$stage <- "Pre-cDC"
    evidence_parts <- c(evidence_parts, paste("前体:", paste(precursor_found, collapse = "+")))
  } else {
    # 默认为未成熟
    if (any(c("Ccr7", "Fscn1") %in% top_genes)) {
      result$stage <- "Immature (低成熟)"
    } else {
      result$stage <- "Early Immature"
    }
  }

  result$stage_evidence <- paste(evidence_parts, collapse = "; ")

  # ==================== 生成标注建议 ====================
  if (!is.na(result$lineage) && result$lineage != "不确定") {
    base_lineage <- gsub(" \\(可能\\)", "", result$lineage)
    if (!is.na(result$stage)) {
      result$suggested_annotation <- paste(result$stage, base_lineage)
    } else {
      result$suggested_annotation <- base_lineage
    }
  } else {
    result$suggested_annotation <- paste("Unknown -", result$stage)
  }

  # ==================== 判断置信度 ====================
  if (result$lineage_score >= 3 && !is.na(result$stage)) {
    result$confidence <- "High"
  } else if (result$lineage_score >= 2) {
    result$confidence <- "Medium"
  } else {
    result$confidence <- "Low"
  }

  return(result)
}

# 分析所有cluster
results <- lapply(names(all_markers), function(cluster_id) {
  analyze_cluster(all_markers[[cluster_id]], cluster_id)
})

results_df <- do.call(rbind, lapply(results, as.data.frame))

# ============================================================================ #
#                          显示结果                                             #
# ============================================================================ #

cat("════════════════════════════════════════════════════════════\n")
cat("【分析结果】标注建议\n")
cat("════════════════════════════════════════════════════════════\n\n")

# 按置信度排序显示
results_df <- results_df[order(results_df$confidence, decreasing = TRUE), ]

for (i in 1:nrow(results_df)) {
  row <- results_df[i, ]

  # 置信度颜色标记
  conf_mark <- switch(as.character(row$confidence),
                      "High" = "[高]",
                      "Medium" = "[中]",
                      "Low" = "[低]")

  cat("────────────────────────────────────────────────────────────\n")
  cat("Cluster", row$cluster, conf_mark, "\n")
  cat("────────────────────────────────────────────────────────────\n")
  cat("  细胞数:     ", row$cell_count, "\n")
  cat("  谱系:       ", row$lineage, "(得分:", row$lineage_score, ")\n")
  cat("  发育阶段:   ", row$stage, "\n")
  if (row$stage_evidence != "") {
    cat("  判断依据:   ", row$stage_evidence, "\n")
  }
  cat("  Top markers:", row$top_markers, "\n")
  cat("  ★ 建议标注: ", row$suggested_annotation, "\n")
  cat("\n")
}

# ============================================================================ #
#                          生成汇总表                                           #
# ============================================================================ #

cat("════════════════════════════════════════════════════════════\n")
cat("【汇总表】\n")
cat("════════════════════════════════════════════════════════════\n\n")

summary_df <- results_df[, c("cluster", "cell_count", "suggested_annotation",
                             "confidence", "top_markers")]
colnames(summary_df) <- c("Cluster", "细胞数", "建议标注", "置信度", "Top Markers")

# 恢复原始顺序
summary_df <- summary_df[order(as.numeric(summary_df$Cluster)), ]
print(summary_df, row.names = FALSE)

# 保存到Excel
output_file <- file.path(output_dir, "Tables", "annotation_suggestions.xlsx")
write.xlsx(summary_df, output_file)
cat("\n汇总表已保存到:", output_file, "\n")

# ============================================================================ #
#                          生成验证图                                           #
# ============================================================================ #

cat("\n")
cat("════════════════════════════════════════════════════════════\n")
cat("【第3步】生成验证图\n")
cat("════════════════════════════════════════════════════════════\n\n")

# 检查关键基因的表达
key_genes <- c(
  # 谱系markers
  "Xcr1", "Clec9a", "Sirpa", "Cd24a",
  # 成熟markers
  "Ccr7", "Fscn1", "Cd86", "Cd40",
  # 增殖markers
  "Mki67", "Top2a"
)

# 过滤存在的基因
existing_genes <- key_genes[key_genes %in% rownames(seurat_obj)]
cat("检测到的关键基因:", paste(existing_genes, collapse = ", "), "\n\n")

if (length(existing_genes) > 0) {
  # DotPlot
  p1 <- DotPlot(seurat_obj, features = existing_genes, group.by = "SCT_clusters") +
    RotatedAxis() +
    labs(title = paste0("关键Marker基因表达 (", analysis_mode, ")"),
         subtitle = "用于验证标注建议")

  # 保存
  ggsave(file.path(output_dir, "Plots", "annotation_validation_dotplot.png"),
         p1, width = 14, height = 8, dpi = 300)
  cat("验证图已保存: annotation_validation_dotplot.png\n")

  print(p1)
}

# ============================================================================ #
#                          文献标准参考                                          #
# ============================================================================ #

cat("\n")
cat("════════════════════════════════════════════════════════════\n")
cat("【文献标准参考】GSE228544 细胞类型定义\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("cDC1 发育轨迹（6个阶段）:\n")
cat("┌─────────────────────────────────────────────────────────────────┐\n")
cat("│ 阶段              │ 关键特征                                    │\n")
cat("├─────────────────────────────────────────────────────────────────┤\n")
cat("│ Pre-cDC1          │ Xcr1+, Clec9a低, Kit可能+, Flt3+           │\n")
cat("│ Proliferating     │ Xcr1+, Mki67高, Top2a高, Pcna高             │\n")
cat("│ Early Immature    │ Xcr1+, Clec9a+, Ccr7低, Cd86低              │\n")
cat("│ Late Immature     │ Xcr1+, Clec9a+, Ccr7开始升高                │\n")
cat("│ Early Mature      │ Xcr1+, Ccr7+, Fscn1+, Cd86+                 │\n")
cat("│ Late Mature       │ Xcr1+, Ccr7高, Relb高, Il12b可能+           │\n")
cat("└─────────────────────────────────────────────────────────────────┘\n\n")

cat("cDC2 发育轨迹（6个阶段）:\n")
cat("┌─────────────────────────────────────────────────────────────────┐\n")
cat("│ 阶段              │ 关键特征                                    │\n")
cat("├─────────────────────────────────────────────────────────────────┤\n")
cat("│ Pre-cDC2          │ Sirpa+, Cd24a+, 低成熟markers               │\n")
cat("│ Proliferating     │ Sirpa+, Mki67高, Top2a高                    │\n")
cat("│ Early Immature    │ Sirpa+, Cd24a+, Irf4+, Ccr7低               │\n")
cat("│ Late Immature     │ Sirpa+, Ccr7开始升高, Cd86低                │\n")
cat("│ Early Mature      │ Sirpa+, Ccr7+, Fscn1+, Cd86+                │\n")
cat("│ Late Mature       │ Sirpa+, Ccr7高, Relb高, Cd40高              │\n")
cat("└─────────────────────────────────────────────────────────────────┘\n\n")

# ============================================================================ #
#                          生成可复制的标注向量                                   #
# ============================================================================ #

cat("════════════════════════════════════════════════════════════\n")
cat("【可复制代码】标注向量模板\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("# 复制以下代码到您的标注脚本中，根据验证结果修改\n\n")
cat("cluster_annotations <- c(\n")

for (i in 1:nrow(summary_df)) {
  row <- summary_df[i, ]
  annotation <- gsub("\"", "'", row$`建议标注`)
  comma <- if(i < nrow(summary_df)) "," else ""
  cat(sprintf('  "%s" = "%s"%s  # %s\n',
              row$Cluster,
              annotation,
              comma,
              row$`Top Markers`))
}
cat(")\n\n")

# ============================================================================ #
#                          下一步指南                                           #
# ============================================================================ #

cat("════════════════════════════════════════════════════════════\n")
cat("【下一步】\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("1. 查看生成的 annotation_validation_dotplot.png\n")
cat("   验证谱系markers（Xcr1 vs Sirpa）和成熟markers（Ccr7, Cd86）\n\n")

cat("2. 对比文献补充材料的标准\n")
cat("   确认您的cluster与文献描述是否一致\n\n")

cat("3. 修改标注向量\n")
cat("   根据验证结果调整建议标注\n\n")

cat("4. 合并相似clusters（如果需要）\n")
cat("   例如：多个Early Immature可以合并\n\n")

cat("5. 应用标注到Seurat对象\n")
cat("   使用 Cell_Type_Annotation_Script.R 中的方法\n\n")

cat("════════════════════════════════════════════════════════════\n")
cat("完成！请根据上述建议进行手动验证和调整。\n")
cat("════════════════════════════════════════════════════════════\n\n")
