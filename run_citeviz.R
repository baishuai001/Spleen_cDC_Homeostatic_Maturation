#!/usr/bin/env Rscript
# CITEViz 启动脚本 - 直接加载服务器上的数据
#
# 使用方法：
#   Rscript run_citeviz.R
# 或在 R 中：
#   source("run_citeviz.R")

library(Seurat)
library(CITEViz)

cat("=================================================\n")
cat("CITEViz 启动脚本\n")
cat("=================================================\n\n")

# 数据路径（服务器本地路径）
DATA_PATH <- "results/cDC1/Robjects/seurat_obj_annotated.rds"

# 检查文件是否存在
if (!file.exists(DATA_PATH)) {
  stop("错误：数据文件不存在 - ", DATA_PATH, "\n",
       "请确认路径是否正确")
}

# 加载数据
cat("正在加载 Seurat 对象...\n")
cat("路径:", DATA_PATH, "\n")

seurat_obj <- readRDS(DATA_PATH)

cat("✓ 加载成功\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数:", nrow(seurat_obj), "\n")
cat("  Assays:", paste(Assays(seurat_obj), collapse=", "), "\n")
cat("  降维:", paste(names(seurat_obj@reductions), collapse=", "), "\n")

# 检查必要的组分
if (!"ADT" %in% Assays(seurat_obj)) {
  warning("警告：未找到 ADT assay，CITEViz 可能无法正常工作")
}

if (length(seurat_obj@reductions) == 0) {
  warning("警告：未找到降维结果，建议先运行 UMAP")
}

# 启动 CITEViz
cat("\n启动 CITEViz...\n")
cat("提示：\n")
cat("  - 应用会在浏览器中自动打开\n")
cat("  - 如果没有自动打开，请访问控制台显示的 URL\n")
cat("  - 按 Ctrl+C 停止应用\n\n")

# 尝试不同的启动函数名（因为不同版本可能不同）
if (exists("run_app", where = "package:CITEViz")) {
  run_app(seurat_object = seurat_obj)
} else if (exists("launchApp", where = "package:CITEViz")) {
  launchApp(seurat_obj)
} else if (exists("CITEViz_app", where = "package:CITEViz")) {
  CITEViz_app(seurat_obj)
} else {
  cat("可用的 CITEViz 函数:\n")
  print(ls("package:CITEViz"))
  stop("未找到 CITEViz 启动函数，请查看上面的可用函数列表")
}
