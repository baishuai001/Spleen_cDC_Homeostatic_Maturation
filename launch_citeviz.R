#!/usr/bin/env Rscript
# ============================================================================
# CITEViz 启动脚本（Seurat v5 兼容版）
# ============================================================================
# 使用方法：
#   1. 直接运行（空启动）： Rscript launch_citeviz.R
#   2. 预加载数据启动：     Rscript launch_citeviz.R --preload
#
# 注意：如果遇到 [object Object] 错误，请先运行：
#   Rscript convert_seurat_v5_to_v4_for_citeviz.R
# ============================================================================

library(CITEViz)
library(Seurat)

cat("\n", strrep("=", 80), "\n")
cat("CITEViz 启动器（Seurat v5 兼容版）\n")
cat(strrep("=", 80), "\n\n")

# ============================================================================
# 1. 检查是否预加载数据
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)
PRELOAD_DATA <- "--preload" %in% args

# v4 兼容文件路径（由 convert_seurat_v5_to_v4_for_citeviz.R 生成）
V4_COMPATIBLE_FILE <- "results/integrated/Robjects/seurat_obj_citeviz_v4_compatible.rds"

# ============================================================================
# 2. 检查文件是否存在
# ============================================================================
if (PRELOAD_DATA) {
  if (!file.exists(V4_COMPATIBLE_FILE)) {
    cat("⚠ 错误：未找到 v4 兼容文件\n")
    cat("  期望路径:", V4_COMPATIBLE_FILE, "\n\n")
    cat("请先运行以下命令生成兼容文件:\n")
    cat("  Rscript convert_seurat_v5_to_v4_for_citeviz.R\n\n")
    stop("缺少必需的数据文件")
  }

  cat("模式：预加载数据\n")
  cat("  加载文件:", V4_COMPATIBLE_FILE, "\n")
} else {
  cat("模式：空启动（需要手动上传文件）\n\n")
  cat("推荐使用的文件:\n")
  cat("  1. ", V4_COMPATIBLE_FILE, " (v4 兼容，推荐)\n")
  cat("  2. results/citeviz/citeviz_fast_5k.rds (5K 细胞)\n")
  cat("  3. results/citeviz/citeviz_medium_8k.rds (8K 细胞)\n\n")

  cat("提示：如果遇到 [object Object] 错误，使用预加载模式:\n")
  cat("  Rscript launch_citeviz.R --preload\n\n")
}

# ============================================================================
# 3. 预加载数据（如果需要）
# ============================================================================
seurat_obj <- NULL

if (PRELOAD_DATA) {
  cat("\n加载 Seurat 对象...\n")

  seurat_obj <- readRDS(V4_COMPATIBLE_FILE)

  cat("  ✓ 加载成功\n")
  cat("  细胞数:", ncol(seurat_obj), "\n")
  cat("  Assay:", paste(names(seurat_obj@assays), collapse=", "), "\n")

  # 验证 Assay 类型
  cat("  Assay 类型:\n")
  for (assay_name in names(seurat_obj@assays)) {
    assay_class <- class(seurat_obj@assays[[assay_name]])[1]
    status <- if (assay_class == "Assay5") "⚠ Assay5" else "✓ Assay"
    cat("    ", assay_name, ":", assay_class, "-", status, "\n")
  }
  cat("\n")
}

# ============================================================================
# 4. 启动 CITEViz
# ============================================================================
cat(strrep("=", 80), "\n")
cat("启动 CITEViz...\n")
cat(strrep("=", 80), "\n\n")

cat("提示:\n")
cat("  - CITEViz 将在浏览器中打开\n")
cat("  - 如果没有预加载数据，需要点击 'Browse...' 上传 .rds 文件\n")
cat("  - 使用 Ctrl+C 停止服务器\n\n")

# 尝试不同的启动函数
launch_success <- FALSE

tryCatch({
  if (PRELOAD_DATA && !is.null(seurat_obj)) {
    # 预加载模式：尝试传入 Seurat 对象
    if (exists("run_app", where = "package:CITEViz")) {
      cat("使用 run_app() 启动...\n\n")
      run_app(seurat_obj)
      launch_success <- TRUE
    } else if (exists("launchApp", where = "package:CITEViz")) {
      cat("使用 launchApp() 启动...\n\n")
      launchApp(seurat_obj)
      launch_success <- TRUE
    } else if (exists("CITEViz", where = "package:CITEViz")) {
      cat("使用 CITEViz() 启动...\n\n")
      CITEViz(seurat_obj)
      launch_success <- TRUE
    }
  } else {
    # 空启动模式
    if (exists("run_app", where = "package:CITEViz")) {
      cat("使用 run_app() 启动...\n\n")
      run_app()
      launch_success <- TRUE
    } else if (exists("launchApp", where = "package:CITEViz")) {
      cat("使用 launchApp() 启动...\n\n")
      launchApp()
      launch_success <- TRUE
    } else if (exists("CITEViz", where = "package:CITEViz")) {
      cat("使用 CITEViz() 启动...\n\n")
      CITEViz()
      launch_success <- TRUE
    }
  }
}, error = function(e) {
  cat("\n✗ 启动失败:", e$message, "\n\n")

  if (grepl("object Object|undefined", e$message, ignore.case = TRUE)) {
    cat("这是 Seurat v5 兼容性错误。建议:\n")
    cat("  1. 确保已运行: Rscript convert_seurat_v5_to_v4_for_citeviz.R\n")
    cat("  2. 使用备用方案: Rscript adt_analysis_seurat.R\n\n")
  }
})

if (!launch_success) {
  cat("\n✗ 错误：未找到 CITEViz 启动函数\n\n")
  cat("可用函数列表:\n")
  print(ls("package:CITEViz"))
  cat("\n")
  stop("请检查 CITEViz 包是否正确安装")
}

cat("\nCITEViz 已关闭\n")
