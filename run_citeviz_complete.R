#!/usr/bin/env Rscript
# ============================================================================
# CITEViz 完整工作流（一键运行）
# ============================================================================
# 自动完成：
#   1. Seurat v5 → v4 格式转换
#   2. 启动 CITEViz（预加载数据）
#
# 使用方法：
#   Rscript run_citeviz_complete.R
# ============================================================================

cat("\n", strrep("=", 80), "\n")
cat("CITEViz 完整工作流\n")
cat(strrep("=", 80), "\n\n")

# ============================================================================
# 步骤 1: 格式转换
# ============================================================================
cat("[步骤 1/2] Seurat v5 → v4 格式转换\n")
cat(strrep("-", 80), "\n")

conversion_result <- tryCatch({
  source("convert_seurat_v5_to_v4_for_citeviz.R", echo = FALSE)
  TRUE
}, error = function(e) {
  cat("\n✗ 转换失败:", e$message, "\n")
  FALSE
})

if (!conversion_result) {
  cat("\n建议：\n")
  cat("  1. 单独运行转换脚本查看详细错误:\n")
  cat("     Rscript convert_seurat_v5_to_v4_for_citeviz.R\n")
  cat("  2. 或使用替代分析方案:\n")
  cat("     Rscript adt_analysis_seurat.R\n\n")
  stop("格式转换失败，无法继续")
}

cat("\n✓ 转换完成\n\n")

# 等待一下让用户看到结果
Sys.sleep(2)

# ============================================================================
# 步骤 2: 启动 CITEViz
# ============================================================================
cat("[步骤 2/2] 启动 CITEViz\n")
cat(strrep("-", 80), "\n\n")

cat("正在启动 CITEViz（预加载模式）...\n")
cat("按 Ctrl+C 可随时停止\n\n")

# 启动 CITEViz（预加载数据）
tryCatch({
  system2("Rscript", args = c("launch_citeviz.R", "--preload"),
          stdout = "", stderr = "")
}, interrupt = function(e) {
  cat("\n\n用户中断，退出 CITEViz\n")
}, error = function(e) {
  cat("\n✗ 启动失败:", e$message, "\n\n")
  cat("建议:\n")
  cat("  1. 手动启动: Rscript launch_citeviz.R --preload\n")
  cat("  2. 检查安装: Rscript check_citeviz_installation.R\n")
  cat("  3. 使用替代方案: Rscript adt_analysis_seurat.R\n\n")
})

cat("\n", strrep("=", 80), "\n")
cat("工作流结束\n")
cat(strrep("=", 80), "\n")
