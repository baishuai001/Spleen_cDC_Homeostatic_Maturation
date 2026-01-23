################################################################################
## 简单的R环境检查脚本
################################################################################

cat("=================================================================\n")
cat("R环境快速检查\n")
cat("=================================================================\n\n")

# 1. R版本
cat("1. R版本:\n")
cat("   ", R.version.string, "\n\n")

# 2. Matrix包
cat("2. Matrix包:\n")
if (requireNamespace("Matrix", quietly = TRUE)) {
  cat("   已安装版本:", as.character(packageVersion("Matrix")), "\n")
  cat("   需要版本: >= 1.6.3\n")

  current_ver <- as.character(packageVersion("Matrix"))
  if (compareVersion(current_ver, "1.6.3") >= 0) {
    cat("   状态: ✓ 版本符合要求\n")
  } else {
    cat("   状态: ✗ 版本过低，需要升级\n")
  }
} else {
  cat("   状态: ✗ 未安装\n")
}
cat("\n")

# 3. Seurat包
cat("3. Seurat包:\n")
tryCatch({
  library(Seurat)
  cat("   已安装版本:", as.character(packageVersion("Seurat")), "\n")
  cat("   状态: ✓ 可以加载\n")
}, error = function(e) {
  cat("   状态: ✗ 无法加载\n")
  cat("   错误:", e$message, "\n")
})
cat("\n")

# 4. 安装路径
cat("4. R包安装路径:\n")
lib_paths <- .libPaths()
for (path in lib_paths) {
  writable <- file.access(path, 2) == 0
  cat("   ", path, ifelse(writable, "(可写)", "(只读)"), "\n")
}
cat("\n")

# 5. 快速建议
cat("=================================================================\n")
cat("快速建议:\n")
cat("=================================================================\n\n")

# 检查Matrix
if (requireNamespace("Matrix", quietly = TRUE)) {
  if (compareVersion(as.character(packageVersion("Matrix")), "1.6.3") < 0) {
    cat("需要升级Matrix包。请在R中运行:\n\n")
    cat('remove.packages("Matrix")\n')
    cat('# 然后退出R (输入 q() 并按回车)\n')
    cat('# 重新启动R后运行:\n')
    cat('install.packages("Matrix")\n\n')

    cat("如果仍然失败，运行:\n")
    cat('source("force_fix_matrix.R")\n\n')
  } else {
    cat("✓ Matrix版本正常，可以使用Seurat进行分析\n\n")
  }
} else {
  cat("需要安装Matrix包:\n")
  cat('install.packages("Matrix")\n\n')
}

cat("=================================================================\n")
