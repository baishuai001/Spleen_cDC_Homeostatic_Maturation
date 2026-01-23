################################################################################
## 最简单的Matrix问题诊断和修复脚本
## 不使用复杂的版本解析，直接尝试安装
################################################################################

cat("=====================================\n")
cat("Matrix包问题快速修复\n")
cat("=====================================\n\n")

cat("步骤 1/5: 显示当前信息\n")
cat("--------------------------------------\n")
cat("R版本:", R.version.string, "\n")

if (requireNamespace("Matrix", quietly = TRUE)) {
  cat("Matrix当前版本:", as.character(packageVersion("Matrix")), "\n")
} else {
  cat("Matrix: 未安装\n")
}
cat("\n")

cat("步骤 2/5: 卸载旧版本\n")
cat("--------------------------------------\n")
try({
  detach("package:Seurat", unload = TRUE)
  detach("package:SeuratObject", unload = TRUE)
  detach("package:Matrix", unload = TRUE)
}, silent = TRUE)

try(remove.packages("Matrix"), silent = TRUE)
cat("已卸载Matrix\n\n")

cat("步骤 3/5: 安装最新版本\n")
cat("--------------------------------------\n")
cat("正在从CRAN安装...\n")
install.packages("Matrix", repos = "https://cloud.r-project.org/", type = "source")
cat("\n")

cat("步骤 4/5: 验证安装\n")
cat("--------------------------------------\n")
library(Matrix)
new_version <- as.character(packageVersion("Matrix"))
cat("新安装的Matrix版本:", new_version, "\n\n")

cat("步骤 5/5: 测试Seurat\n")
cat("--------------------------------------\n")
test_result <- tryCatch({
  library(Seurat)
  cat("✓ Seurat加载成功!\n")
  cat("  Seurat版本:", as.character(packageVersion("Seurat")), "\n")
  TRUE
}, error = function(e) {
  cat("✗ Seurat加载失败\n")
  cat("  错误:", e$message, "\n")
  FALSE
})

cat("\n=====================================\n")
if (test_result) {
  cat("✓✓✓ 修复成功! ✓✓✓\n")
  cat("=====================================\n\n")
  cat("可以开始使用CITEViz了:\n")
  cat("source('10.script_CITEViz_ADT_analysis.R')\n\n")
} else {
  cat("修复未完成\n")
  cat("=====================================\n\n")
  cat("请尝试以下方法:\n\n")
  cat("方法1: 从源代码安装（可能需要编译工具）\n")
  cat('install.packages("Matrix", type = "source")\n\n')

  cat("方法2: 尝试binary安装\n")
  cat('install.packages("Matrix", type = "binary")\n\n')

  cat("方法3: 尝试both类型\n")
  cat('install.packages("Matrix", type = "both")\n\n')

  cat("方法4: 升级R到最新版本（最可靠）\n")
  cat("下载: https://cran.r-project.org/\n\n")

  cat("方法5: 使用旧版本Seurat（临时方案）\n")
  cat("如果Matrix无法升级，可以安装兼容的旧版本Seurat\n\n")
}
