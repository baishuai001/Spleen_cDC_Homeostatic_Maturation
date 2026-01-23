################################################################################
## 快速修复Matrix版本问题
## 这是一个精简版的修复脚本，只解决Matrix版本问题
################################################################################

cat("正在修复Matrix版本问题...\n\n")

# 1. 卸载当前Matrix包
cat("1. 卸载旧版本Matrix...\n")
try(detach("package:Matrix", unload = TRUE), silent = TRUE)
try(remove.packages("Matrix"), silent = TRUE)

# 2. 安装最新版本Matrix
cat("2. 安装最新版本Matrix...\n")
install.packages("Matrix", dependencies = TRUE, repos = "https://cloud.r-project.org/")

# 3. 验证版本
cat("3. 验证安装...\n")
library(Matrix)
cat("✓ Matrix版本:", as.character(packageVersion("Matrix")), "\n\n")

# 4. 测试Seurat
cat("4. 测试Seurat加载...\n")
library(Seurat)
cat("✓ Seurat加载成功!\n\n")

cat("=================================================================\n")
cat("修复完成! 你现在可以正常使用Seurat了。\n")
cat("=================================================================\n")
