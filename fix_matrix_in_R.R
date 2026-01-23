################################################################################
## 在R中直接运行此代码块来修复Matrix版本问题
## 复制粘贴下面的全部代码到R控制台
################################################################################

cat("=================================================================\n")
cat("开始修复Matrix版本问题 (当前版本太低)\n")
cat("=================================================================\n\n")

# 1. 尝试卸载所有相关包
cat("步骤1: 卸载相关包...\n")
packages_to_remove <- c("Seurat", "SeuratObject", "Matrix")
for (pkg in packages_to_remove) {
  if (pkg %in% rownames(installed.packages())) {
    cat("  - 卸载", pkg, "...\n")
    try(detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE), silent = TRUE)
    try(remove.packages(pkg), silent = TRUE)
  }
}

cat("\n步骤2: 安装最新版本的Matrix...\n")
install.packages("Matrix", repos = "https://cloud.r-project.org/", dependencies = TRUE)

cat("\n步骤3: 重新安装Seurat相关包...\n")
install.packages("SeuratObject", repos = "https://cloud.r-project.org/", dependencies = TRUE)
install.packages("Seurat", repos = "https://cloud.r-project.org/", dependencies = TRUE)

cat("\n=================================================================\n")
cat("安装完成！正在验证...\n")
cat("=================================================================\n\n")

# 验证
library(Matrix)
cat("✓ Matrix 版本:", as.character(packageVersion("Matrix")), "\n")

library(SeuratObject)
cat("✓ SeuratObject 版本:", as.character(packageVersion("SeuratObject")), "\n")

library(Seurat)
cat("✓ Seurat 版本:", as.character(packageVersion("Seurat")), "\n")

cat("\n=================================================================\n")
cat("✓✓✓ 修复成功! 现在可以使用Seurat了 ✓✓✓\n")
cat("=================================================================\n\n")

cat("下一步: 运行CITEViz分析\n")
cat("source('10.script_CITEViz_ADT_analysis.R')\n\n")
