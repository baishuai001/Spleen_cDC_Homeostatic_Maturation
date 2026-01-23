################################################################################
## 修复R包依赖问题
## 解决Matrix版本过低导致的Seurat加载失败
################################################################################

cat("=================================================================\n")
cat("正在修复R包依赖问题...\n")
cat("=================================================================\n\n")

################################################################################
## 步骤1: 更新Matrix包
################################################################################

cat("步骤1: 检查当前Matrix包版本...\n")
if (requireNamespace("Matrix", quietly = TRUE)) {
  current_version <- packageVersion("Matrix")
  cat("当前Matrix版本:", as.character(current_version), "\n")
  cat("需要的版本: >= 1.6.3\n\n")
}

cat("正在更新Matrix包到最新版本...\n")
# 先卸载旧版本
try(detach("package:Matrix", unload = TRUE), silent = TRUE)

# 安装最新版本的Matrix
install.packages("Matrix", dependencies = TRUE, repos = "https://cloud.r-project.org/")

cat("Matrix包更新完成!\n")
cat("新版本:", as.character(packageVersion("Matrix")), "\n\n")

################################################################################
## 步骤2: 重新安装Seurat相关包
################################################################################

cat("步骤2: 重新安装Seurat相关包...\n")

# 卸载可能有问题的包
packages_to_reinstall <- c("SeuratObject", "Seurat")

for (pkg in packages_to_reinstall) {
  if (pkg %in% rownames(installed.packages())) {
    cat("卸载", pkg, "...\n")
    try(remove.packages(pkg), silent = TRUE)
  }
}

# 重新安装SeuratObject
cat("\n安装SeuratObject...\n")
install.packages("SeuratObject", dependencies = TRUE, repos = "https://cloud.r-project.org/")

# 重新安装Seurat
cat("\n安装Seurat...\n")
install.packages("Seurat", dependencies = TRUE, repos = "https://cloud.r-project.org/")

################################################################################
## 步骤3: 验证安装
################################################################################

cat("\n=================================================================\n")
cat("步骤3: 验证安装...\n")
cat("=================================================================\n\n")

# 测试Matrix
cat("测试Matrix包...\n")
library(Matrix)
cat("✓ Matrix版本:", as.character(packageVersion("Matrix")), "\n\n")

# 测试SeuratObject
cat("测试SeuratObject包...\n")
library(SeuratObject)
cat("✓ SeuratObject版本:", as.character(packageVersion("SeuratObject")), "\n\n")

# 测试Seurat
cat("测试Seurat包...\n")
library(Seurat)
cat("✓ Seurat版本:", as.character(packageVersion("Seurat")), "\n\n")

################################################################################
## 步骤4: 检查CITEViz相关依赖
################################################################################

cat("=================================================================\n")
cat("步骤4: 检查CITEViz相关依赖...\n")
cat("=================================================================\n\n")

# CITEViz所需的其他包
required_packages <- c("dplyr", "ggplot2", "shiny", "plotly", "DT")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("安装", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
  } else {
    cat("✓", pkg, "已安装 (版本:", as.character(packageVersion(pkg)), ")\n")
  }
}

################################################################################
## 步骤5: 安装CITEViz
################################################################################

cat("\n=================================================================\n")
cat("步骤5: 安装CITEViz...\n")
cat("=================================================================\n\n")

if (!requireNamespace("devtools", quietly = TRUE)) {
  cat("安装devtools...\n")
  install.packages("devtools", dependencies = TRUE, repos = "https://cloud.r-project.org/")
}

library(devtools)

# 检查CITEViz是否已安装
if (!requireNamespace("CITEViz", quietly = TRUE)) {
  cat("安装CITEViz from GitHub...\n")
  devtools::install_github("maxsonBraunLab/CITE-Viz", dependencies = TRUE)
  cat("✓ CITEViz安装完成!\n")
} else {
  cat("✓ CITEViz已安装 (版本:", as.character(packageVersion("CITEViz")), ")\n")
  cat("\n如需更新CITEViz，运行:\n")
  cat("devtools::install_github('maxsonBraunLab/CITE-Viz', force = TRUE)\n")
}

################################################################################
## 完成
################################################################################

cat("\n=================================================================\n")
cat("✓ 所有依赖包已成功安装和验证!\n")
cat("=================================================================\n\n")

cat("已安装的关键包版本:\n")
cat("- R版本:", R.version.string, "\n")
cat("- Matrix:", as.character(packageVersion("Matrix")), "\n")
cat("- Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("- SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
if (requireNamespace("CITEViz", quietly = TRUE)) {
  cat("- CITEViz:", as.character(packageVersion("CITEViz")), "\n")
}

cat("\n你现在可以运行CITEViz分析脚本了:\n")
cat("source('10.script_CITEViz_ADT_analysis.R')\n\n")

cat("如果仍然遇到问题，请尝试:\n")
cat("1. 重启R会话: .rs.restartR() (RStudio) 或 q() 然后重新启动R\n")
cat("2. 更新所有包: update.packages(ask = FALSE)\n")
cat("3. 检查R版本是否 >= 4.0.0\n\n")
