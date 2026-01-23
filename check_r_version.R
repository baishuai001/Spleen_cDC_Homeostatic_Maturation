################################################################################
## 检查R版本和包兼容性
################################################################################

cat("=================================================================\n")
cat("R环境诊断报告\n")
cat("=================================================================\n\n")

# R版本信息
cat("【R版本信息】\n")
cat("完整版本:", R.version.string, "\n")
cat("主版本号:", R.version$major, "\n")
cat("次版本号:", R.version$minor, "\n")
r_version <- as.numeric(paste(R.version$major, R.version$minor, sep = "."))
cat("版本号:", r_version, "\n\n")

# 版本兼容性检查
cat("【版本兼容性分析】\n")
if (r_version >= 4.3) {
  cat("✓ R版本优秀 (>= 4.3) - 支持所有最新R包\n")
} else if (r_version >= 4.2) {
  cat("✓ R版本良好 (>= 4.2) - 支持Matrix 1.6.3+\n")
} else if (r_version >= 4.0) {
  cat("⚠ R版本较旧 (4.0-4.1) - 可能需要从源代码编译Matrix\n")
  cat("  建议升级到R 4.2或更高版本\n")
} else {
  cat("✗ R版本过旧 (< 4.0) - 不支持最新版本的Seurat和Matrix\n")
  cat("  强烈建议升级到R 4.3或更高版本\n")
}
cat("\n")

# 当前已安装的包版本
cat("【已安装包版本】\n")
packages <- c("Matrix", "SeuratObject", "Seurat", "dplyr", "ggplot2")
for (pkg in packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- as.character(packageVersion(pkg))
    cat(sprintf("%-15s: %s\n", pkg, version))
  } else {
    cat(sprintf("%-15s: 未安装\n", pkg))
  }
}
cat("\n")

# Matrix版本特定检查
cat("【Matrix包详细信息】\n")
if (requireNamespace("Matrix", quietly = TRUE)) {
  matrix_version <- as.character(packageVersion("Matrix"))
  cat("当前版本:", matrix_version, "\n")
  cat("需要版本: >= 1.6.3\n")

  if (compareVersion(matrix_version, "1.6.3") >= 0) {
    cat("状态: ✓ 版本满足要求\n")
  } else {
    cat("状态: ✗ 版本过低\n")
    version_diff <- compareVersion(matrix_version, "1.6.3")
    cat("版本差距:", abs(version_diff), "个版本\n")
  }
} else {
  cat("Matrix未安装\n")
}
cat("\n")

# 库路径信息
cat("【R包安装路径】\n")
lib_paths <- .libPaths()
for (i in seq_along(lib_paths)) {
  cat(sprintf("%d. %s\n", i, lib_paths[i]))
  cat(sprintf("   可写: %s\n", file.access(lib_paths[i], 2) == 0))
}
cat("\n")

# CRAN镜像设置
cat("【CRAN镜像设置】\n")
repos <- getOption("repos")
cat("当前CRAN镜像:", repos["CRAN"], "\n\n")

# 系统信息
cat("【系统信息】\n")
cat("操作系统:", Sys.info()["sysname"], "\n")
cat("架构:", Sys.info()["machine"], "\n")
cat("节点名:", Sys.info()["nodename"], "\n\n")

# 建议
cat("=================================================================\n")
cat("【诊断建议】\n")
cat("=================================================================\n\n")

if (r_version < 4.2) {
  cat("🔴 紧急建议: 升级R版本\n")
  cat("   当前R版本不支持Matrix 1.6.3+\n")
  cat("   请从 https://cran.r-project.org/ 下载R 4.3或更高版本\n\n")

  cat("升级步骤:\n")
  cat("1. 下载最新R版本: https://cran.r-project.org/\n")
  cat("2. 安装新版本R\n")
  cat("3. 重新安装所有R包\n")
  cat("4. 运行: install.packages(c('Matrix', 'Seurat'))\n\n")

} else {
  matrix_version <- as.character(packageVersion("Matrix"))
  if (compareVersion(matrix_version, "1.6.3") < 0) {
    cat("🟡 建议: 更新Matrix包\n")
    cat("   你的R版本支持Matrix 1.6.3+，但当前安装的版本过低\n\n")

    cat("更新方法:\n")
    cat("1. 退出R\n")
    cat("2. 重新启动R\n")
    cat("3. 运行: remove.packages('Matrix')\n")
    cat("4. 运行: install.packages('Matrix', type='both')\n")
    cat("5. 如果仍然失败，运行: source('force_fix_matrix.R')\n\n")
  } else {
    cat("🟢 一切正常!\n")
    cat("   你的R环境配置良好，可以运行CITEViz分析\n\n")
  }
}

cat("如需进一步帮助，请提供此诊断报告的输出\n\n")
