#!/usr/bin/env Rscript
# CITEViz 安装检查脚本

cat("=================================================\n")
cat("CITEViz 安装检查\n")
cat("=================================================\n\n")

# 1. 检查包是否已安装
cat("1. 检查 CITEViz 包是否已安装...\n")
if ("CITEViz" %in% installed.packages()[, "Package"]) {
  cat("  ✓ CITEViz 已安装\n")

  # 显示版本信息
  version <- packageVersion("CITEViz")
  cat("  版本:", as.character(version), "\n")

  # 显示安装位置
  path <- find.package("CITEViz")
  cat("  位置:", path, "\n")

} else {
  cat("  ✗ CITEViz 未安装\n")
  cat("\n安装方法:\n")
  cat("  devtools::install_github('maxsonBraunLab/CITEViz')\n\n")
  stop("请先安装 CITEViz")
}

# 2. 尝试加载包
cat("\n2. 尝试加载 CITEViz...\n")
tryCatch({
  library(CITEViz)
  cat("  ✓ CITEViz 加载成功\n")
}, error = function(e) {
  cat("  ✗ 加载失败:", e$message, "\n")
  stop("CITEViz 加载失败")
})

# 3. 检查包的依赖
cat("\n3. 检查依赖包...\n")
required_packages <- c("Seurat", "shiny", "ggplot2", "plotly", "DT", "shinydashboard")

for (pkg in required_packages) {
  if (pkg %in% installed.packages()[, "Package"]) {
    cat("  ✓", pkg, "\n")
  } else {
    cat("  ✗", pkg, "未安装\n")
  }
}

# 4. 检查可用的函数
cat("\n4. 检查 CITEViz 的主要函数...\n")

# 列出所有导出的函数
exported_functions <- ls("package:CITEViz")
cat("  导出的函数数量:", length(exported_functions), "\n")

cat("\n  可用函数列表:\n")
for (func in exported_functions) {
  cat("    -", func, "\n")
}

# 检查关键函数
key_functions <- c("run_app", "launchApp", "CITEViz")
cat("\n  关键函数检查:\n")
for (func in key_functions) {
  if (exists(func, where = "package:CITEViz")) {
    cat("    ✓", func, "\n")
  } else {
    cat("    ✗", func, "不存在\n")
  }
}

# 5. 检查 Seurat 版本兼容性
cat("\n5. 检查 Seurat 版本兼容性...\n")
seurat_version <- packageVersion("Seurat")
cat("  Seurat 版本:", as.character(seurat_version), "\n")

if (seurat_version >= "5.0.0") {
  cat("  ⚠ 警告：Seurat v5 可能与 CITEViz 不完全兼容\n")
  cat("  CITEViz 可能是基于 Seurat v4 开发的\n")
} else {
  cat("  ✓ Seurat v4，应该兼容\n")
}

# 6. 尝试获取帮助文档
cat("\n6. 检查帮助文档...\n")
tryCatch({
  help_file <- help("CITEViz", package = "CITEViz")
  if (!is.null(help_file)) {
    cat("  ✓ 帮助文档可用\n")
  }
}, error = function(e) {
  cat("  ⚠ 帮助文档不可用\n")
})

# 7. 检查示例数据
cat("\n7. 检查包内容...\n")
pkg_path <- find.package("CITEViz")

# 检查 R 目录
r_dir <- file.path(pkg_path, "R")
if (dir.exists(r_dir)) {
  r_files <- list.files(r_dir, pattern = "\\.R$")
  cat("  R 脚本数量:", length(r_files), "\n")
}

# 检查 inst 目录（shiny app）
inst_dir <- file.path(pkg_path, "inst")
if (dir.exists(inst_dir)) {
  cat("  ✓ inst 目录存在\n")

  # 检查 shiny app
  app_dirs <- list.dirs(inst_dir, recursive = FALSE)
  if (length(app_dirs) > 0) {
    cat("  Shiny 应用目录:\n")
    for (dir in app_dirs) {
      cat("    -", basename(dir), "\n")
    }
  }
}

# 8. 测试启动函数（不实际启动）
cat("\n8. 测试启动函数签名...\n")

# 检查主启动函数的参数
if (exists("run_app", where = "package:CITEViz")) {
  cat("  run_app 函数存在\n")

  # 获取函数参数
  tryCatch({
    args <- formals(CITEViz::run_app)
    cat("  函数参数:\n")
    for (arg_name in names(args)) {
      cat("    -", arg_name, "\n")
    }
  }, error = function(e) {
    cat("  无法获取函数参数\n")
  })
}

# 总结
cat("\n=================================================\n")
cat("检查总结\n")
cat("=================================================\n\n")

cat("安装状态: ✓ CITEViz 已安装并可加载\n")
cat("Seurat 版本:", as.character(seurat_version), "\n")

if (seurat_version >= "5.0.0") {
  cat("\n⚠ 重要提醒:\n")
  cat("  - 您使用的是 Seurat v5\n")
  cat("  - CITEViz 可能不完全支持 Seurat v5\n")
  cat("  - 这可能是之前报错的根本原因\n")
  cat("\n建议:\n")
  cat("  1. 使用修复后的 adt_analysis_seurat.R（推荐）\n")
  cat("  2. 或降级到 Seurat v4（不推荐）\n")
}

cat("\n启动方法:\n")
if (exists("run_app", where = "package:CITEViz")) {
  cat("  library(CITEViz)\n")
  cat("  run_app()  # 然后上传 RDS 文件\n")
} else if (exists("launchApp", where = "package:CITEViz")) {
  cat("  library(CITEViz)\n")
  cat("  launchApp()  # 然后上传 RDS 文件\n")
} else {
  cat("  未找到标准启动函数\n")
  cat("  请查看包文档\n")
}

cat("\n如需查看更多信息:\n")
cat("  help(package = 'CITEViz')\n")
cat("  browseVignettes('CITEViz')\n")
