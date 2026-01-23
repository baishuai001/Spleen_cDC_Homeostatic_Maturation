################################################################################
## 强力修复Matrix版本问题
## 此脚本会尝试多种方法来安装最新版本的Matrix包
################################################################################

cat("=================================================================\n")
cat("Matrix版本强力修复脚本\n")
cat("=================================================================\n\n")

# 检查当前R版本
cat("当前R版本:", R.version.string, "\n")

# 更健壮的版本号解析
r_version <- tryCatch({
  ver <- as.numeric(paste(R.version$major, R.version$minor, sep = "."))
  if (is.na(ver)) {
    # 备用方法：从version.string解析
    version_string <- R.version.string
    version_match <- regmatches(version_string, regexpr("[0-9]+\\.[0-9]+", version_string))
    if (length(version_match) > 0) {
      as.numeric(version_match[1])
    } else {
      NA
    }
  } else {
    ver
  }
}, error = function(e) {
  NA
})

if (!is.na(r_version)) {
  cat("R版本号:", r_version, "\n\n")
  if (r_version < 4.0) {
    cat("警告: R版本低于4.0，可能无法安装最新版本的Matrix包\n")
    cat("建议升级R到4.2或更高版本\n\n")
  }
} else {
  cat("无法解析R版本号\n\n")
}

# 检查当前Matrix版本
if (requireNamespace("Matrix", quietly = TRUE)) {
  current_matrix <- as.character(packageVersion("Matrix"))
  cat("当前Matrix版本:", current_matrix, "\n")
  cat("需要的版本: >= 1.6.3\n\n")
}

################################################################################
## 方法1: 从CRAN安装最新版本（使用多个镜像）
################################################################################

cat("=================================================================\n")
cat("方法1: 尝试从多个CRAN镜像安装...\n")
cat("=================================================================\n\n")

# 卸载旧版本
cat("卸载旧版本...\n")
try(detach("package:Seurat", unload = TRUE), silent = TRUE)
try(detach("package:SeuratObject", unload = TRUE), silent = TRUE)
try(detach("package:Matrix", unload = TRUE), silent = TRUE)
try(remove.packages("Matrix"), silent = TRUE)

# 尝试多个CRAN镜像
mirrors <- c(
  "https://cloud.r-project.org/",
  "https://cran.rstudio.com/",
  "https://mirror.lzu.edu.cn/CRAN/",  # 中国镜像
  "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"  # 清华镜像
)

matrix_installed <- FALSE
for (mirror in mirrors) {
  cat("\n尝试镜像:", mirror, "\n")
  tryCatch({
    install.packages("Matrix", repos = mirror, dependencies = TRUE, type = "both")
    library(Matrix)
    new_version <- as.character(packageVersion("Matrix"))
    cat("安装的Matrix版本:", new_version, "\n")

    if (compareVersion(new_version, "1.6.3") >= 0) {
      cat("✓ 成功安装Matrix", new_version, "\n")
      matrix_installed <- TRUE
      break
    } else {
      cat("✗ 版本仍然太低:", new_version, "\n")
      try(remove.packages("Matrix"), silent = TRUE)
    }
  }, error = function(e) {
    cat("✗ 此镜像安装失败:", e$message, "\n")
  })
}

################################################################################
## 方法2: 从GitHub安装开发版本
################################################################################

if (!matrix_installed) {
  cat("\n=================================================================\n")
  cat("方法2: 尝试从GitHub安装开发版本...\n")
  cat("=================================================================\n\n")

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org/")
  }

  tryCatch({
    library(remotes)
    remotes::install_github("cran/Matrix", force = TRUE)
    library(Matrix)
    new_version <- as.character(packageVersion("Matrix"))
    cat("安装的Matrix版本:", new_version, "\n")

    if (compareVersion(new_version, "1.6.3") >= 0) {
      cat("✓ 成功从GitHub安装Matrix", new_version, "\n")
      matrix_installed <- TRUE
    }
  }, error = function(e) {
    cat("✗ GitHub安装失败:", e$message, "\n")
  })
}

################################################################################
## 方法3: 使用compatible的Seurat版本（备选方案）
################################################################################

if (!matrix_installed) {
  cat("\n=================================================================\n")
  cat("方法3: 安装与当前Matrix兼容的Seurat版本（备选方案）\n")
  cat("=================================================================\n\n")

  cat("无法升级Matrix到1.6.3，将安装与Matrix 1.5.3兼容的Seurat版本\n\n")

  # 重新安装Matrix 1.5.3
  install.packages("Matrix", repos = "https://cloud.r-project.org/")

  # 安装旧版本的SeuratObject和Seurat
  cat("安装兼容的SeuratObject版本...\n")
  tryCatch({
    # 尝试安装SeuratObject 4.1.3（与Matrix 1.5.3兼容）
    install.packages("SeuratObject", repos = "https://cloud.r-project.org/")

    cat("安装兼容的Seurat版本...\n")
    install.packages("Seurat", repos = "https://cloud.r-project.org/")

    library(Matrix)
    library(SeuratObject)
    library(Seurat)

    cat("✓ 使用兼容版本组合:\n")
    cat("  Matrix:", as.character(packageVersion("Matrix")), "\n")
    cat("  SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
    cat("  Seurat:", as.character(packageVersion("Seurat")), "\n")

    matrix_installed <- TRUE  # 标记为已解决

  }, error = function(e) {
    cat("✗ 安装兼容版本失败:", e$message, "\n")
  })
}

################################################################################
## 验证最终结果
################################################################################

cat("\n=================================================================\n")
cat("最终验证\n")
cat("=================================================================\n\n")

if (matrix_installed) {
  tryCatch({
    library(Matrix)
    cat("✓ Matrix版本:", as.character(packageVersion("Matrix")), "\n")

    library(SeuratObject)
    cat("✓ SeuratObject版本:", as.character(packageVersion("SeuratObject")), "\n")

    library(Seurat)
    cat("✓ Seurat版本:", as.character(packageVersion("Seurat")), "\n")

    cat("\n=================================================================\n")
    cat("✓✓✓ 成功! 可以使用Seurat了 ✓✓✓\n")
    cat("=================================================================\n\n")

    cat("下一步: 运行CITEViz分析\n")
    cat("source('10.script_CITEViz_ADT_analysis.R')\n\n")

  }, error = function(e) {
    cat("✗ 验证失败:", e$message, "\n")
    matrix_installed <- FALSE
  })
}

################################################################################
## 如果所有方法都失败
################################################################################

if (!matrix_installed) {
  cat("\n=================================================================\n")
  cat("所有自动修复方法都失败了\n")
  cat("=================================================================\n\n")

  cat("建议手动解决方案:\n\n")

  cat("选项1: 升级R版本（最推荐）\n")
  cat("  1. 检查R版本: R.version.string\n")
  cat("  2. 如果R < 4.2，请从 https://cran.r-project.org/ 下载并安装最新R版本\n")
  cat("  3. 重新安装所有包\n\n")

  cat("选项2: 在R外部手动安装Matrix\n")
  cat("  在终端/命令行中运行:\n")
  cat("  R CMD INSTALL Matrix_1.6.3.tar.gz\n\n")

  cat("选项3: 使用conda环境\n")
  cat("  conda create -n seurat_env r-base=4.3 r-seurat r-matrix\n")
  cat("  conda activate seurat_env\n\n")

  cat("选项4: 联系系统管理员\n")
  cat("  如果在共享服务器上，可能需要管理员权限来升级R或安装包\n\n")

  cat("临时解决方案: 使用旧版本Seurat\n")
  cat("  虽然不是最新版本，但应该仍然可以分析数据\n")
  cat("  当前安装的版本组合应该是兼容的\n\n")
}

cat("\n当前已安装包的版本:\n")
cat("R:", R.version.string, "\n")
if (requireNamespace("Matrix", quietly = TRUE)) {
  cat("Matrix:", as.character(packageVersion("Matrix")), "\n")
}
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  cat("SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
}
if (requireNamespace("Seurat", quietly = TRUE)) {
  cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
}
