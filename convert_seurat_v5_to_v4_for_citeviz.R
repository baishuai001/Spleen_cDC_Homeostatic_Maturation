#!/usr/bin/env Rscript
# ============================================================================
# Seurat v5 → v4 格式转换脚本（专为 CITEViz 设计）
# ============================================================================
# 目的：将 Seurat v5 对象（Assay5 结构）强制转换为 v3/v4 格式（Assay 结构）
# 原因：CITEViz 是基于旧版 Seurat 开发的，无法识别新的 Assay5 对象
#
# 错误症状：
# - Error: [object Object]
# - Error: 缺少参数"..4"
# - ggplot2 内部崩溃
#
# 解决方案：使用 Seurat::as.Seurat() 进行深度格式转换
# ============================================================================

library(Seurat)
library(SeuratObject)

cat("\n", strrep("=", 80), "\n")
cat("Seurat v5 → v4 格式转换（CITEViz 兼容）\n")
cat(strrep("=", 80), "\n\n")

# ============================================================================
# 1. 配置参数
# ============================================================================
INPUT_FILE <- "results/integrated/Robjects/seurat_obj_annotated.rds"
OUTPUT_FILE <- "results/integrated/Robjects/seurat_obj_citeviz_v4_compatible.rds"

# 是否降采样（推荐：用于大数据集，加快 CITEViz 加载速度）
DOWNSAMPLE <- TRUE
MAX_CELLS <- 10000  # 最多保留的细胞数

# ============================================================================
# 2. 检查 Seurat 版本
# ============================================================================
seurat_version <- packageVersion("Seurat")
cat("当前 Seurat 版本:", as.character(seurat_version), "\n")

if (seurat_version < "5.0.0") {
  cat("  ℹ 您的 Seurat 是 v4，无需转换\n")
  cat("  直接复制文件...\n")
  file.copy(INPUT_FILE, OUTPUT_FILE, overwrite = TRUE)
  quit(save = "no", status = 0)
}

cat("  ✓ 检测到 Seurat v5，需要转换\n\n")

# ============================================================================
# 3. 加载数据
# ============================================================================
cat("加载 Seurat 对象...\n")
cat("  文件:", INPUT_FILE, "\n")

seurat_obj <- readRDS(INPUT_FILE)

cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数:", nrow(seurat_obj), "\n")
cat("  检测到的 Assay:", paste(names(seurat_obj@assays), collapse=", "), "\n\n")

# ============================================================================
# 4. 检查当前 Assay 类型
# ============================================================================
cat("检查 Assay 类型...\n")
assay_types <- sapply(names(seurat_obj@assays), function(assay_name) {
  class(seurat_obj@assays[[assay_name]])[1]
})

for (assay_name in names(assay_types)) {
  cat("  ", assay_name, ":", assay_types[assay_name], "\n")
}

if (all(assay_types != "Assay5")) {
  cat("\n  ℹ 所有 Assay 已经是 v3/v4 格式，无需转换\n")
} else {
  cat("\n  ⚠ 检测到 Assay5 对象，开始转换...\n\n")
}

# ============================================================================
# 5. 降采样（可选，但强烈推荐）
# ============================================================================
if (DOWNSAMPLE && ncol(seurat_obj) > MAX_CELLS) {
  cat("降采样细胞...\n")
  cat("  原始细胞数:", ncol(seurat_obj), "\n")
  cat("  目标细胞数:", MAX_CELLS, "\n")

  set.seed(42)
  sampled_cells <- sample(colnames(seurat_obj), MAX_CELLS)
  seurat_obj <- seurat_obj[, sampled_cells]

  cat("  ✓ 降采样完成，当前细胞数:", ncol(seurat_obj), "\n\n")
}

# ============================================================================
# 6. 核心转换：Assay5 → Assay（v3/v4 格式）
# ============================================================================
cat("执行格式转换...\n")

# 方法 1：逐个 Assay 转换（最彻底）
for (assay_name in names(seurat_obj@assays)) {
  current_assay <- seurat_obj@assays[[assay_name]]

  if (inherits(current_assay, "Assay5")) {
    cat("  转换", assay_name, "（Assay5 → Assay）...\n")

    # 使用 as() 进行类型转换
    tryCatch({
      # 关键：将 Assay5 对象转换为 Assay 对象
      seurat_obj@assays[[assay_name]] <- as(current_assay, "Assay")
      cat("    ✓ 转换成功\n")
    }, error = function(e) {
      cat("    ✗ 转换失败:", e$message, "\n")
      cat("    尝试备用方法...\n")

      # 备用方法：手动重建 Assay 对象
      tryCatch({
        # 提取数据层
        counts <- LayerData(current_assay, layer = "counts")
        data <- LayerData(current_assay, layer = "data")

        # 创建新的 v3/v4 Assay 对象
        new_assay <- CreateAssayObject(counts = counts)
        new_assay <- SetAssayData(new_assay, slot = "data", new.data = data)

        # 如果有 scale.data，也复制过来
        if ("scale.data" %in% names(current_assay@layers)) {
          scale_data <- LayerData(current_assay, layer = "scale.data")
          new_assay <- SetAssayData(new_assay, slot = "scale.data", new.data = scale_data)
        }

        seurat_obj@assays[[assay_name]] <- new_assay
        cat("    ✓ 备用方法成功\n")
      }, error = function(e2) {
        cat("    ✗✗ 备用方法也失败:", e2$message, "\n")
        cat("    保留原始 Assay（可能导致 CITEViz 无法使用）\n")
      })
    })
  } else {
    cat("  ", assay_name, "已经是 v3/v4 格式，跳过\n")
  }
}

cat("\n")

# ============================================================================
# 7. 验证转换结果
# ============================================================================
cat("验证转换结果...\n")
assay_types_after <- sapply(names(seurat_obj@assays), function(assay_name) {
  class(seurat_obj@assays[[assay_name]])[1]
})

conversion_success <- TRUE
for (assay_name in names(assay_types_after)) {
  status <- if (assay_types_after[assay_name] == "Assay5") "✗ 仍是 Assay5" else "✓ 已转换为 Assay"
  cat("  ", assay_name, ":", assay_types_after[assay_name], "-", status, "\n")

  if (assay_types_after[assay_name] == "Assay5") {
    conversion_success <- FALSE
  }
}

cat("\n")

# ============================================================================
# 8. 额外的兼容性处理
# ============================================================================
cat("执行额外兼容性处理...\n")

# 确保 ADT 数据可访问
if ("ADT" %in% names(seurat_obj@assays)) {
  cat("  检查 ADT 数据...\n")

  tryCatch({
    # 尝试访问 ADT 数据（使用旧 API）
    adt_data <- GetAssayData(seurat_obj, assay = "ADT", slot = "data")
    cat("    ✓ ADT 数据可访问，行数:", nrow(adt_data), "\n")
  }, error = function(e) {
    cat("    ⚠ ADT 数据访问失败:", e$message, "\n")
  })
}

# 确保降维结果存在
if (length(seurat_obj@reductions) > 0) {
  cat("  检查降维结果:", paste(names(seurat_obj@reductions), collapse=", "), "\n")
} else {
  cat("  ⚠ 未找到降维结果，CITEViz 可能无法绘图\n")
}

cat("\n")

# ============================================================================
# 9. 保存转换后的对象
# ============================================================================
cat("保存转换后的对象...\n")
cat("  输出文件:", OUTPUT_FILE, "\n")

# 创建输出目录
dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)

# 保存
saveRDS(seurat_obj, OUTPUT_FILE)

cat("  ✓ 保存成功\n\n")

# ============================================================================
# 10. 总结
# ============================================================================
cat(strrep("=", 80), "\n")
cat("转换完成！\n")
cat(strrep("=", 80), "\n\n")

if (conversion_success) {
  cat("✓ 所有 Assay 已成功转换为 v3/v4 格式\n")
  cat("✓ 该文件应该可以被 CITEViz 正常读取\n\n")
  cat("下一步：运行以下命令启动 CITEViz\n")
  cat("  Rscript launch_citeviz.R\n\n")
  cat("或者在 R 中运行：\n")
  cat("  library(CITEViz)\n")
  cat("  seurat_obj <- readRDS('", OUTPUT_FILE, "')\n", sep = "")
  cat("  run_app(seurat_obj)  # 或者 launchApp(seurat_obj)\n\n")
} else {
  cat("⚠ 部分 Assay 未能成功转换\n")
  cat("⚠ CITEViz 可能仍然无法正常工作\n\n")
  cat("建议：\n")
  cat("  1. 检查 CITEViz 是否支持 Seurat v5\n")
  cat("  2. 考虑使用 adt_analysis_seurat.R 作为替代方案\n")
  cat("  3. 降级到 Seurat v4（不推荐）\n\n")
}

cat("转换文件信息：\n")
cat("  文件大小:", file.size(OUTPUT_FILE) / 1024^2, "MB\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  Assay:", paste(names(seurat_obj@assays), collapse=", "), "\n")
