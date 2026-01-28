#!/usr/bin/env Rscript
# ============================================================================
# Seurat 对象自动检测脚本
# ============================================================================
# 目的：在项目中查找可用的 Seurat 对象文件
# 用途：在运行 CITEViz 之前检查是否有可用的数据
# ============================================================================

cat("\n", strrep("=", 80), "\n")
cat("Seurat 对象文件检测\n")
cat(strrep("=", 80), "\n\n")

# ============================================================================
# 1. 定义可能的文件位置
# ============================================================================
possible_files <- c(
  # 标准输出位置
  "results/integrated/Robjects/seurat_obj_annotated.rds",
  "results/integrated/Robjects/seurat_obj_for_annotation.rds",
  "results/integrated/seurat_obj_annotated.rds",

  # CITEViz 相关
  "results/integrated/Robjects/seurat_obj_citeviz_compatible.rds",
  "results/integrated/Robjects/seurat_obj_citeviz_v4_compatible.rds",
  "results/citeviz/citeviz_fast_5k.rds",
  "results/citeviz/citeviz_medium_8k.rds",
  "results/citeviz/citeviz_full_compatible.rds",

  # 可能的备份位置
  "seurat_obj_annotated.rds",
  "seurat_obj.rds",
  "data/seurat_obj_annotated.rds"
)

# ============================================================================
# 2. 检查文件是否存在
# ============================================================================
cat("检查预定义位置...\n")

found_files <- c()
for (file_path in possible_files) {
  if (file.exists(file_path)) {
    file_size <- file.size(file_path) / 1024^2  # 转换为 MB
    cat("  ✓ 找到:", file_path, sprintf("(%.1f MB)", file_size), "\n")
    found_files <- c(found_files, file_path)
  }
}

if (length(found_files) == 0) {
  cat("  ✗ 未在预定义位置找到文件\n\n")
}

# ============================================================================
# 3. 全局搜索 .rds 文件
# ============================================================================
cat("\n执行全局搜索 .rds 文件...\n")
cat("  （这可能需要一些时间）\n")

# 使用 system2 执行 find 命令
rds_files <- tryCatch({
  find_output <- system2("find",
                         args = c(".", "-name", "*.rds", "-type", "f",
                                  "!", "-path", "*/.git/*",
                                  "!", "-path", "*/renv/*"),
                         stdout = TRUE, stderr = FALSE)
  if (length(find_output) > 0) {
    find_output
  } else {
    character(0)
  }
}, error = function(e) {
  character(0)
})

if (length(rds_files) > 0) {
  cat("\n找到以下 .rds 文件:\n")

  # 检查每个文件，尝试判断是否是 Seurat 对象
  for (rds_file in rds_files) {
    file_size <- file.size(rds_file) / 1024^2

    # 文件大小过滤（Seurat 对象通常 > 10 MB）
    if (file_size > 10) {
      cat("  ", rds_file, sprintf("(%.1f MB)", file_size), "\n")

      # 尝试快速检查是否是 Seurat 对象
      is_seurat <- tryCatch({
        obj <- readRDS(rds_file)
        inherits(obj, "Seurat")
      }, error = function(e) {
        FALSE
      })

      if (is_seurat) {
        cat("    → ✓ 确认为 Seurat 对象\n")
        if (!rds_file %in% found_files) {
          found_files <- c(found_files, rds_file)
        }
      }
    }
  }
} else {
  cat("  ✗ 未找到任何 .rds 文件\n")
}

# ============================================================================
# 4. 总结和建议
# ============================================================================
cat("\n", strrep("=", 80), "\n")
cat("检测结果总结\n")
cat(strrep("=", 80), "\n\n")

if (length(found_files) > 0) {
  cat("✓ 找到", length(found_files), "个可用的 Seurat 对象文件\n\n")

  cat("推荐使用（按优先级排序）:\n")

  # 按优先级排序
  priority_order <- c(
    "seurat_obj_citeviz_v4_compatible.rds",
    "seurat_obj_citeviz_compatible.rds",
    "seurat_obj_annotated.rds",
    "seurat_obj.rds"
  )

  for (priority_file in priority_order) {
    matching_files <- grep(priority_file, found_files, value = TRUE)
    if (length(matching_files) > 0) {
      cat("\n优先级", which(priority_order == priority_file), ":\n")
      for (f in matching_files) {
        cat("  ", f, "\n")
      }
    }
  }

  # 显示其他文件
  other_files <- found_files[!grepl(paste(priority_order, collapse="|"), found_files)]
  if (length(other_files) > 0) {
    cat("\n其他可用文件:\n")
    for (f in other_files) {
      cat("  ", f, "\n")
    }
  }

  cat("\n下一步操作:\n")
  cat("  1. 转换为 v4 兼容格式:\n")
  cat("     修改 convert_seurat_v5_to_v4_for_citeviz.R 中的 INPUT_FILE 路径\n")
  cat("     然后运行: Rscript convert_seurat_v5_to_v4_for_citeviz.R\n\n")
  cat("  2. 或者使用 Seurat 原生分析:\n")
  cat("     修改 adt_analysis_seurat.R 中的输入文件路径\n")
  cat("     然后运行: Rscript adt_analysis_seurat.R\n\n")

} else {
  cat("✗ 未找到任何 Seurat 对象文件\n\n")

  cat("这意味着您还没有运行前面的分析步骤。\n\n")

  cat("必需步骤:\n")
  cat("  1. 运行下游分析脚本生成 Seurat 对象:\n")
  cat("     在 R 中运行: rmarkdown::render('4.downstream_analysis.Rmd')\n")
  cat("     或使用 knitr: knitr::knit('4.downstream_analysis.Rmd')\n\n")

  cat("  2. 该脚本会生成以下文件:\n")
  cat("     results/integrated/Robjects/seurat_obj_annotated.rds\n\n")

  cat("  3. 生成后即可运行 CITEViz 相关工具\n\n")

  cat("快速检查:\n")
  cat("  - 查看是否有 .Rmd 文件: ls -lh *.Rmd\n")
  cat("  - 查看项目 README: cat README.md\n\n")
}

cat(strrep("=", 80), "\n")
cat("检测完成\n")
cat(strrep("=", 80), "\n")
