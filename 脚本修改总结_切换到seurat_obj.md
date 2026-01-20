# 轨迹分析脚本修改总结 - 切换到 seurat_obj_annotated.rds

## 📋 修改概述

根据用户要求，将两个轨迹分析脚本从使用专用文件（`data_for_monocle3.rds`）切换到使用**通用 Seurat 对象**（`seurat_obj_annotated.rds`）。

---

## ✅ 已完成的修改

### 1. 数据加载逻辑（5.monocle3_trajectory_analysis.Rmd）

#### 修改前：
```r
# 只搜索 data_for_monocle3.rds
path <- file.path("data", mode, "Robjects", "data_for_monocle3.rds")
monocle3_data <- readRDS(input_file)
```

#### 修改后：
```r
# 优先搜索 seurat_obj_annotated.rds（新方案）⭐
path <- file.path("data", mode, "Robjects", "seurat_obj_annotated.rds")

# 备用：data_for_monocle3.rds（向后兼容）
if (!found) {
  path <- file.path("data", mode, "Robjects", "data_for_monocle3.rds")
}

# 智能加载
if (grepl("seurat_obj_annotated", file_basename)) {
  seurat_obj <- readRDS(input_file)  # 新方案
} else {
  monocle3_data <- readRDS(input_file)  # 旧方案
}
```

### 2. CDS 对象构建

#### 修改前：
```r
# 手动从 monocle3_data 构建
expression_matrix <- monocle3_data$normalized
cds <- new_cell_data_set(expression_matrix, cell_metadata, gene_metadata)
```

#### 修改后：
```r
# 新方案：使用 SeuratWrappers 转换 ⭐
if (exists("seurat_obj")) {
  library(SeuratWrappers)
  cds <- as.cell_data_set(seurat_obj)
} else {
  # 旧方案：手动构建（向后兼容）
  cds <- new_cell_data_set(...)
}
```

### 3. UMAP 命名问题修复 🔧

#### 问题：
Monocle3 的 `cluster_cells()` 和 `plot_cells()` 函数的 `reduction_method` 参数只接受固定值：
- `"UMAP"`
- `"tSNE"`
- `"PCA"`
- `"LSI"`
- `"Aligned"`

不接受自定义名称如 `"WNN_UMAP"`、`"ADT_UMAP"` 等。

#### 解决方案：
```r
# 1. 检测可用的 UMAP 空间
has_wnn_umap <- "wnn.umap" %in% names(reducedDims(cds))
has_adt_umap <- "adt.umap" %in% names(reducedDims(cds))
has_rna_umap <- "umap" %in% names(reducedDims(cds))

# 2. 根据用户配置选择UMAP（参数：preferred_umap_for_trajectory）
if (preferred_umap_for_trajectory == "wnn" && has_wnn_umap) {
  selected_umap_source <- "wnn.umap"
  selected_umap_label <- "WNN"
}

# 3. 关键：将选择的 UMAP 复制为标准名称 "UMAP" ⭐
reducedDims(cds)[["UMAP"]] <- reducedDims(cds)[[selected_umap_source]]

# 4. 同时保留原始 UMAP（用于对比可视化）
reducedDims(cds)[["RNA_UMAP"]] <- ...
reducedDims(cds)[["ADT_UMAP"]] <- ...
reducedDims(cds)[["WNN_UMAP"]] <- ...

# 5. 在所有 Monocle3 函数中使用标准名称
cds <- cluster_cells(cds, reduction_method = "UMAP")  # ✅ 正确
# NOT: reduction_method = "WNN_UMAP"  # ❌ 会报错
```

### 4. 聚类函数修复

#### 修改前：
```r
cds <- cluster_cells(
  cds,
  reduction_method = selected_umap_name,  # 可能是 "WNN_UMAP" - 会报错❌
  ...
)
```

#### 修改后：
```r
cds <- cluster_cells(
  cds,
  reduction_method = "UMAP",  # 固定使用标准名称 ✅
  ...
)
```

---

## ⚠️ 剩余问题和待修复部分

### 1. 可视化部分（5.monocle3_trajectory_analysis.Rmd）

**问题**：脚本中有 ADT UMAP、WNN UMAP、三空间对比图的代码块，它们使用了 `plot_cells()` 函数，该函数也要求 `reduction_method` 是标准名称。

**当前状态**：
- 主要轨迹可视化（基于选择的 UMAP）：✅ 已修复
- ADT UMAP 轨迹图：❌ 会报错（第874-920行）
- WNN UMAP 轨迹图：❌ 会报错（第922-968行）
- 三空间对比图：❌ 会报错（第970-1000行）

**解决方案**：

#### 方案A：删除这些代码块（最简单）
只保留主要的轨迹图（基于用户选择的 UMAP），删除其他空间的单独可视化。

#### 方案B：使用自定义绘图函数（保留功能）
不使用 `plot_cells()`，而是使用自定义的 `plot_trajectory_on_umap()` 函数（该函数不调用 Monocle3 API，直接用 ggplot2 绘图）。

```r
# 示例：使用自定义函数绘制多个 UMAP
plots_list <- list()

if ("RNA_UMAP" %in% names(reducedDims(cds))) {
  plots_list[[1]] <- plot_trajectory_on_umap(cds, "RNA_UMAP", "pseudotime")
}
if ("ADT_UMAP" %in% names(reducedDims(cds))) {
  plots_list[[2]] <- plot_trajectory_on_umap(cds, "ADT_UMAP", "pseudotime")
}
if ("WNN_UMAP" %in% names(reducedDims(cds))) {
  plots_list[[3]] <- plot_trajectory_on_umap(cds, "WNN_UMAP", "pseudotime")
}

combined <- wrap_plots(plots_list, ncol = 3)
```

**注意**：自定义函数绘制的图**不会显示轨迹曲线**（principal graph），只显示伪时间分布。如果需要轨迹曲线，必须使用 Monocle3 的 `plot_cells()`，而这要求 `reduction_method = "UMAP"`（标准名称）。

#### 方案C：为每个 UMAP 创建临时 CDS（最完整，但复杂）
为每个 UMAP 空间创建一个临时的 CDS 对象，将该 UMAP 重命名为 "UMAP"，然后调用 `plot_cells()`：

```r
# 绘制 WNN UMAP 轨迹（带曲线）
if ("WNN_UMAP" %in% names(reducedDims(cds))) {
  # 创建临时 CDS
  cds_temp <- cds
  reducedDims(cds_temp)[["UMAP"]] <- reducedDims(cds)[["WNN_UMAP"]]

  # 使用标准名称绘图
  p <- plot_cells(cds_temp, reduction_method = "UMAP", ...)

  # 保存
  ggsave("trajectory_wnn_umap.png", p)
}
```

### 2. 整合脚本（5.trajectory_analysis_integrated.Rmd）

**状态**：尚未修改

**需要修改**：
- 数据加载逻辑（切换到 `seurat_obj_annotated.rds`）
- Monocle3 CDS 构建
- Slingshot 输入准备（从 Seurat 对象提取）
- UMAP 命名问题（与单一脚本相同）

---

## 🎯 推荐行动方案

### 立即行动（用户选择）

**选项1：简化版本（推荐，快速可用）⭐**
- 删除会报错的 ADT/WNN/三空间对比可视化代码块
- 只保留基于用户选择的主要 UMAP 轨迹图
- **优势**：立即可用，无报错
- **劣势**：失去多空间对比功能

**选项2：完整版本（耗时，功能完整）**
- 使用方案B或C修复所有可视化代码块
- 保留多空间对比功能
- **优势**：功能完整，可对比三种 UMAP
- **劣势**：需要更多时间修改和测试

### 对于整合脚本

**建议**：先修复单一脚本（5.monocle3），验证可用后，再同样修改整合脚本（5.integrated）。

---

## 📝 使用说明（基于当前修改）

### 确认路径正确

输入文件应该在：
```
data/cDC1/Robjects/seurat_obj_annotated.rds  ← 新方案（优先）
data/cDC2/Robjects/seurat_obj_annotated.rds
data/integrated/Robjects/seurat_obj_annotated.rds

# 或（向后兼容）
data/cDC1/Robjects/data_for_monocle3.rds  ← 旧方案（备用）
```

### 运行脚本

```r
# 1. 打开 5.monocle3_trajectory_analysis.Rmd

# 2. 设置参数
default_analysis_mode <- "cDC1"
preferred_umap_for_trajectory <- "wnn"  # 可选: "wnn", "adt", "rna"

# 3. Knit脚本
# 脚本会自动：
#   - 搜索 seurat_obj_annotated.rds
#   - 使用 SeuratWrappers 转换
#   - 将选择的 UMAP 复制为标准名称 "UMAP"
#   - 运行轨迹分析
```

### 预期行为

```
搜索 Seurat 对象...
  检查: data/cDC1/Robjects/seurat_obj_annotated.rds ... ✓ 找到

✓ 输入文件: data/cDC1/Robjects/seurat_obj_annotated.rds
  检测到的模式: cDC1

使用 SeuratWrappers 转换 Seurat 对象...
✓ Seurat 对象已转换为 CellDataSet
  Cells: 5000
  Genes: 20000

整理 UMAP 降维空间...
检测到的 UMAP 空间:
  ✓ RNA UMAP
  ✓ ADT UMAP
  ✓ WNN UMAP

选择用于轨迹构建: WNN UMAP
  → 已将 wnn.umap 复制为标准 'UMAP' 名称

聚类细胞（基于选择的 UMAP空间： WNN UMAP）
✓ 聚类完成
```

---

## 🔧 技术要点

### SeuratWrappers 转换

```r
# SeuratWrappers 的 as.cell_data_set() 会自动：
# 1. 转换表达矩阵
# 2. 转换元数据
# 3. 继承降维结果（使用 Seurat 的命名，如 "umap"、"wnn.umap"）
# 4. 保留聚类信息

library(SeuratWrappers)
cds <- as.cell_data_set(seurat_obj)

# 检查降维名称
names(reducedDims(cds))
# [1] "PCA"      "umap"     "adt.umap" "wnn.umap"
```

### UMAP 命名映射

| Seurat | Monocle3 (after conversion) | 标准化后 |
|--------|----------------------------|---------|
| `umap` | `umap` | → `RNA_UMAP` |
| `adt_umap` | `adt.umap` | → `ADT_UMAP` |
| `wnn_umap` | `wnn.umap` | → `WNN_UMAP` |
| 用户选择 | 复制 | → `UMAP` (标准名称) ⭐ |

### 关键代码段

```r
# 检测和重命名 UMAP
if (preferred_umap_for_trajectory == "wnn" && has_wnn_umap) {
  selected_umap_source <- if("wnn.umap" %in% names(reducedDims(cds))) "wnn.umap" else "wnn_umap"
  selected_umap_label <- "WNN"
}

# 复制为标准名称（Monocle3 要求）
reducedDims(cds)[["UMAP"]] <- reducedDims(cds)[[selected_umap_source]]

# 保留原始（用于对比）
reducedDims(cds)[["WNN_UMAP"]] <- reducedDims(cds)[[selected_umap_source]]
```

---

## ✅ 测试清单

修改完成后，请测试：

- [ ] 脚本能找到 `seurat_obj_annotated.rds`
- [ ] SeuratWrappers 成功转换
- [ ] UMAP 空间选择正确（根据 `preferred_umap_for_trajectory`）
- [ ] `cluster_cells()` 不报错
- [ ] `learn_graph()` 正常运行
- [ ] 主要轨迹图生成成功
- [ ] （可选）多空间对比图正常

---

## 📚 相关文档

- **Seurat対Monocle3**: https://satijalab.org/seurat/articles/conversion_vignette.html
- **SeuratWrappers**: https://github.com/satijalab/seurat-wrappers
- **Monocle3 文档**: https://cole-trapnell-lab.github.io/monocle3/

---

**创建时间**：2026-01-20
**修改脚本**：
- 5.monocle3_trajectory_analysis.Rmd（部分完成）
- 5.trajectory_analysis_integrated.Rmd（待修改）

**下一步**：用户选择修复方案（简化版 vs 完整版）
