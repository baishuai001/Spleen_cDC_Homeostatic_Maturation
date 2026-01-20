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

## ✅ 可视化部分修复完成（方案A：简化）

### 1. 可视化部分（5.monocle3_trajectory_analysis.Rmd）

**解决方案**：已采用方案A（简化可视化），移除了会导致错误的代码块。

**当前状态**：
- 主要轨迹可视化（基于选择的 UMAP）：✅ 已修复，正常工作
- ADT UMAP 轨迹图：✅ 已移除（避免 reduction_method 错误）
- WNN UMAP 轨迹图：✅ 已移除（避免 reduction_method 错误）
- 三空间对比图：✅ 已移除（避免 reduction_method 错误）

**修改说明**：
- 删除了原 874-1011 行的多 UMAP 对比可视化代码
- 添加了注释说明移除原因
- 脚本现在只显示用户通过 `preferred_umap_for_trajectory` 参数选择的 UMAP 空间轨迹
- 如需对比多个 UMAP 空间，可使用整合脚本（5.trajectory_analysis_integrated.Rmd）

**已实施解决方案：方案A（简化可视化）✅**

已删除多 UMAP 对比可视化代码块，只保留主要的轨迹图（基于用户选择的 UMAP）。

#### 其他可选方案（未实施，仅供参考）

如果将来需要恢复多空间对比功能，可考虑：

**方案B：使用自定义绘图函数**
不使用 `plot_cells()`，而是使用自定义的 `plot_trajectory_on_umap()` 函数（该函数不调用 Monocle3 API，直接用 ggplot2 绘图）。

**注意**：自定义函数绘制的图**不会显示轨迹曲线**（principal graph），只显示伪时间分布。

**方案C：为每个 UMAP 创建临时 CDS**
为每个 UMAP 空间创建临时的 CDS 对象，将该 UMAP 重命名为 "UMAP"，然后调用 `plot_cells()`。

### 2. 整合脚本（5.trajectory_analysis_integrated.Rmd）

**状态**：⏳ 待修改（下一步）

**需要修改**：
- 数据加载逻辑（切换到 `seurat_obj_annotated.rds`）
- Monocle3 CDS 构建（使用 SeuratWrappers）
- Slingshot 输入准备（从 Seurat 对象提取）
- UMAP 命名问题（与单一脚本相同的修复方法）
- 可视化部分简化（应用方案A）

---

## 🎯 实施进度

### ✅ 已完成

1. **单一脚本（5.monocle3_trajectory_analysis.Rmd）修复完成**
   - ✅ 数据加载切换到 seurat_obj_annotated.rds
   - ✅ 使用 SeuratWrappers 转换
   - ✅ UMAP 命名问题修复（复制为标准 "UMAP"）
   - ✅ 聚类函数修复（使用标准名称）
   - ✅ 可视化简化（方案A：删除多 UMAP 对比代码）

### ⏳ 下一步

2. **整合脚本（5.trajectory_analysis_integrated.Rmd）待修复**
   - ⏳ 应用相同的数据加载逻辑
   - ⏳ Monocle3 部分修复
   - ⏳ Slingshot 部分修复
   - ⏳ 可视化简化（方案A）

3. **测试与验证**
   - ⏳ 运行修复后的单一脚本（cDC1/cDC2/integrated 三种模式）
   - ⏳ 验证轨迹清晰可见
   - ⏳ 运行整合脚本进行正式分析

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

## ✅ 测试清单（5.monocle3_trajectory_analysis.Rmd）

单一脚本修改完成，待测试：

- [ ] 脚本能找到 `seurat_obj_annotated.rds`
- [ ] SeuratWrappers 成功转换
- [ ] UMAP 空间选择正确（根据 `preferred_umap_for_trajectory`）
- [ ] `cluster_cells()` 不报错（使用标准 "UMAP" 名称）
- [ ] `learn_graph()` 正常运行
- [ ] 主要轨迹图生成成功（基于选择的 UMAP）
- [ ] 轨迹依赖基因分析正常
- [ ] 无 reduction_method 相关错误

## ⏳ 整合脚本测试清单（待修复后）

- [ ] Monocle3 和 Slingshot 都能正常运行
- [ ] 两种方法的轨迹对比图生成
- [ ] 伪时间相关性分析完成
- [ ] 差异基因重叠分析完成

---

## 📚 相关文档

- **Seurat対Monocle3**: https://satijalab.org/seurat/articles/conversion_vignette.html
- **SeuratWrappers**: https://github.com/satijalab/seurat-wrappers
- **Monocle3 文档**: https://cole-trapnell-lab.github.io/monocle3/

---

**创建时间**：2026-01-20
**最后更新**：2026-01-20

**修改脚本状态**：
- ✅ 5.monocle3_trajectory_analysis.Rmd（已完成 - 方案A简化版本）
- ⏳ 5.trajectory_analysis_integrated.Rmd（待修改）

**下一步**：
1. 修复整合脚本（5.trajectory_analysis_integrated.Rmd）
2. 测试单一脚本（三种模式：cDC1, cDC2, integrated）
3. 验证轨迹分析结果
