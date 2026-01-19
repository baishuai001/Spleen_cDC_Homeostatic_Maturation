# 下游分析脚本 WNN UMAP 完整支持总结

## 📅 修改日期
2026-01-19

## 🎯 修改目标

为 `4.downstream_analysis.Rmd` 的**所有可视化模块**添加 WNN UMAP 支持，确保用户可以在多模态整合空间查看所有分析结果。

---

## ✅ 完成的修改

### 1. 输入文件优先级修复 ⭐

**位置：** 第221-226行

**问题：** 脚本优先加载 `*_pre_wnn.rds`（WNN分析前的对象，无WNN UMAP），导致后续可视化缺少WNN UMAP。

**解决方案：** 修改文件加载优先级，优先加载 `*_final.rds`（WNN分析后的对象，包含WNN UMAP）。

```r
# 优先加载final.rds（包含WNN UMAP）
if (file.exists(input_file_final)) {
  input_file <- input_file_final
  cat("  ✅ 找到 final.rds（包含WNN UMAP）:", input_file, "\n")
} else if (file.exists(input_file_prewnn)) {
  input_file <- input_file_prewnn
  cat("  ⚠️  未找到 final.rds，使用 pre_wnn.rds（无WNN UMAP）:", input_file, "\n")
}
```

**影响：** 🔑 这是所有WNN UMAP功能的基础，确保后续所有可视化都能访问WNN UMAP。

---

### 2. WNN UMAP 作者注释可视化 ⭐

**位置：** 第1633-1697行

**功能：** 在WNN UMAP降维空间上显示作者的原始细胞类型注释。

#### 2.1 单图：WNN UMAP + 作者注释

```r
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  p_wnn_author <- DimPlot(seurat_obj,
                          reduction = "wnn_umap",
                          group.by = "author_annotation",
                          label = TRUE,
                          repel = TRUE,
                          label.size = 4)

  # 输出文件: author_annotation_wnn_umap.png/.pdf
}
```

**效果：**
- 在WNN-combined Clustering UMAP上显示真实细胞类型名称（如 Pre-cDC1s, Mature cDC1s）
- 而非仅显示聚类数字（0, 1, 2, 3...）

#### 2.2 三联对比图

```r
# 对比三种视图：
# - WNN Clusters（算法聚类）
# - Author Annotation（作者注释）
# - Detailed Annotation（详细注释）
p_wnn_clusters | p_wnn_annot | p_wnn_detailed

# 输出文件: wnn_umap_triple_comparison.png/.pdf
```

**用途：**
- 评估WNN聚类与作者注释的一致性
- 发现WNN算法分得更细或合并的情况

---

### 3. 多模态 Marker 综合可视化 - WNN UMAP 支持

**位置：** 第2028-2071行

**功能：** 为阶段标志基因的FeaturePlot添加WNN UMAP版本。

#### 修改内容

**原始代码：** 仅生成RNA UMAP版本
```r
p_feature <- FeaturePlot(seurat_obj,
                         features = key_genes,
                         reduction = "umap", ...)  # 仅RNA UMAP
```

**修改后：** 同时生成RNA UMAP和WNN UMAP版本
```r
# --- RNA UMAP版本 ---
p_feature <- FeaturePlot(seurat_obj,
                         features = key_genes,
                         reduction = "umap", ...)
plot_and_save(p_feature, "multimodal_featureplot_stage_markers", ...)

# --- WNN UMAP版本（新增）---
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  p_feature_wnn <- FeaturePlot(seurat_obj,
                               features = key_genes,
                               reduction = "wnn_umap", ...)
  plot_and_save(p_feature_wnn, "multimodal_featureplot_stage_markers_wnn", ...)
}
```

#### 输出文件

| 文件名 | 降维空间 | 描述 |
|--------|---------|------|
| `multimodal_featureplot_stage_markers.png` | RNA UMAP | 阶段标志基因（RNA空间） |
| `multimodal_featureplot_stage_markers_wnn.png` | WNN UMAP | 阶段标志基因（WNN整合空间）⭐ |

---

### 4. 扩展基因组可视化 - 基因 FeaturePlot WNN 支持

**位置：** 第2191-2251行

**功能：** 为8组功能基因的FeaturePlot添加WNN UMAP版本。

#### 8组功能基因

1. **调控基因** (7个): Cd274, Pdcd1lg2, Cd200, Fas, Aldh1a2, Socs1, Socs2
2. **迁移基因** (8个): Ccr7, Myo1g, Cxcl16, Adam8, Icam1, Fscn1, Marcks, Marcksl1
3. **标记基因** (7个): H2-Ab1, Itgax, Itgam, Itgae, Xcr1, Sirpa, Cd24a
4. **成熟基因** (5个): Cd80, Cd86, Cd40, Relb, Cd83
5. **TH2反应基因** (7个): Il4ra, Il4i1, Ccl17, Ccl22, Tnfrsf4, Stat6, Bcl2l1
6. **转录因子** (15个): Sfpi1, Zbtb46, Irf8, Nfil3, Id2, Batf3, Irf4, Cebpa, Cebpb, Zeb2, Klf4, Tbx21, Rorc, Notch2, Nr4a3
7. **发育门控** (12个): Itgax, H2-Ab1, Flt3, Kit, Sirpa, Csf1r, Siglech, Ly6c1, Itgam, Spn, Fcgr1, Ly6d

**总计：61个基因**

#### 修改内容

**原始代码：** 仅生成RNA UMAP版本
```r
for (group_name in names(genes_extended_groups)) {
  p_feature_group <- FeaturePlot(seurat_obj,
                                 features = genes_to_plot,
                                 reduction = "umap", ...)  # 仅RNA UMAP
  plot_and_save(p_feature_group, paste0("extended_genes_featureplot_", filename_safe), ...)
}
```

**修改后：** 每组基因生成RNA UMAP和WNN UMAP两个版本
```r
for (group_name in names(genes_extended_groups)) {
  # --- RNA UMAP版本 ---
  p_feature_group <- FeaturePlot(seurat_obj,
                                 features = genes_to_plot,
                                 reduction = "umap", ...)
  plot_and_save(p_feature_group, paste0("extended_genes_featureplot_", filename_safe), ...)

  # --- WNN UMAP版本（新增）---
  if ("wnn_umap" %in% names(seurat_obj@reductions)) {
    p_feature_group_wnn <- FeaturePlot(seurat_obj,
                                       features = genes_to_plot,
                                       reduction = "wnn_umap", ...)
    plot_and_save(p_feature_group_wnn,
                  paste0("extended_genes_featureplot_", filename_safe, "_wnn"), ...)
  }
}
```

#### 输出文件（每组2个文件）

**RNA UMAP版本：**
1. `extended_genes_featureplot_调控基因.png`
2. `extended_genes_featureplot_迁移基因.png`
3. `extended_genes_featureplot_标记基因.png`
4. `extended_genes_featureplot_成熟基因.png`
5. `extended_genes_featureplot_TH2反应基因.png`
6. `extended_genes_featureplot_转录因子.png`
7. `extended_genes_featureplot_发育门控.png`

**WNN UMAP版本（新增）：**
1. `extended_genes_featureplot_调控基因_wnn.png` ⭐
2. `extended_genes_featureplot_迁移基因_wnn.png` ⭐
3. `extended_genes_featureplot_标记基因_wnn.png` ⭐
4. `extended_genes_featureplot_成熟基因_wnn.png` ⭐
5. `extended_genes_featureplot_TH2反应基因_wnn.png` ⭐
6. `extended_genes_featureplot_转录因子_wnn.png` ⭐
7. `extended_genes_featureplot_发育门控_wnn.png` ⭐

---

### 5. 扩展基因组可视化 - ADT蛋白 FeaturePlot WNN 支持

**位置：** 第2276-2330行

**功能：** 为ADT发育门控蛋白的FeaturePlot添加WNN UMAP版本。

#### ADT蛋白列表（12个）

adt_CD11c, adt_IA-IE, adt_CD135, adt_CD117, adt_CD172a, adt_CD115, adt_SiglecH, adt_Ly6C, adt_CD11b, adt_CD43, adt_CD64, adt_Ly6D

#### 修改内容

**原始代码：** 仅生成RNA UMAP版本
```r
for (i in seq(1, length(adt_available), by = 9)) {
  adt_batch <- adt_available[i:min(i+8, length(adt_available))]

  p_feature_adt <- FeaturePlot(seurat_obj,
                               features = adt_batch,
                               reduction = "umap", ...)  # 仅RNA UMAP
  plot_and_save(p_feature_adt, paste0("adt_proteins_featureplot_batch", ceiling(i/9)), ...)
}
```

**修改后：** 每批次生成RNA UMAP和WNN UMAP两个版本
```r
for (i in seq(1, length(adt_available), by = 9)) {
  adt_batch <- adt_available[i:min(i+8, length(adt_available))]

  # --- RNA UMAP版本 ---
  p_feature_adt <- FeaturePlot(seurat_obj,
                               features = adt_batch,
                               reduction = "umap",
                               cols = c("lightgrey", "blue"), ...)
  plot_and_save(p_feature_adt, paste0("adt_proteins_featureplot_batch", ceiling(i/9)), ...)

  # --- WNN UMAP版本（新增）---
  if ("wnn_umap" %in% names(seurat_obj@reductions)) {
    p_feature_adt_wnn <- FeaturePlot(seurat_obj,
                                     features = adt_batch,
                                     reduction = "wnn_umap",
                                     cols = c("lightgrey", "blue"), ...)
    plot_and_save(p_feature_adt_wnn,
                  paste0("adt_proteins_featureplot_batch", ceiling(i/9), "_wnn"), ...)
  }
}
```

#### 输出文件

**RNA UMAP版本：**
- `adt_proteins_featureplot_batch1.png`（前9个蛋白）
- `adt_proteins_featureplot_batch2.png`（后3个蛋白，如果有）

**WNN UMAP版本（新增）：**
- `adt_proteins_featureplot_batch1_wnn.png` ⭐
- `adt_proteins_featureplot_batch2_wnn.png` ⭐

**颜色方案：** 蓝色渐变 (lightgrey → blue)，与基因的红色渐变形成对比。

---

### 6. 下游分析数据准备 - WNN UMAP 支持

**位置：** 第2344-2458行

**功能：** 为轨迹分析工具（Monocle3、CellOracle等）准备的数据添加WNN UMAP坐标。

#### 6.1 Monocle3 数据

**位置：** 第2344-2356行

```r
monocle3_data <- list(
  counts = ...,
  metadata = ...,
  umap = as.data.frame(Embeddings(seurat_obj, "umap")),
  wnn_umap = if ("wnn_umap" %in% names(seurat_obj@reductions)) {  # 新增
               as.data.frame(Embeddings(seurat_obj, "wnn_umap"))
             } else {
               NULL
             }
)
```

**输出文件：** `results/{mode}/Trajectory/monocle3_data.rds`

#### 6.2 CellOracle 数据

**位置：** 第2381-2393行

```r
celloracle_data <- list(
  ...,
  umap_coords = as.data.frame(Embeddings(seurat_obj, "umap")),
  wnn_umap_coords = if ("wnn_umap" %in% names(seurat_obj@reductions)) {  # 新增
                      as.data.frame(Embeddings(seurat_obj, "wnn_umap"))
                    } else {
                      NULL
                    }
)
```

**输出文件：** `results/{mode}/CellOracle/celloracle_data.rds`

#### 6.3 CSV导出

**位置：** 第2442-2458行

```r
# WNN UMAP坐标导出（新增）
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  wnn_umap_df <- as.data.frame(Embeddings(seurat_obj, "wnn_umap"))
  wnn_umap_df$cell_id <- rownames(wnn_umap_df)
  wnn_umap_df$cell_type <- seurat_obj$author_annotation[rownames(wnn_umap_df)]
  wnn_umap_df$cluster <- seurat_obj$WNN_clusters[rownames(wnn_umap_df)]

  write.csv(wnn_umap_df,
            file.path(output.dir, "Exported_Data", "wnn_umap_coordinates.csv"),
            row.names = FALSE)
}
```

**输出文件：** `results/{mode}/Exported_Data/wnn_umap_coordinates.csv`

**CSV内容：**
- `wnnUMAP_1`, `wnnUMAP_2`: WNN UMAP坐标
- `cell_id`: 细胞ID
- `cell_type`: 作者注释
- `cluster`: WNN聚类

---

## 📊 完整输出文件清单

### 1. 注释可视化（`results/{mode}/Plots/Annotation/`）

| 文件名 | 降维空间 | 描述 |
|--------|---------|------|
| `author_annotation_umap.png` | RNA UMAP | 作者注释（RNA空间） |
| `author_annotation_wnn_umap.png` ⭐ | WNN UMAP | 作者注释（WNN整合空间） |
| `wnn_umap_triple_comparison.png` ⭐ | WNN UMAP | WNN三联对比图 |

### 2. 标志基因可视化（`results/{mode}/Plots/Markers/`）

#### 多模态标志基因

| 文件名 | 降维空间 | 描述 |
|--------|---------|------|
| `multimodal_featureplot_stage_markers.png` | RNA UMAP | 阶段标志基因 |
| `multimodal_featureplot_stage_markers_wnn.png` ⭐ | WNN UMAP | 阶段标志基因（WNN） |

#### 扩展基因组（8组功能基因）

**RNA UMAP版本：**
- `extended_genes_dotplot_8groups.png`（DotPlot）
- `extended_genes_featureplot_{组名}.png`（7个FeaturePlot）

**WNN UMAP版本（新增）：**
- `extended_genes_featureplot_调控基因_wnn.png` ⭐
- `extended_genes_featureplot_迁移基因_wnn.png` ⭐
- `extended_genes_featureplot_标记基因_wnn.png` ⭐
- `extended_genes_featureplot_成熟基因_wnn.png` ⭐
- `extended_genes_featureplot_TH2反应基因_wnn.png` ⭐
- `extended_genes_featureplot_转录因子_wnn.png` ⭐
- `extended_genes_featureplot_发育门控_wnn.png` ⭐

#### ADT蛋白可视化

**RNA UMAP版本：**
- `adt_proteins_dotplot.png`（DotPlot）
- `adt_proteins_featureplot_batch{N}.png`（FeaturePlot）

**WNN UMAP版本（新增）：**
- `adt_proteins_featureplot_batch1_wnn.png` ⭐
- `adt_proteins_featureplot_batch2_wnn.png` ⭐（如果有）

### 3. 下游分析数据（`results/{mode}/`）

| 文件路径 | 包含WNN UMAP | 用途 |
|---------|-------------|------|
| `Trajectory/monocle3_data.rds` | ✅ | Monocle3轨迹分析 |
| `CellOracle/celloracle_data.rds` | ✅ | CellOracle基因调控网络 |
| `Exported_Data/wnn_umap_coordinates.csv` | ✅ | WNN UMAP坐标导出 |

---

## 🎨 可视化参数总结

| 可视化类型 | RNA UMAP | WNN UMAP | 颜色方案 | 布局 |
|-----------|---------|---------|---------|------|
| 作者注释 DimPlot | ✅ | ✅ | 分类颜色 | 12x10 |
| WNN三联对比 | - | ✅ | 分类颜色 | 24x8 |
| 基因 FeaturePlot | ✅ | ✅ | lightgrey→red | 3x3 |
| ADT FeaturePlot | ✅ | ✅ | lightgrey→blue | 3x3 |
| 基因 DotPlot | ✅ | - | RdYlBu | 自适应 |
| ADT DotPlot | ✅ | - | RdYlBu | 自适应 |

---

## 🔍 关键技术特性

### 1. 条件检测 WNN UMAP 可用性

所有 WNN UMAP 可视化都包含条件检测：

```r
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  # 生成 WNN UMAP 可视化
} else {
  cat("⚠️ 当前对象无 WNN UMAP 降维结果\n")
}
```

**优势：**
- 兼容有/无 WNN UMAP 的数据集
- 避免运行时错误
- 提供清晰的用户反馈

### 2. 文件命名规范

| 原文件名 | WNN 版本文件名 | 规则 |
|---------|---------------|------|
| `xxx.png` | `xxx_wnn.png` | 添加 `_wnn` 后缀 |
| `xxx_batch1.png` | `xxx_batch1_wnn.png` | 在批次号后添加 `_wnn` |

**优势：**
- 易于识别和对比
- 保持文件组织一致性

### 3. 标题和副标题标注

**RNA UMAP：**
```r
title = "基因表达 - RNA UMAP"
subtitle = "RNA降维空间 | ..."
```

**WNN UMAP：**
```r
title = "基因表达 - WNN UMAP"
subtitle = "WNN整合空间 | ..."
```

**优势：**
- 用户可立即区分降维空间
- 避免混淆

---

## 📋 使用指南

### 运行脚本

```r
# 1. 设置分析模式
analysis_mode <- "cDC1"  # 或 "cDC2" 或 "integrated"

# 2. 运行脚本
rmarkdown::render("4.downstream_analysis.Rmd")
```

### 验证 WNN UMAP 支持

```r
# 读取输出的 Seurat 对象
seurat_obj <- readRDS("results/cDC1/Robjects/seurat_obj_annotated.rds")

# 检查是否包含 WNN UMAP
"wnn_umap" %in% names(seurat_obj@reductions)  # 应该返回 TRUE

# 查看降维结果
names(seurat_obj@reductions)
# [1] "pca"      "umap"     "adt_pca"  "adt_umap" "wnn_umap"

# 可视化对比
library(Seurat)
library(ggplot2)

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "author_annotation")
p2 <- DimPlot(seurat_obj, reduction = "wnn_umap", group.by = "author_annotation")
p1 | p2
```

### 使用下游分析数据

#### Monocle3 轨迹分析

```r
# 读取 Monocle3 数据
monocle3_data <- readRDS("results/cDC1/Trajectory/monocle3_data.rds")

# 访问 WNN UMAP 坐标
wnn_coords <- monocle3_data$wnn_umap
head(wnn_coords)

# 可视化
library(ggplot2)
ggplot(wnn_coords, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = monocle3_data$metadata$author_annotation)) +
  labs(title = "cDC1 Developmental Trajectory on WNN UMAP",
       color = "Cell Type")
```

#### CellOracle 基因调控网络

```python
# Python 中使用
import pandas as pd
import pickle

# 读取 CellOracle 数据
with open("results/cDC1/CellOracle/celloracle_data.rds", "rb") as f:
    celloracle_data = pickle.load(f)

# 访问 WNN UMAP 坐标
wnn_umap = celloracle_data['wnn_umap_coords']

# 可视化
import matplotlib.pyplot as plt
plt.scatter(wnn_umap.iloc[:, 0], wnn_umap.iloc[:, 1],
            c=celloracle_data['metadata']['author_annotation'])
plt.title("cDC1 on WNN UMAP")
plt.xlabel("wnnUMAP_1")
plt.ylabel("wnnUMAP_2")
plt.show()
```

---

## ⚠️ 注意事项

### 1. 输入文件要求

**必须使用 `*_final.rds` 才有 WNN UMAP：**

| 文件名 | WNN UMAP | 来源 |
|--------|---------|------|
| `cDC1_final.rds` ✅ | 有 | CITE-seq_Complete_Analysis_v3.Rmd Section 8 |
| `cDC1_pre_wnn.rds` ❌ | 无 | CITE-seq_Complete_Analysis_v3.Rmd Section 5 |

脚本已自动优先加载 `*_final.rds`，但如果该文件不存在，会回退到 `*_pre_wnn.rds` 并跳过 WNN UMAP 可视化。

### 2. 作者注释来源

作者注释从缓存文件加载：`results/Cache/{mode}_01_subset_extracted.rds`

确保以下缓存文件存在：
- cDC1: `results/Cache/cDC1_01_subset_extracted.rds`
- cDC2: `results/Cache/cDC2_01_subset_extracted.rds`
- integrated: `results/Cache/integrated_01_subset_extracted.rds`

### 3. 三种 UMAP 的区别

| UMAP类型 | 降维空间 | 反映的信息 | 适用场景 |
|---------|---------|-----------|---------|
| **RNA UMAP** (`umap`) | 基于基因表达 | 转录组状态 | 基因表达模式分析 |
| **ADT UMAP** (`adt_umap`) | 基于蛋白表达 | 蛋白表型 | 表面标志蛋白分析 |
| **WNN UMAP** (`wnn_umap`) | 整合 RNA+ADT | 多模态整合 | 综合评估细胞状态 ⭐ |

**推荐使用场景：**
- **RNA UMAP**: 评估基因表达驱动的细胞分群
- **WNN UMAP**: 评估综合多模态信息的细胞分群（更准确）⭐
- **对比两者**: 发现RNA和蛋白不一致的细胞亚群

---

## 🎯 核心改进总结

### 修改位置（按章节）

| 章节 | 行号 | 修改内容 | WNN UMAP 支持 |
|------|------|---------|--------------|
| 输入文件加载 | 221-226 | 优先加载 final.rds | ✅ 基础 |
| 可视化作者注释 | 1633-1697 | WNN UMAP 作者注释可视化 | ✅ 新增 |
| 多模态 Marker 可视化 | 2028-2071 | 阶段标志基因 FeaturePlot | ✅ 新增 |
| 扩展基因组可视化（基因） | 2191-2251 | 8组基因 FeaturePlot | ✅ 新增 |
| 扩展基因组可视化（ADT） | 2276-2330 | ADT蛋白 FeaturePlot | ✅ 新增 |
| Monocle3 数据准备 | 2344-2356 | 添加 wnn_umap 字段 | ✅ 新增 |
| CellOracle 数据准备 | 2381-2393 | 添加 wnn_umap_coords 字段 | ✅ 新增 |
| CSV 导出 | 2442-2458 | 导出 WNN UMAP 坐标 | ✅ 新增 |

### 关键功能

1. ✅ **输入文件优先级修复** - 确保加载包含 WNN UMAP 的对象
2. ✅ **WNN UMAP 作者注释可视化** - 在 WNN-combined Clustering 上显示真实细胞类型
3. ✅ **所有 FeaturePlot 支持 WNN UMAP** - 基因、标志基因、ADT蛋白全覆盖
4. ✅ **下游分析工具 WNN UMAP 支持** - Monocle3、CellOracle 可使用 WNN UMAP
5. ✅ **WNN UMAP 坐标导出** - CSV 格式方便外部工具使用

---

## 🚀 预期效果

运行修改后的脚本，您将获得：

### 1. 在 WNN-combined Clustering 图上：
- ✅ 可以看到 WNN_clusters（0,1,2,3...）
- ✅ 可以看到作者注释（Pre-cDC1s, Mature 等）
- ✅ 可以对比两者的一致性

### 2. 所有基因/蛋白可视化都有 WNN UMAP 版本：
- ✅ 阶段标志基因
- ✅ 8组功能基因（61个基因）
- ✅ ADT蛋白（12个蛋白）

### 3. 下游分析工具可使用 WNN UMAP：
- ✅ Monocle3 轨迹分析
- ✅ CellOracle 基因调控网络
- ✅ 其他自定义分析（通过 CSV 导出）

---

## 📚 相关文档

- `WNN_作者注释分析报告.md` - WNN 注释功能的可行性分析
- `WNN_UMAP_作者注释可视化说明.md` - WNN UMAP 作者注释功能详细说明
- `4.downstream_analysis_modifications_summary.md` - 下游分析脚本修改总结（包括8组基因和ADT蛋白）

---

## ✨ 总结

**所有可视化模块现已完整支持 WNN UMAP！**

- ✅ 作者注释可视化
- ✅ 多模态标志基因可视化
- ✅ 扩展基因组可视化（基因 + ADT蛋白）
- ✅ 下游分析数据准备

**核心优势：**
1. **完整性** - 涵盖所有主要可视化模块
2. **一致性** - 统一的条件检测和文件命名规范
3. **可用性** - 自动检测 WNN UMAP 可用性，兼容不同数据集
4. **可扩展性** - 下游分析工具可直接使用 WNN UMAP 数据

**用户可以自由选择在 RNA UMAP 或 WNN UMAP 上查看任何分析结果！** 🎉
