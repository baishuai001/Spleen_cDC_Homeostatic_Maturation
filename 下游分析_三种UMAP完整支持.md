# 下游分析准备 - 三种UMAP空间完整支持

## 📅 修改日期
2026-01-19

## 🎯 修改目标

为"# 🔥 下游分析准备"模块添加**ADT UMAP**支持，确保所有下游分析工具都能使用**三种UMAP空间**：

1. **RNA UMAP** - 基于基因表达的降维空间
2. **ADT UMAP** - 基于蛋白表达的降维空间 ⭐
3. **WNN UMAP** - 整合RNA+ADT的多模态空间

---

## 📋 修改内容

### 1️⃣ Monocle3 拟时序分析数据

**位置：** 第2452-2498行

**功能：** 为Monocle3添加三种UMAP坐标，支持在不同降维空间进行轨迹分析。

#### 修改前

```r
monocle3_data <- list(
  counts = ...,
  normalized = ...,
  metadata = ...,
  umap = Embeddings(seurat_obj, "umap"),  # 仅RNA UMAP
  cell_type = ...,
  pca = ...
)

# 只有WNN UMAP
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  monocle3_data$wnn_umap <- Embeddings(seurat_obj, "wnn_umap")
}
```

#### 修改后

```r
monocle3_data <- list(
  counts = ...,
  normalized = ...,
  metadata = ...,
  umap = Embeddings(seurat_obj, "umap"),  # RNA UMAP
  cell_type = ...,
  pca = ...
)

# 🆕 添加 ADT UMAP
if ("adt_umap" %in% names(seurat_obj@reductions)) {
  monocle3_data$adt_umap <- Embeddings(seurat_obj, "adt_umap")
  cat("  ✓ 已添加 ADT UMAP 坐标\n")
}

# 添加 WNN UMAP
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  monocle3_data$wnn_umap <- Embeddings(seurat_obj, "wnn_umap")
  cat("  ✓ 已添加 WNN UMAP 坐标\n")
}

# 动态显示包含的内容
content_list <- c("counts, normalized, metadata, umap, cell_type, pca")
if ("adt_umap" %in% names(monocle3_data)) {
  content_list <- paste0(content_list, ", adt_umap")
}
if ("wnn_umap" %in% names(monocle3_data)) {
  content_list <- paste0(content_list, ", wnn_umap")
}
cat("  包含:", content_list, "\n\n")
```

#### 输出文件

**文件名：** `results/{mode}/Robjects/data_for_monocle3.rds`

**包含内容：**
- `counts` - 原始计数矩阵
- `normalized` - 归一化表达矩阵
- `metadata` - 细胞元数据
- `umap` - RNA UMAP坐标
- `adt_umap` - **ADT UMAP坐标** ⭐（新增）
- `wnn_umap` - WNN UMAP坐标
- `cell_type` - 细胞类型注释
- `pca` - PCA结果

---

### 2️⃣ CellOracle 基因调控网络分析数据

**位置：** 第2492-2537行

**功能：** 为CellOracle添加三种UMAP坐标，支持在不同降维空间进行基因调控网络分析和可视化。

#### 修改前

```r
celloracle_data <- list(
  expression_matrix = ...,
  metadata = ...,
  umap_coords = as.data.frame(Embeddings(seurat_obj, "umap")),  # 仅RNA UMAP
  cell_type = ...,
  variable_genes = ...
)

# 只有WNN UMAP
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  celloracle_data$wnn_umap_coords <- as.data.frame(Embeddings(seurat_obj, "wnn_umap"))
}
```

#### 修改后

```r
celloracle_data <- list(
  expression_matrix = ...,
  metadata = ...,
  umap_coords = as.data.frame(Embeddings(seurat_obj, "umap")),  # RNA UMAP
  cell_type = ...,
  variable_genes = ...
)

# 🆕 添加 ADT UMAP
if ("adt_umap" %in% names(seurat_obj@reductions)) {
  celloracle_data$adt_umap_coords <- as.data.frame(Embeddings(seurat_obj, "adt_umap"))
  cat("  ✓ 已添加 ADT UMAP 坐标\n")
}

# 添加 WNN UMAP
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  celloracle_data$wnn_umap_coords <- as.data.frame(Embeddings(seurat_obj, "wnn_umap"))
  cat("  ✓ 已添加 WNN UMAP 坐标\n")
}

# 动态显示包含的内容
content_list <- "expression_matrix, metadata, umap_coords, cell_type, variable_genes"
if ("adt_umap_coords" %in% names(celloracle_data)) {
  content_list <- paste0(content_list, ", adt_umap_coords")
}
if ("wnn_umap_coords" %in% names(celloracle_data)) {
  content_list <- paste0(content_list, ", wnn_umap_coords")
}
cat("  包含:", content_list, "\n\n")
```

#### 输出文件

**文件名：** `results/{mode}/Robjects/data_for_celloracle.rds`

**包含内容：**
- `expression_matrix` - 归一化表达矩阵
- `metadata` - 细胞元数据
- `umap_coords` - RNA UMAP坐标（DataFrame格式）
- `adt_umap_coords` - **ADT UMAP坐标** ⭐（新增）
- `wnn_umap_coords` - WNN UMAP坐标
- `cell_type` - 细胞类型注释
- `variable_genes` - 高变基因列表

---

### 3️⃣ CSV 通用格式导出

**位置：** 第2554-2607行

**功能：** 导出三种UMAP坐标的CSV文件，方便在Python、Excel等工具中使用。

#### 修改前

```r
# RNA UMAP 坐标
write.csv(umap_df, "umap_coordinates.csv")

# 只有WNN UMAP
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  write.csv(wnn_umap_df, "wnn_umap_coordinates.csv")
}
```

#### 修改后

```r
# RNA UMAP 坐标
umap_df <- as.data.frame(Embeddings(seurat_obj, "umap"))
umap_df$cell_id <- rownames(umap_df)
umap_df$cell_type <- seurat_obj$author_annotation
umap_df$cluster <- seurat_obj$SCT_clusters
write.csv(umap_df, "umap_coordinates.csv")
cat("  ✓ umap_coordinates.csv 已保存 (RNA UMAP)\n")

# 🆕 ADT UMAP 坐标（新增）
if ("adt_umap" %in% names(seurat_obj@reductions)) {
  adt_umap_df <- as.data.frame(Embeddings(seurat_obj, "adt_umap"))
  adt_umap_df$cell_id <- rownames(adt_umap_df)
  adt_umap_df$cell_type <- seurat_obj$author_annotation
  adt_umap_df$cluster <- seurat_obj$SCT_clusters
  write.csv(adt_umap_df, "adt_umap_coordinates.csv")
  cat("  ✓ adt_umap_coordinates.csv 已保存 (ADT UMAP)\n")
}

# WNN UMAP 坐标
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  wnn_umap_df <- as.data.frame(Embeddings(seurat_obj, "wnn_umap"))
  wnn_umap_df$cell_id <- rownames(wnn_umap_df)
  wnn_umap_df$cell_type <- seurat_obj$author_annotation
  wnn_umap_df$cluster <- seurat_obj$SCT_clusters
  wnn_umap_df$WNN_cluster <- seurat_obj$WNN_clusters  # WNN聚类
  write.csv(wnn_umap_df, "wnn_umap_coordinates.csv")
  cat("  ✓ wnn_umap_coordinates.csv 已保存 (WNN UMAP)\n")
}
```

#### 输出文件

**1. RNA UMAP坐标：** `results/{mode}/Tables/umap_coordinates.csv`

| UMAP_1 | UMAP_2 | cell_id | cell_type | cluster |
|--------|--------|---------|-----------|---------|
| -2.34 | 5.67 | AAACCCAA... | Pre-cDC1s | 0 |
| 3.21 | -1.45 | AAACCCAB... | Mature cDC1 | 2 |

**2. ADT UMAP坐标：** `results/{mode}/Tables/adt_umap_coordinates.csv` ⭐

| adtUMAP_1 | adtUMAP_2 | cell_id | cell_type | cluster |
|-----------|-----------|---------|-----------|---------|
| -3.12 | 4.56 | AAACCCAA... | Pre-cDC1s | 0 |
| 2.78 | -2.34 | AAACCCAB... | Mature cDC1 | 2 |

**3. WNN UMAP坐标：** `results/{mode}/Tables/wnn_umap_coordinates.csv`

| wnnUMAP_1 | wnnUMAP_2 | cell_id | cell_type | cluster | WNN_cluster |
|-----------|-----------|---------|-----------|---------|-------------|
| -2.89 | 5.12 | AAACCCAA... | Pre-cDC1s | 0 | 0 |
| 3.01 | -1.89 | AAACCCAB... | Mature cDC1 | 2 | 3 |

---

## 📂 完整输出文件结构

运行修改后的脚本，会在 `results/{mode}/` 目录下生成：

```
results/cDC1/
├── Robjects/
│   ├── seurat_obj_annotated.rds          # 带注释的完整Seurat对象
│   ├── data_for_monocle3.rds             # Monocle3数据（含3种UMAP）✅
│   ├── data_for_celloracle.rds           # CellOracle数据（含3种UMAP）✅
│   ├── data_for_sctenifoldknk.rds        # scTenifoldKnk数据
│   └── gene_lists_reference.rds          # 基因列表参考
│
└── Tables/
    ├── cell_metadata.csv                  # 细胞元数据
    ├── umap_coordinates.csv               # RNA UMAP坐标
    ├── adt_umap_coordinates.csv           # 🆕 ADT UMAP坐标 ⭐
    └── wnn_umap_coordinates.csv           # WNN UMAP坐标
```

---

## 🎯 使用场景示例

### 场景 1：Monocle3 拟时序分析 - 选择不同UMAP空间

```r
library(monocle3)

# 加载数据
data <- readRDS("results/cDC1/Robjects/data_for_monocle3.rds")

# 创建Monocle3对象
cds <- new_cell_data_set(
  expression_data = data$counts,
  cell_metadata = data$metadata
)

# 方案1：使用RNA UMAP（基因表达驱动的轨迹）
reducedDims(cds)$UMAP <- data$umap
plot_cells(cds, color_cells_by = "author_annotation",
           label_groups_by_cluster = FALSE) +
  labs(title = "cDC1 Trajectory on RNA UMAP")

# 方案2：使用ADT UMAP（蛋白表达驱动的轨迹）⭐
reducedDims(cds)$UMAP <- data$adt_umap
plot_cells(cds, color_cells_by = "author_annotation",
           label_groups_by_cluster = FALSE) +
  labs(title = "cDC1 Trajectory on ADT UMAP")

# 方案3：使用WNN UMAP（整合RNA+蛋白的轨迹，最准确）
reducedDims(cds)$UMAP <- data$wnn_umap
plot_cells(cds, color_cells_by = "author_annotation",
           label_groups_by_cluster = FALSE) +
  labs(title = "cDC1 Trajectory on WNN UMAP")

# 推断轨迹
cds <- learn_graph(cds)
cds <- order_cells(cds)

# 可视化拟时间
plot_cells(cds, color_cells_by = "pseudotime")
```

---

### 场景 2：CellOracle 基因调控网络 - 在不同UMAP空间可视化

```python
import celloracle as co
import pandas as pd
import pickle
import matplotlib.pyplot as plt

# 加载数据
with open("results/cDC1/Robjects/data_for_celloracle.rds", "rb") as f:
    data = pickle.load(f)

# 创建CellOracle对象
oracle = co.Oracle()
oracle.adata = sc.AnnData(
    X=data['expression_matrix'],
    obs=data['metadata']
)

# 方案1：在RNA UMAP上可视化GRN
oracle.adata.obsm['X_umap'] = data['umap_coords'].values
fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.umap(oracle.adata, color='author_annotation', ax=ax,
           title='GRN on RNA UMAP')

# 方案2：在ADT UMAP上可视化GRN ⭐
oracle.adata.obsm['X_umap'] = data['adt_umap_coords'].values
fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.umap(oracle.adata, color='author_annotation', ax=ax,
           title='GRN on ADT UMAP')

# 方案3：在WNN UMAP上可视化GRN（最准确）
oracle.adata.obsm['X_umap'] = data['wnn_umap_coords'].values
fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.umap(oracle.adata, color='author_annotation', ax=ax,
           title='GRN on WNN UMAP')

# 构建基因调控网络
oracle.get_cluster_specific_TFdict_from_Links()
oracle.fit_GRN_for_simulation()

# 模拟Irf8敲除
oracle.simulate_shift(perturb_condition={"Irf8": 0.0})
oracle.calculate_embedding_shift(sigma_corr=0.05)

# 在WNN UMAP上可视化扰动效果
oracle.plot_quiver(color="author_annotation", scale=5)
```

---

### 场景 3：Python自定义分析 - 对比三种UMAP

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 读取三种UMAP坐标
rna_umap = pd.read_csv("results/cDC1/Tables/umap_coordinates.csv")
adt_umap = pd.read_csv("results/cDC1/Tables/adt_umap_coordinates.csv")
wnn_umap = pd.read_csv("results/cDC1/Tables/wnn_umap_coordinates.csv")

# 三联对比图
fig, axes = plt.subplots(1, 3, figsize=(24, 8))

# RNA UMAP
sns.scatterplot(data=rna_umap, x='UMAP_1', y='UMAP_2',
                hue='cell_type', s=10, alpha=0.7, ax=axes[0])
axes[0].set_title("RNA UMAP (Gene Expression Space)", fontsize=16)
axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# ADT UMAP
sns.scatterplot(data=adt_umap, x='adtUMAP_1', y='adtUMAP_2',
                hue='cell_type', s=10, alpha=0.7, ax=axes[1])
axes[1].set_title("ADT UMAP (Protein Expression Space) ⭐", fontsize=16)
axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# WNN UMAP
sns.scatterplot(data=wnn_umap, x='wnnUMAP_1', y='wnnUMAP_2',
                hue='cell_type', s=10, alpha=0.7, ax=axes[2])
axes[2].set_title("WNN UMAP (Integrated Space)", fontsize=16)
axes[2].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("三种UMAP对比图.png", dpi=300)
plt.show()
```

**分析要点：**
- 比较同一细胞在三种UMAP上的位置
- 发现RNA和蛋白信息不一致的细胞亚群
- 评估WNN整合的效果

---

## 🔍 三种UMAP在下游分析中的选择建议

### 1. Monocle3 拟时序分析

| 研究目标 | 推荐UMAP | 原因 |
|---------|---------|------|
| 基因表达驱动的发育轨迹 | RNA UMAP | 最直接反映转录组变化 |
| 蛋白表型驱动的发育轨迹 | ADT UMAP ⭐ | 表面标志蛋白定义的轨迹 |
| 综合评估发育轨迹 | WNN UMAP | 整合RNA+蛋白，最准确 |
| 对比不同轨迹 | 三种都用 | 发现RNA-蛋白差异 |

**推荐：优先使用WNN UMAP，参考RNA和ADT UMAP**

---

### 2. CellOracle 基因调控网络

| 研究目标 | 推荐UMAP | 原因 |
|---------|---------|------|
| 基因调控网络可视化 | RNA UMAP | 基因→基因调控 |
| 蛋白-基因关联分析 | ADT UMAP ⭐ | 蛋白表型相关的基因调控 |
| 综合GRN分析 | WNN UMAP | 整合信息，细胞分群更准确 |
| 扰动分析可视化 | WNN UMAP | 综合评估扰动效果 |

**推荐：优先使用WNN UMAP进行GRN分析和扰动模拟**

---

### 3. 自定义Python/R分析

| 分析类型 | 推荐UMAP | 原因 |
|---------|---------|------|
| 基因表达可视化 | RNA UMAP | 最直接 |
| 蛋白表达可视化 | ADT UMAP ⭐⭐⭐ | **蛋白在蛋白空间最准确** |
| 细胞分群可视化 | WNN UMAP | 最可靠的分群 |
| 三种UMAP对比 | 三种都用 | 发现多模态不一致性 |

**推荐：根据具体分析目标选择，蛋白分析必用ADT UMAP**

---

## 💡 关键技术特性

### 1. 条件检测UMAP可用性

所有UMAP添加都包含条件检测：

```r
# ADT UMAP条件检测
if ("adt_umap" %in% names(seurat_obj@reductions)) {
  # 添加ADT UMAP坐标
  monocle3_data$adt_umap <- Embeddings(seurat_obj, "adt_umap")
  cat("  ✓ 已添加 ADT UMAP 坐标\n")
}

# WNN UMAP条件检测
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  # 添加WNN UMAP坐标
  monocle3_data$wnn_umap <- Embeddings(seurat_obj, "wnn_umap")
  cat("  ✓ 已添加 WNN UMAP 坐标\n")
}
```

**优势：**
- 兼容有/无ADT数据的Seurat对象
- 兼容有/无WNN分析的对象
- 不会因缺少降维结果而报错

---

### 2. 动态内容显示

```r
# Monocle3数据动态显示
content_list <- c("counts, normalized, metadata, umap, cell_type, pca")
if ("adt_umap" %in% names(monocle3_data)) {
  content_list <- paste0(content_list, ", adt_umap")
}
if ("wnn_umap" %in% names(monocle3_data)) {
  content_list <- paste0(content_list, ", wnn_umap")
}
cat("  包含:", content_list, "\n\n")
```

**输出示例：**

```
【2/5】为 Monocle3 拟时序分析准备数据...
  ✓ 已添加 ADT UMAP 坐标
  ✓ 已添加 WNN UMAP 坐标
  ✓ data_for_monocle3.rds 已保存
  包含: counts, normalized, metadata, umap, cell_type, pca, adt_umap, wnn_umap
```

**优势：**
- 清楚显示数据包含哪些UMAP
- 便于用户确认数据完整性

---

### 3. 文件命名规范

| UMAP类型 | CSV文件名 | 列名前缀 |
|---------|----------|---------|
| RNA UMAP | `umap_coordinates.csv` | `UMAP_1`, `UMAP_2` |
| ADT UMAP ⭐ | `adt_umap_coordinates.csv` | `adtUMAP_1`, `adtUMAP_2` |
| WNN UMAP | `wnn_umap_coordinates.csv` | `wnnUMAP_1`, `wnnUMAP_2` |

**优势：**
- 清晰标识UMAP类型
- 便于批量处理和对比
- 避免混淆

---

## ⚠️ 注意事项

### 1. UMAP可用性检查

在使用下游分析数据前，先检查包含哪些UMAP：

```r
# 读取Monocle3数据
data <- readRDS("results/cDC1/Robjects/data_for_monocle3.rds")

# 检查可用的UMAP
names(data)
# [1] "counts"      "normalized"  "metadata"    "umap"
# [5] "cell_type"   "analysis_mode" "pca"       "adt_umap"
# [9] "wnn_umap"

# 检查ADT UMAP
"adt_umap" %in% names(data)  # TRUE - 可以使用

# 检查WNN UMAP
"wnn_umap" %in% names(data)  # TRUE - 可以使用
```

---

### 2. UMAP存在的前提条件

| UMAP类型 | 存在条件 |
|---------|---------|
| RNA UMAP | ✅ 总是存在（Seurat标准流程） |
| ADT UMAP | ⚠️ 需要CITE-seq数据 + 运行ADT降维分析 |
| WNN UMAP | ⚠️ 需要CITE-seq数据 + 运行WNN多模态整合 |

**如何确保有ADT和WNN UMAP：**

在 `CITE-seq_Complete_Analysis_v3.Rmd` 中：
- Section 4: ADT降维分析 → 生成 `adt_umap`
- Section 8: WNN多模态整合 → 生成 `wnn_umap`
- 使用 `*_final.rds` 输入文件（包含完整降维结果）

---

### 3. 输入文件要求

**必须使用 `*_final.rds` 才有完整的三种UMAP：**

| 文件名 | RNA UMAP | ADT UMAP | WNN UMAP |
|--------|---------|---------|---------|
| `cDC1_final.rds` ✅ | ✅ | ✅ | ✅ |
| `cDC1_pre_wnn.rds` ⚠️ | ✅ | ✅ | ❌ |
| `cDC1_03_clustered.rds` ⚠️ | ✅ | ❌ | ❌ |

脚本已自动优先加载 `*_final.rds`（第221-226行）。

---

## 📊 修改统计

| 模块 | 行号 | 修改内容 | 新增字段 |
|------|------|---------|---------|
| Monocle3数据 | 2473-2498 | 添加ADT UMAP + 动态显示 | `adt_umap` ⭐ |
| CellOracle数据 | 2515-2537 | 添加ADT UMAP + 动态显示 | `adt_umap_coords` ⭐ |
| CSV导出 | 2563-2607 | 添加ADT UMAP CSV文件 | `adt_umap_coordinates.csv` ⭐ |

---

## ✨ 总结

### 核心改进

✅ **Monocle3数据** - 支持三种UMAP轨迹分析
✅ **CellOracle数据** - 支持三种UMAP GRN可视化
✅ **CSV导出** - 导出三种UMAP坐标文件
✅ **动态显示** - 自动检测并显示包含的UMAP
✅ **兼容性** - 兼容有/无ADT、WNN数据的对象

### 关键优势

1. **完整性** - 所有下游分析工具都能使用三种UMAP
2. **灵活性** - 用户可根据分析目标选择最佳UMAP
3. **准确性** - ADT UMAP为蛋白相关分析提供最佳空间
4. **对比性** - 可对比三种UMAP发现多模态差异

### 使用价值

**对于研究者：**
- 在拟时序分析中选择最合适的降维空间
- 在基因调控网络分析中综合RNA和蛋白信息
- 在自定义分析中对比三种UMAP的差异

**对于CITE-seq数据：**
- 充分利用多模态数据优势
- ADT UMAP为蛋白驱动的分析提供最佳视角 ⭐
- WNN UMAP为综合分析提供最准确结果

---

## 📚 相关文档

- `三种UMAP空间完整支持总结.md` - FeaturePlot可视化的三种UMAP支持
- `下游分析准备模块_详细说明.md` - 下游分析模块详细说明
- `下游分析_WNN_UMAP_完整支持总结.md` - 之前的WNN UMAP支持总结

---

**现在，所有下游分析工具都能在RNA、ADT、WNN三种UMAP空间中进行分析！** 🎉
