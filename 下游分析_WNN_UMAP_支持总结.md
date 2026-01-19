# 下游分析 WNN UMAP 支持 - 完整总结

## 📅 更新日期
2026-01-19

## 🎯 更新目标

在所有涉及UMAP可视化和下游分析的地方，全面添加**WNN UMAP支持**，确保：
1. 可以在WNN UMAP上显示作者注释
2. 所有下游分析工具都能使用WNN UMAP坐标
3. 为轨迹分析等准备的RDS文件包含WNN UMAP

---

## ✅ 完成的修改

### 1. 修复输入文件优先级（第221-226行）

**问题**：
- 原先优先加载 `*_pre_wnn.rds`（WNN分析前，无wnn_umap）
- 导致无法在WNN UMAP上显示作者注释

**修复**：
```r
# 修改前 ❌
if (file.exists(input_file_prewnn)) {
  input_file <- input_file_prewnn  # 优先pre_wnn
}

# 修改后 ✅
if (file.exists(input_file_final)) {
  input_file <- input_file_final  # 优先final (包含WNN UMAP)
}
```

**效果**：
- ✅ 现在优先加载包含WNN UMAP的 `*_final.rds`
- ✅ 确保WNN UMAP相关功能可用

---

### 2. 添加WNN UMAP上的作者注释可视化（第1633-1697行）

**新增内容**：

#### 2.1 单图：WNN UMAP + 作者注释
```r
DimPlot(seurat_obj, reduction = "wnn_umap",
        group.by = "author_annotation",
        label = TRUE)
```
- 文件：`author_annotation_wnn_umap.png/.pdf`
- 功能：在WNN UMAP上直接显示作者的细胞类型注释

#### 2.2 三联对比图
```r
p_wnn_clusters | p_wnn_annot | p_wnn_detailed
```
- 文件：`wnn_umap_triple_comparison.png/.pdf`
- 功能：对比三种视图
  - WNN Clusters（数字0,1,2...）
  - **Author Annotation**（细胞类型名称）⭐
  - Detailed Annotation（详细分类）

**效果**：
- ✅ 可以在WNN-combined Clustering图上看到作者注释
- ✅ 对比WNN聚类与作者注释的一致性

---

### 3. Monocle3数据准备 - 添加WNN UMAP（第2344-2356行）

**新增内容**：
```r
# 如果有 WNN UMAP 结果
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  monocle3_data$wnn_umap <- Embeddings(seurat_obj, "wnn_umap")
  cat("  ✓ 已添加 WNN UMAP 坐标\n")
}
```

**文件**：`data_for_monocle3.rds`

**内容变化**：
```
修改前：
- counts, normalized, metadata, umap, cell_type, pca

修改后：
- counts, normalized, metadata, umap, wnn_umap ✅, cell_type, pca
```

**用途**：
```r
# Monocle3拟时序分析可以使用WNN UMAP
monocle3_data <- readRDS('data_for_monocle3.rds')

# 方法1：使用RNA UMAP
reducedDims(cds)$UMAP <- monocle3_data$umap

# 方法2：使用WNN UMAP（推荐）✨
if (!is.null(monocle3_data$wnn_umap)) {
  reducedDims(cds)$UMAP <- monocle3_data$wnn_umap
}
```

---

### 4. CellOracle数据准备 - 添加WNN UMAP（第2381-2393行）

**新增内容**：
```r
# 如果有 WNN UMAP 结果
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  celloracle_data$wnn_umap_coords <- as.data.frame(Embeddings(seurat_obj, "wnn_umap"))
  cat("  ✓ 已添加 WNN UMAP 坐标\n")
}
```

**文件**：`data_for_celloracle.rds`

**内容变化**：
```
修改前：
- expression_matrix, metadata, umap_coords, cell_type, variable_genes

修改后：
- expression_matrix, metadata, umap_coords, wnn_umap_coords ✅, cell_type, variable_genes
```

**用途**：
```python
# CellOracle基因调控网络可以使用WNN UMAP
import anndata

# 添加RNA UMAP
adata.obsm['X_umap'] = celloracle_data['umap_coords'].values

# 添加WNN UMAP（如果有）✨
if 'wnn_umap_coords' in celloracle_data:
    adata.obsm['X_wnn_umap'] = celloracle_data['wnn_umap_coords'].values
    # 使用WNN UMAP进行可视化
    adata.obsm['X_umap'] = adata.obsm['X_wnn_umap']
```

---

### 5. 通用CSV导出 - 添加WNN UMAP坐标（第2442-2458行）

**新增内容**：
```r
# 🆕 WNN UMAP 坐标（如果有）
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  wnn_umap_df <- as.data.frame(Embeddings(seurat_obj, "wnn_umap"))
  wnn_umap_df$cell_id <- rownames(wnn_umap_df)
  if (has_author_annotation) {
    wnn_umap_df$cell_type <- seurat_obj$author_annotation
  }
  wnn_umap_df$cluster <- seurat_obj$SCT_clusters
  if ("WNN_clusters" %in% colnames(seurat_obj@meta.data)) {
    wnn_umap_df$WNN_cluster <- seurat_obj$WNN_clusters
  }

  write.csv(wnn_umap_df,
            file.path(output_dir, "Tables", "wnn_umap_coordinates.csv"),
            row.names = FALSE)
}
```

**新增文件**：`wnn_umap_coordinates.csv`

**文件结构**：
```csv
cell_id,wnnUMAP_1,wnnUMAP_2,cell_type,cluster,WNN_cluster
CELL_1,-2.5,1.3,Pre-cDC1s,0,0
CELL_2,-2.4,1.5,Pre-cDC1s,0,0
CELL_3,3.2,-1.8,Mature cDC1s,5,7
...
```

**用途**：
- 通用的WNN UMAP坐标文件
- 可用于Python、R等任意工具
- 包含细胞类型注释和聚类信息

---

### 6. 下游分析使用指南更新（第2506-2596行）

**新增内容**：

#### 6.1 Monocle3示例更新
```r
# 方法2：使用 WNN UMAP（推荐，如果有）
if (!is.null(monocle3_data$wnn_umap)) {
  reducedDims(cds)$UMAP <- monocle3_data$wnn_umap
  cat('使用 WNN UMAP（多模态整合）\n')
}
```

#### 6.2 CellOracle示例更新
```python
# 添加 WNN UMAP（如果有）
if 'wnn_umap_coords' in celloracle_data:
    adata.obsm['X_wnn_umap'] = celloracle_data['wnn_umap_coords'].values
    # 使用WNN UMAP进行可视化
    adata.obsm['X_umap'] = adata.obsm['X_wnn_umap']
```

#### 6.3 新增第4节：通用WNN UMAP坐标使用
```r
# 从CSV读取WNN UMAP坐标
wnn_umap <- read.csv('wnn_umap_coordinates.csv')

# 可视化
library(ggplot2)
ggplot(wnn_umap, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = cell_type)) +
  geom_point(size = 0.5) +
  theme_classic() +
  labs(title = 'WNN UMAP - Author Annotation')
```

---

## 📊 输出文件清单

运行修改后的脚本，将在以下位置生成新文件：

### 可视化图片（`results/{mode}/Plots/Annotation/`）

**RNA UMAP（原有）**：
1. `author_annotation_umap.png/.pdf`
2. `author_annotation_detailed_umap.png/.pdf`
3. `clustering_vs_annotation_comparison.png/.pdf`

**WNN UMAP（新增）**：
4. ✅ **`author_annotation_wnn_umap.png/.pdf`** - WNN UMAP + 作者注释
5. ✅ **`wnn_umap_triple_comparison.png/.pdf`** - 三联对比图

### 下游分析数据（`results/{mode}/Robjects/`）

**更新的文件**：
1. ✅ **`data_for_monocle3.rds`** - 新增wnn_umap字段
2. ✅ **`data_for_celloracle.rds`** - 新增wnn_umap_coords字段

### 通用CSV（`results/{mode}/Tables/`）

**新增文件**：
3. ✅ **`wnn_umap_coordinates.csv`** - WNN UMAP坐标+注释

---

## 🎯 三种UMAP的对比

| UMAP类型 | 降维空间 | 反映的信息 | 适用场景 |
|---------|---------|-----------|----------|
| **RNA UMAP** | 基于基因表达（SCT） | 转录组状态 | 基因表达驱动的分析 |
| **ADT UMAP** | 基于蛋白表达 | 蛋白表型 | 表面标记驱动的分析 |
| **WNN UMAP** ⭐ | 整合RNA+ADT | 多模态整合 | **推荐用于综合分析** |

### 为什么优先使用WNN UMAP？

1. **多模态整合**：同时考虑基因和蛋白信息
2. **更稳健**：减少单一模态的技术噪音
3. **生物学意义更强**：整合转录和蛋白表型
4. **适合CITE-seq**：充分利用多模态数据优势

---

## 🔧 使用方法

### 在RStudio中运行

```r
# 设置分析模式
analysis_mode <- "cDC1"  # 或 "cDC2" 或 "integrated"

# 运行脚本
rmarkdown::render("4.downstream_analysis.Rmd")
```

### 预期输出（控制台）

```
🔍 检查输入文件...
  ✅ 找到 final.rds（包含WNN UMAP）: results/Robjects/GSE228544_cDC1_final.rds

...

📊 可视化作者注释...

【WNN UMAP可视化】  ✅
✓ WNN UMAP作者注释可视化完成

...

【2/5】为 Monocle3 拟时序分析准备数据...
  ✓ 已添加 WNN UMAP 坐标  ✅
  ✓ data_for_monocle3.rds 已保存
  包含: counts, normalized, metadata, umap, wnn_umap, cell_type, pca

【3/5】为 CellOracle 基因调控网络分析准备数据...
  ✓ 已添加 WNN UMAP 坐标  ✅
  ✓ data_for_celloracle.rds 已保存
  包含: expression_matrix, metadata, umap_coords, wnn_umap_coords, cell_type, variable_genes

【5/5】导出通用格式数据...
  ✓ cell_metadata.csv 已保存
  ✓ umap_coordinates.csv 已保存 (RNA UMAP)
  ✓ wnn_umap_coordinates.csv 已保存 (WNN UMAP)  ✅
```

---

## ⚠️ 注意事项

### 1. WNN UMAP的前提条件

WNN UMAP只在完整运行了WNN分析的数据中存在：

**有WNN UMAP** ✅：
- `*_final.rds`（来自 `CITE-seq_Complete_Analysis_v3.Rmd` Section 8）

**无WNN UMAP** ❌：
- `*_pre_wnn.rds`（WNN分析前的对象）

### 2. 自动降级处理

脚本会自动检测WNN UMAP是否存在：

```r
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  # 使用WNN UMAP
} else {
  cat("⚠️ 当前对象无WNN UMAP降维结果\n")
  # 降级使用RNA UMAP
}
```

### 3. 文件优先级

修改后的加载顺序：
1. ✅ 优先：`*_final.rds`（包含WNN UMAP）
2. 备选：`*_pre_wnn.rds`（无WNN UMAP，会提示警告）

### 4. 下游分析工具支持

| 工具 | RNA UMAP | WNN UMAP | 推荐 |
|------|---------|----------|------|
| Monocle3 | ✅ | ✅ | WNN UMAP |
| CellOracle | ✅ | ✅ | WNN UMAP |
| scTenifoldKnk | ✅ | N/A | RNA |
| 通用可视化 | ✅ | ✅ | WNN UMAP |

---

## 📝 代码示例

### 示例1：Monocle3拟时序分析（使用WNN UMAP）

```r
library(monocle3)

# 加载数据
monocle3_data <- readRDS("results/cDC1/Robjects/data_for_monocle3.rds")

# 创建CellDataSet
cds <- new_cell_data_set(
  expression_data = monocle3_data$counts,
  cell_metadata = monocle3_data$metadata
)

# 使用WNN UMAP（推荐）
if (!is.null(monocle3_data$wnn_umap)) {
  reducedDims(cds)$UMAP <- monocle3_data$wnn_umap
  cat("✓ 使用WNN UMAP进行轨迹分析\n")
} else {
  reducedDims(cds)$UMAP <- monocle3_data$umap
  cat("⚠️ 使用RNA UMAP（无WNN UMAP）\n")
}

# 学习轨迹
cds <- learn_graph(cds)
cds <- order_cells(cds)

# 可视化（在WNN UMAP空间）
plot_cells(cds, color_cells_by = "pseudotime",
           label_cell_groups = FALSE)

plot_cells(cds, color_cells_by = "cell_type",
           label_cell_groups = FALSE)
```

### 示例2：直接使用WNN UMAP坐标可视化

```r
library(ggplot2)

# 读取WNN UMAP坐标
wnn_umap <- read.csv("results/cDC1/Tables/wnn_umap_coordinates.csv")

# 按作者注释染色
ggplot(wnn_umap, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = cell_type)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_classic() +
  labs(title = "WNN UMAP - Author Annotation",
       x = "WNN UMAP 1", y = "WNN UMAP 2") +
  theme(legend.position = "right")

# 按WNN聚类染色
ggplot(wnn_umap, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = as.factor(WNN_cluster))) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_classic() +
  labs(title = "WNN UMAP - WNN Clusters",
       x = "WNN UMAP 1", y = "WNN UMAP 2",
       color = "WNN Cluster")
```

### 示例3：对比RNA UMAP vs WNN UMAP

```r
library(patchwork)

# RNA UMAP
rna_umap <- read.csv("results/cDC1/Tables/umap_coordinates.csv")
p1 <- ggplot(rna_umap, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.5) +
  theme_classic() +
  labs(title = "RNA UMAP")

# WNN UMAP
wnn_umap <- read.csv("results/cDC1/Tables/wnn_umap_coordinates.csv")
p2 <- ggplot(wnn_umap, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = cell_type)) +
  geom_point(size = 0.5) +
  theme_classic() +
  labs(title = "WNN UMAP")

# 并排对比
p1 | p2
```

---

## 🎉 总结

### 修改的文件
- `4.downstream_analysis.Rmd`

### 修改的代码块
1. ✅ 第221-226行：修复输入文件优先级
2. ✅ 第1633-1697行：添加WNN UMAP作者注释可视化
3. ✅ 第2344-2356行：Monocle3添加wnn_umap
4. ✅ 第2381-2393行：CellOracle添加wnn_umap_coords
5. ✅ 第2442-2458行：CSV导出wnn_umap坐标
6. ✅ 第2506-2596行：更新下游分析使用指南

### 新增的功能
- ✅ WNN UMAP上的作者注释可视化
- ✅ WNN UMAP三联对比图
- ✅ Monocle3支持WNN UMAP
- ✅ CellOracle支持WNN UMAP
- ✅ WNN UMAP坐标CSV导出
- ✅ 完整的WNN UMAP使用指南

### 影响的下游分析
- ✅ 拟时序分析（Monocle3/Slingshot）
- ✅ 基因调控网络（CellOracle）
- ✅ 轨迹推断（任意工具）
- ✅ 通用可视化

---

**现在所有涉及UMAP的分析都支持WNN UMAP了！** 🎊
