# 三种UMAP空间完整支持总结

## 📅 修改日期
2026-01-19

## 🎯 修改目标

为 `4.downstream_analysis.Rmd` 的所有FeaturePlot可视化添加**三种UMAP空间**的完整支持：

1. **RNA UMAP** - 基于基因表达的降维空间
2. **ADT UMAP** - 基于蛋白表达的降维空间 ⭐
3. **WNN UMAP** - 整合RNA+ADT的多模态空间

---

## 🧬 三种UMAP空间的区别

### 1. RNA UMAP (`umap`)

**降维依据：** 基于高变基因的表达量

**反映信息：**
- 转录组状态
- 基因表达模式
- 细胞的转录异质性

**适用场景：**
- 分析基因表达驱动的细胞分群
- 研究转录调控网络
- 识别基因表达标志物

**示例：** 两群细胞在RNA水平高变基因表达不同，在RNA UMAP上会分开

---

### 2. ADT UMAP (`adt_umap`) ⭐

**降维依据：** 基于ADT（抗体标记蛋白）的表达量

**反映信息：**
- 细胞表面蛋白表型
- 蛋白标志物表达
- 细胞的蛋白表达异质性

**适用场景：**
- 分析表面标志蛋白驱动的细胞分群
- 研究蛋白表达模式
- 识别蛋白标志物

**特别重要：** 对于**ADT蛋白可视化**，ADT UMAP是最准确的降维空间！

**示例：** cDC1和cDC2在蛋白标志物（如CD172a、XCR1）表达不同，在ADT UMAP上会清晰分开

---

### 3. WNN UMAP (`wnn_umap`)

**降维依据：** 整合RNA高变基因 + ADT蛋白的加权最近邻（Weighted Nearest Neighbor）

**反映信息：**
- 综合转录组和蛋白组的细胞状态
- 多模态整合后的细胞异质性
- 更全面的细胞身份信息

**适用场景：**
- 多模态综合分析（最准确）
- 当RNA和蛋白信息不完全一致时
- 需要最可靠的细胞分群

**优势：** 整合了两种信息来源，比单独的RNA或ADT UMAP更准确！

**示例：** 某些细胞在RNA UMAP上混在一起，但ADT蛋白表达不同，WNN UMAP会综合两者信息将它们正确分开

---

## 🔄 三种UMAP的对比

| 特性 | RNA UMAP | ADT UMAP ⭐ | WNN UMAP |
|------|---------|-----------|----------|
| **数据来源** | 基因表达 | 蛋白表达 | RNA+ADT整合 |
| **信息维度** | 转录组 | 蛋白组 | 多模态 |
| **准确性（单一模态）** | 高（基因） | 高（蛋白） | - |
| **准确性（综合）** | 中 | 中 | 最高 ⭐ |
| **适用于基因可视化** | ✅ 最佳 | ✅ 可用 | ✅ 推荐 |
| **适用于蛋白可视化** | ✅ 可用 | ✅ 最佳 ⭐ | ✅ 推荐 |
| **适用于细胞分群** | ✅ 可用 | ✅ 可用 | ✅ 最佳 ⭐ |

---

## 📊 修改内容详解

### 1️⃣ 多模态Marker综合可视化

**位置：** 第2042-2081行

**功能：** 为阶段标志基因FeaturePlot添加三种UMAP支持

#### 修改前
```r
# 只有 RNA UMAP 版本
p_feature <- FeaturePlot(seurat_obj, features = key_genes,
                         reduction = "umap", ...)
```

#### 修改后
```r
# --- RNA UMAP版本 ---
p_feature <- FeaturePlot(seurat_obj, features = key_genes,
                         reduction = "umap",
                         cols = c("lightgrey", "red"), ...)
plot_and_save(p_feature, "multimodal_featureplot_stage_markers", ...)

# --- ADT UMAP版本（如果有）⭐ ---
if ("adt_umap" %in% names(seurat_obj@reductions)) {
  p_feature_adt <- FeaturePlot(seurat_obj, features = key_genes,
                               reduction = "adt_umap",
                               cols = c("lightgrey", "red"), ...)
  plot_and_save(p_feature_adt, "multimodal_featureplot_stage_markers_adt", ...)
}

# --- WNN UMAP版本（如果有）---
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  p_feature_wnn <- FeaturePlot(seurat_obj, features = key_genes,
                               reduction = "wnn_umap",
                               cols = c("lightgrey", "red"), ...)
  plot_and_save(p_feature_wnn, "multimodal_featureplot_stage_markers_wnn", ...)
}
```

#### 输出文件（每组3个）

| 文件名 | UMAP空间 | 颜色 | 说明 |
|--------|---------|------|------|
| `multimodal_featureplot_stage_markers.png` | RNA UMAP | 红色 | 基因在RNA空间的表达 |
| `multimodal_featureplot_stage_markers_adt.png` ⭐ | ADT UMAP | 红色 | 基因在蛋白空间的表达 |
| `multimodal_featureplot_stage_markers_wnn.png` | WNN UMAP | 红色 | 基因在整合空间的表达 |

---

### 2️⃣ 扩展基因组可视化 - 基因FeaturePlot

**位置：** 第2207-2292行

**功能：** 为8组功能基因添加三种UMAP支持

#### 8组功能基因

1. 调控基因 (7个)
2. 迁移基因 (8个)
3. 标记基因 (7个)
4. 成熟基因 (5个)
5. TH2反应基因 (7个)
6. 转录因子 (15个)
7. 发育门控 (12个)

**总计：61个基因**

#### 修改后代码模式

```r
for (group_name in names(genes_extended_groups)) {
  genes_to_plot <- head(genes_in_group, 9)

  # --- RNA UMAP版本 ---
  p_feature_group <- FeaturePlot(..., reduction = "umap", ...)
  plot_and_save(p_feature_group, paste0("extended_genes_featureplot_", filename_safe), ...)

  # --- ADT UMAP版本 ⭐ ---
  if ("adt_umap" %in% names(seurat_obj@reductions)) {
    p_feature_group_adt <- FeaturePlot(..., reduction = "adt_umap", ...)
    plot_and_save(p_feature_group_adt, paste0("extended_genes_featureplot_", filename_safe, "_adt"), ...)
  }

  # --- WNN UMAP版本 ---
  if ("wnn_umap" %in% names(seurat_obj@reductions)) {
    p_feature_group_wnn <- FeaturePlot(..., reduction = "wnn_umap", ...)
    plot_and_save(p_feature_group_wnn, paste0("extended_genes_featureplot_", filename_safe, "_wnn"), ...)
  }
}
```

#### 输出文件（每组3个，共21个文件）

**以"调控基因"为例：**

| 文件名 | UMAP空间 | 颜色 |
|--------|---------|------|
| `extended_genes_featureplot_调控基因.png` | RNA UMAP | 红色 |
| `extended_genes_featureplot_调控基因_adt.png` ⭐ | ADT UMAP | 红色 |
| `extended_genes_featureplot_调控基因_wnn.png` | WNN UMAP | 红色 |

**所有8组基因都有3个版本，总计 8 × 3 = 24 个FeaturePlot文件！**

---

### 3️⃣ 扩展基因组可视化 - ADT蛋白FeaturePlot

**位置：** 第2317-2380行

**功能：** 为ADT蛋白FeaturePlot添加三种UMAP支持

#### ADT蛋白列表（12个）

- adt_CD11c, adt_IA-IE, adt_CD135, adt_CD117
- adt_CD172a, adt_CD115, adt_SiglecH, adt_Ly6C
- adt_CD11b, adt_CD43, adt_CD64, adt_Ly6D

#### 修改前

```r
cat("【4/4】生成ADT蛋白分别FeaturePlot（RNA UMAP + WNN UMAP）...\n")
# 只有 RNA UMAP 和 WNN UMAP 版本
```

#### 修改后

```r
cat("【4/4】生成ADT蛋白分别FeaturePlot（RNA UMAP + ADT UMAP + WNN UMAP）...\n")

for (i in seq(1, length(adt_available), by = 9)) {
  adt_batch <- adt_available[i:min(i+8, length(adt_available))]

  # --- RNA UMAP版本 ---
  p_feature_adt_rna <- FeaturePlot(..., reduction = "umap", cols = c("lightgrey", "blue"), ...)
  plot_and_save(p_feature_adt_rna, paste0("adt_proteins_featureplot_batch", ceiling(i/9)), ...)

  # --- ADT UMAP版本 ⭐（最重要！）---
  if ("adt_umap" %in% names(seurat_obj@reductions)) {
    p_feature_adt_adt <- FeaturePlot(..., reduction = "adt_umap", cols = c("lightgrey", "blue"), ...)
    plot_and_save(p_feature_adt_adt, paste0("adt_proteins_featureplot_batch", ceiling(i/9), "_adt"), ...)
  }

  # --- WNN UMAP版本 ---
  if ("wnn_umap" %in% names(seurat_obj@reductions)) {
    p_feature_adt_wnn <- FeaturePlot(..., reduction = "wnn_umap", cols = c("lightgrey", "blue"), ...)
    plot_and_save(p_feature_adt_wnn, paste0("adt_proteins_featureplot_batch", ceiling(i/9), "_wnn"), ...)
  }
}
```

**重要说明：**
- ADT蛋白使用**蓝色渐变** (`lightgrey → blue`)
- 基因使用**红色渐变** (`lightgrey → red`)
- 这样可以一眼区分基因和蛋白可视化

#### 输出文件（每批次3个）

| 文件名 | UMAP空间 | 颜色 | 说明 |
|--------|---------|------|------|
| `adt_proteins_featureplot_batch1.png` | RNA UMAP | 蓝色 | 蛋白在RNA空间的表达 |
| `adt_proteins_featureplot_batch1_adt.png` ⭐⭐⭐ | ADT UMAP | 蓝色 | **蛋白在蛋白空间的表达（最佳）** |
| `adt_proteins_featureplot_batch1_wnn.png` | WNN UMAP | 蓝色 | 蛋白在整合空间的表达 |

---

## 📁 完整输出文件清单

### 标志基因可视化（`results/{mode}/Plots/Markers/`）

#### 1. 多模态标志基因（3个文件）

| 文件名 | UMAP | 颜色 |
|--------|------|------|
| `multimodal_featureplot_stage_markers.png` | RNA | 红 |
| `multimodal_featureplot_stage_markers_adt.png` ⭐ | ADT | 红 |
| `multimodal_featureplot_stage_markers_wnn.png` | WNN | 红 |

#### 2. 扩展基因组 - 8组基因（24个文件）

**每组基因都有3个版本：**

| 基因组 | RNA UMAP | ADT UMAP ⭐ | WNN UMAP |
|--------|---------|-----------|----------|
| 调控基因 | `extended_genes_featureplot_调控基因.png` | `extended_genes_featureplot_调控基因_adt.png` | `extended_genes_featureplot_调控基因_wnn.png` |
| 迁移基因 | `extended_genes_featureplot_迁移基因.png` | `extended_genes_featureplot_迁移基因_adt.png` | `extended_genes_featureplot_迁移基因_wnn.png` |
| 标记基因 | `extended_genes_featureplot_标记基因.png` | `extended_genes_featureplot_标记基因_adt.png` | `extended_genes_featureplot_标记基因_wnn.png` |
| 成熟基因 | `extended_genes_featureplot_成熟基因.png` | `extended_genes_featureplot_成熟基因_adt.png` | `extended_genes_featureplot_成熟基因_wnn.png` |
| TH2反应基因 | `extended_genes_featureplot_TH2反应基因.png` | `extended_genes_featureplot_TH2反应基因_adt.png` | `extended_genes_featureplot_TH2反应基因_wnn.png` |
| 转录因子 | `extended_genes_featureplot_转录因子.png` | `extended_genes_featureplot_转录因子_adt.png` | `extended_genes_featureplot_转录因子_wnn.png` |
| 发育门控 | `extended_genes_featureplot_发育门控.png` | `extended_genes_featureplot_发育门控_adt.png` | `extended_genes_featureplot_发育门控_wnn.png` |

#### 3. 扩展基因组 - ADT蛋白（6个文件）

| 批次 | RNA UMAP | ADT UMAP ⭐⭐⭐ | WNN UMAP |
|------|---------|--------------|----------|
| Batch 1 | `adt_proteins_featureplot_batch1.png` | `adt_proteins_featureplot_batch1_adt.png` | `adt_proteins_featureplot_batch1_wnn.png` |
| Batch 2 | `adt_proteins_featureplot_batch2.png` | `adt_proteins_featureplot_batch2_adt.png` | `adt_proteins_featureplot_batch2_wnn.png` |

**总计：33个FeaturePlot文件！**

---

## 🎨 颜色方案

| 可视化类型 | RNA UMAP | ADT UMAP | WNN UMAP | 颜色渐变 |
|-----------|---------|---------|---------|---------|
| **基因表达** | ✅ | ✅ | ✅ | lightgrey → **red** 🔴 |
| **ADT蛋白** | ✅ | ✅ | ✅ | lightgrey → **blue** 🔵 |

**设计理念：**
- 红色代表基因（RNA）
- 蓝色代表蛋白（ADT）
- 一眼区分可视化类型

---

## 🔍 关键技术特性

### 1. 条件检测降维空间

```r
# RNA UMAP（必定存在）
p_rna <- FeaturePlot(seurat_obj, reduction = "umap", ...)

# ADT UMAP（条件检测）
if ("adt_umap" %in% names(seurat_obj@reductions)) {
  p_adt <- FeaturePlot(seurat_obj, reduction = "adt_umap", ...)
}

# WNN UMAP（条件检测）
if ("wnn_umap" %in% names(seurat_obj@reductions)) {
  p_wnn <- FeaturePlot(seurat_obj, reduction = "wnn_umap", ...)
}
```

**优势：**
- 兼容有/无ADT数据的Seurat对象
- 兼容有/无WNN分析的对象
- 不会因缺少降维结果而报错

---

### 2. 文件命名规范

| 原文件名 | ADT版本 | WNN版本 |
|---------|---------|---------|
| `xxx.png` | `xxx_adt.png` ⭐ | `xxx_wnn.png` |
| `xxx_batch1.png` | `xxx_batch1_adt.png` ⭐ | `xxx_batch1_wnn.png` |

**优势：**
- 清晰标识降维空间类型
- 方便批量查找和对比
- 保持文件组织一致性

---

### 3. 标题和副标题标注

**RNA UMAP：**
```r
title = "基因表达 (RNA UMAP)"
subtitle = "Expression level in RNA space"
```

**ADT UMAP：**
```r
title = "基因表达 (ADT UMAP)"
subtitle = "Expression level in ADT protein space"
```

**WNN UMAP：**
```r
title = "基因表达 (WNN UMAP)"
subtitle = "Expression level in WNN-integrated space"
```

**ADT蛋白在ADT UMAP特别标注：**
```r
title = "ADT蛋白表达（发育门控） - ADT UMAP ⭐"
subtitle = "批次1 | 9个蛋白 | ADT蛋白降维空间（最佳）"
```

---

## 📚 使用场景推荐

### 场景 1: 基因表达可视化

**推荐优先级：**
1. **WNN UMAP** ⭐⭐⭐ - 综合RNA+蛋白信息，最准确
2. **RNA UMAP** ⭐⭐ - 反映转录组层面
3. **ADT UMAP** ⭐ - 参考对比

**示例：** 查看转录因子Irf8在cDC1发育中的表达模式
- WNN UMAP: 最准确地显示Irf8在各发育阶段的表达
- RNA UMAP: 显示转录组水平的Irf8表达
- ADT UMAP: 看Irf8表达是否与蛋白表型相关

---

### 场景 2: ADT蛋白表达可视化 ⭐⭐⭐

**推荐优先级：**
1. **ADT UMAP** ⭐⭐⭐ - **蛋白在蛋白空间，最准确！**
2. **WNN UMAP** ⭐⭐ - 综合信息，也很好
3. **RNA UMAP** ⭐ - 参考对比

**示例：** 查看CD172a（cDC2标志蛋白）的表达模式
- **ADT UMAP**: 最准确显示CD172a蛋白表达在cDC2上高，cDC1上低
- WNN UMAP: 综合RNA+蛋白信息显示CD172a
- RNA UMAP: 可能不够准确（因为是蛋白，不是RNA）

---

### 场景 3: 对比不同降维空间

**研究问题：** RNA表达和蛋白表达是否一致？

**操作：** 同时查看同一个基因/蛋白在三种UMAP上的分布

**示例：** 查看H2-Ab1（MHC-II）
- RNA UMAP: 显示H2-Ab1 mRNA表达
- ADT UMAP: 显示IA-IE（MHC-II）蛋白表达
- 对比: 发现某些细胞mRNA高但蛋白低（可能是转录后调控）

---

## ⚠️ 注意事项

### 1. 降维结果的可用性

| UMAP类型 | 何时存在 |
|---------|---------|
| RNA UMAP (`umap`) | ✅ 总是存在（Seurat标准流程） |
| ADT UMAP (`adt_umap`) | ⚠️ 仅CITE-seq数据且运行了ADT降维 |
| WNN UMAP (`wnn_umap`) | ⚠️ 仅CITE-seq数据且运行了WNN分析 |

**检查方法：**
```r
# 查看所有降维结果
names(seurat_obj@reductions)
# [1] "pca"      "umap"     "adt_pca"  "adt_umap" "wnn_umap"

# 检查ADT UMAP
"adt_umap" %in% names(seurat_obj@reductions)  # TRUE or FALSE

# 检查WNN UMAP
"wnn_umap" %in% names(seurat_obj@reductions)  # TRUE or FALSE
```

---

### 2. 输入文件要求

**必须使用 `*_final.rds` 才有完整的降维结果：**

| 文件名 | RNA UMAP | ADT UMAP | WNN UMAP |
|--------|---------|---------|---------|
| `cDC1_final.rds` ✅ | ✅ | ✅ | ✅ |
| `cDC1_pre_wnn.rds` ⚠️ | ✅ | ✅ | ❌ |

脚本已自动优先加载 `*_final.rds`（第221-226行）。

---

### 3. 基因vs蛋白可视化的最佳实践

#### 查看基因表达：
```r
# 推荐顺序：WNN UMAP > RNA UMAP > ADT UMAP
FeaturePlot(seurat_obj, features = "Irf8", reduction = "wnn_umap")  # 最佳
FeaturePlot(seurat_obj, features = "Irf8", reduction = "umap")      # 很好
FeaturePlot(seurat_obj, features = "Irf8", reduction = "adt_umap")  # 参考
```

#### 查看蛋白表达：
```r
# 推荐顺序：ADT UMAP > WNN UMAP > RNA UMAP
FeaturePlot(seurat_obj, features = "adt_CD172a", reduction = "adt_umap")  # 最佳 ⭐⭐⭐
FeaturePlot(seurat_obj, features = "adt_CD172a", reduction = "wnn_umap")  # 很好
FeaturePlot(seurat_obj, features = "adt_CD172a", reduction = "umap")      # 参考
```

---

## 🎯 核心改进总结

### 修改统计

| 章节 | 行号 | 修改内容 | 新增UMAP支持 |
|------|------|---------|-------------|
| 多模态Marker可视化 | 2042-2081 | 阶段标志基因FeaturePlot | ✅ RNA + ADT ⭐ + WNN |
| 扩展基因组可视化（基因） | 2207-2292 | 8组基因FeaturePlot | ✅ RNA + ADT ⭐ + WNN |
| 扩展基因组可视化（蛋白） | 2317-2380 | ADT蛋白FeaturePlot | ✅ RNA + ADT ⭐⭐⭐ + WNN |

### 新增文件统计

| 可视化类型 | 原有文件 | 新增ADT版本 | 新增WNN版本 | 总文件数 |
|-----------|---------|-----------|-----------|---------|
| 多模态标志基因 | 1 | +1 ⭐ | +1 | 3 |
| 8组功能基因 | 8 | +8 ⭐ | +8 | 24 |
| ADT蛋白（2批次） | 2 | +2 ⭐⭐⭐ | +2 | 6 |
| **总计** | **11** | **+11** | **+11** | **33** |

---

## ✨ 总结

### 核心成果

✅ **完整支持三种UMAP空间** - RNA、ADT、WNN
✅ **33个FeaturePlot文件** - 覆盖所有基因和蛋白
✅ **自动条件检测** - 兼容不同数据集
✅ **清晰文件命名** - `_adt` 和 `_wnn` 后缀
✅ **颜色区分** - 红色基因、蓝色蛋白
✅ **最佳实践标注** - ADT蛋白在ADT UMAP标记为"最佳"

### 关键优势

1. **全面性** - 每个FeaturePlot都有3种UMAP版本可选
2. **准确性** - 为不同可视化类型推荐最佳UMAP空间
3. **灵活性** - 用户可以对比不同UMAP空间的结果
4. **一致性** - 统一的命名规范和颜色方案
5. **可靠性** - 条件检测确保不会因缺少降维结果而报错

### 使用价值

**对于研究者：**
- 可以从RNA、蛋白、整合三个角度查看数据
- 发现RNA-蛋白表达的一致性或差异
- 选择最适合的降维空间进行可视化

**对于CITE-seq数据：**
- 充分利用多模态数据的优势
- ADT UMAP为蛋白可视化提供最佳视角 ⭐
- WNN UMAP为综合分析提供最准确的细胞分群

---

## 📖 相关文档

- `下游分析_WNN_UMAP_完整支持总结.md` - 之前的WNN UMAP支持总结
- `下游分析准备模块_详细说明.md` - 下游分析数据准备说明
- `WNN_UMAP_作者注释可视化说明.md` - WNN UMAP作者注释功能

---

**现在，用户可以在RNA、ADT、WNN三种UMAP空间中自由查看所有基因和蛋白的表达模式！** 🎉
