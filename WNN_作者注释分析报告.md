# WNN 分析后添加作者注释的可行性分析

## 📋 问题概述

您的问题是：**能否将WNN分析后的聚类结果，加上作者的原始注释？**

**答案：✅ 完全可以，而且数据已经准备好了，只需要修改代码即可。**

---

## 🔍 当前代码流程分析

### 1️⃣ Section 2.2: 提取原始数据和作者注释

**文件位置：** `CITE-seq_Complete_Analysis_v3.Rmd` 第663-717行

**功能：** 从原始 Seurat 对象中提取数据

```r
subset_data <- load_or_compute(
  cache_file = paste0(config$cache_prefix, "_01_subset_extracted.rds"),
  compute_fn = function() {
    # 读取原始 Seurat 对象
    subset_seurat <- readRDS(rds_path)

    # 保存 metadata（包含作者的所有注释！）
    return(list(
      rawDataRNA = rawDataRNA_filtered,
      rawDataADT = rawDataADT_filtered,
      cells = cells_to_keep,
      metadata = subset_seurat@meta.data[cells_to_keep, ],  # ✅ 包含作者注释
      original_idents = Idents(subset_seurat)[cells_to_keep],
      cell_types = all_idents
    ))
  }
)
```

**关键点：**
- ✅ `metadata` 字段包含了原始 Seurat 对象的所有 metadata
- ✅ 原始 Seurat 对象（如 `Spleen.cDC1_subset_final.rds`）包含作者的注释字段：
  - `final_annotation2021`（最终版注释）
  - `Annotation_detailed`（详细版注释）
  - `Annotation`（基础版注释）
  - 等等

---

### 2️⃣ 第819行: 提取到全局变量

```r
metadata_original <- subset_data$metadata
```

**说明：**
- ✅ `metadata_original` 包含了所有作者注释
- ✅ 这个变量在整个脚本中都可用

---

### 3️⃣ Section 5: 创建 Seurat 对象（**问题所在！**）

**文件位置：** 第1252-1335行

**当前代码：**

```r
seurat_obj <- load_or_compute(
  cache_file = paste0(config$cache_prefix, "_04_seurat_clustered.rds"),
  compute_fn = function() {
    # 创建 Seurat 对象
    seurat_obj <- CreateSeuratObject(counts = rawDataFiltered, project = config$project_name)

    # ⚠️ 问题：只添加了 original_cluster
    if (exists("metadata_original", envir = .GlobalEnv) && !is.null(metadata_original)) {
      if ("original_cluster" %in% colnames(metadata_original)) {
        seurat_obj$original_cluster <- metadata_original$original_cluster[colnames(seurat_obj)]
      }
    }

    # ❌ 缺失：没有添加其他作者注释字段！
    # 如 final_annotation2021, Annotation_detailed 等

    # ... 后续的 SCTransform, UMAP, 聚类等
    return(seurat_obj)
  }
)
```

**问题：**
- ❌ 只检查和添加了 `original_cluster` 字段
- ❌ 没有添加作者的其他重要注释字段（`final_annotation2021`、`Annotation_detailed` 等）

---

### 4️⃣ Section 8: WNN 分析

**文件位置：** 第2302-2397行

```r
seurat_obj_adt_rna_wnn <- load_or_compute(
  cache_file = paste0(config$cache_prefix, "_07_wnn_integrated.rds"),
  compute_fn = function() {
    # ... WNN 融合
    seurat_obj_clean <- FindMultiModalNeighbors(...)
    seurat_obj_clean <- RunUMAP(..., reduction.name = "wnn_umap")
    seurat_obj_clean <- FindClusters(...)
    seurat_obj_clean$WNN_clusters <- Idents(seurat_obj_clean)

    return(seurat_obj_clean)
  }
)
```

**说明：**
- ✅ WNN 对象 `seurat_obj_adt_rna_wnn` 继承了前面 Seurat 对象的所有 metadata
- ❌ 但由于 Section 5 没有添加作者注释，所以 WNN 对象也没有

---

### 5️⃣ 最终保存

**文件位置：** 第3741行

```r
final_rds <- paste0(output.dir, "Robjects/", config$project_name, "_final.rds")
saveRDS(seurat_obj_adt_rna_wnn, final_rds)
```

---

## ✅ 解决方案

### 方案一：修改 Section 5，在创建 Seurat 对象时添加所有作者注释（**推荐**）

**修改位置：** `CITE-seq_Complete_Analysis_v3.Rmd` 第1265-1270行

**原始代码：**

```r
# 添加原始细胞类型信息
if (exists("metadata_original", envir = .GlobalEnv) && !is.null(metadata_original)) {
  if ("original_cluster" %in% colnames(metadata_original)) {
    seurat_obj$original_cluster <- metadata_original$original_cluster[colnames(seurat_obj)]
  }
}
```

**修改后的代码：**

```r
# 添加原始作者注释信息（完整版）
if (exists("metadata_original", envir = .GlobalEnv) && !is.null(metadata_original)) {

  # 定义需要添加的作者注释字段
  annotation_fields <- c(
    "final_annotation2021",      # 最终版注释（优先）
    "Annotation_detailed",        # 详细版注释
    "Annotation",                 # 基础版注释
    "annotated_clusters",         # 注释后的聚类
    "minimalist_annotation",      # 简化版注释
    "original_cluster"            # 原始聚类
  )

  # 检查哪些字段在 metadata_original 中存在
  available_fields <- intersect(annotation_fields, colnames(metadata_original))

  if (length(available_fields) > 0) {
    cat("  添加作者注释字段:\n")

    # 添加所有可用的注释字段
    for (field in available_fields) {
      seurat_obj[[field]] <- metadata_original[[field]][colnames(seurat_obj)]
      cat("    ✓", field, "\n")
    }

    # 创建标准化的注释字段（用于后续分析）
    if ("final_annotation2021" %in% available_fields) {
      seurat_obj$author_annotation <- metadata_original$final_annotation2021[colnames(seurat_obj)]
      cat("  ✓ 使用 final_annotation2021 作为主要注释 (author_annotation)\n")
    } else if ("Annotation" %in% available_fields) {
      seurat_obj$author_annotation <- metadata_original$Annotation[colnames(seurat_obj)]
      cat("  ✓ 使用 Annotation 作为主要注释 (author_annotation)\n")
    }

    if ("Annotation_detailed" %in% available_fields) {
      seurat_obj$author_annotation_detailed <- metadata_original$Annotation_detailed[colnames(seurat_obj)]
      cat("  ✓ 添加详细版注释 (author_annotation_detailed)\n")
    }
  } else {
    cat("  ⚠️ 未在 metadata_original 中找到作者注释字段\n")
  }
}
```

**优势：**
- ✅ 一次性添加所有作者注释字段
- ✅ 自动检测可用字段，适应不同的数据源
- ✅ 创建标准化的 `author_annotation` 和 `author_annotation_detailed` 字段
- ✅ WNN 对象会自动继承这些字段

---

### 方案二：在 WNN 分析后添加注释（备选方案）

如果不想修改主脚本，可以在 WNN 分析完成后再添加：

**在 Section 8 WNN 分析的 compute_fn 中添加：**

```r
seurat_obj_adt_rna_wnn <- load_or_compute(
  cache_file = paste0(config$cache_prefix, "_07_wnn_integrated.rds"),
  compute_fn = function() {
    # ... 原有的 WNN 融合代码 ...

    # ========== 新增：添加作者注释 ==========
    if (exists("metadata_original", envir = .GlobalEnv) && !is.null(metadata_original)) {

      # 确保细胞顺序一致
      common_cells <- intersect(colnames(seurat_obj_clean), rownames(metadata_original))

      # 添加作者注释
      if ("final_annotation2021" %in% colnames(metadata_original)) {
        seurat_obj_clean$author_annotation <- metadata_original[common_cells, "final_annotation2021"]
        cat("  ✓ 添加 author_annotation (final_annotation2021)\n")
      }

      if ("Annotation_detailed" %in% colnames(metadata_original)) {
        seurat_obj_clean$author_annotation_detailed <- metadata_original[common_cells, "Annotation_detailed"]
        cat("  ✓ 添加 author_annotation_detailed\n")
      }
    }

    return(seurat_obj_clean)
  }
)
```

---

### 方案三：事后补救（如果已经运行完）

如果您已经运行完 WNN 分析，可以事后添加注释：

```r
# 读取已保存的 WNN 对象
seurat_wnn <- readRDS("results/cDC1/Robjects/cDC1_final.rds")

# 读取原始注释
subset_data <- readRDS("results/Cache/cDC1_01_subset_extracted.rds")
metadata_original <- subset_data$metadata

# 添加注释
common_cells <- intersect(colnames(seurat_wnn), rownames(metadata_original))

if ("final_annotation2021" %in% colnames(metadata_original)) {
  seurat_wnn$author_annotation <- metadata_original[common_cells, "final_annotation2021"]
}

if ("Annotation_detailed" %in% colnames(metadata_original)) {
  seurat_wnn$author_annotation_detailed <- metadata_original[common_cells, "Annotation_detailed"]
}

# 重新保存
saveRDS(seurat_wnn, "results/cDC1/Robjects/cDC1_final_with_annotations.rds")
```

---

## 🎯 推荐实施步骤

### 第1步：验证 metadata_original 中确实有作者注释

```r
# 在 RStudio 中运行
subset_data <- readRDS("results/Cache/cDC1_01_subset_extracted.rds")
metadata_original <- subset_data$metadata

# 查看包含哪些列
colnames(metadata_original)

# 查看作者注释的分布
table(metadata_original$final_annotation2021)
table(metadata_original$Annotation_detailed)
```

### 第2步：选择合适的方案修改脚本

- **如果还没运行 WNN 分析**：使用方案一（修改 Section 5）
- **如果已经运行完**：使用方案三（事后补救）

### 第3步：验证结果

```r
# 读取 WNN 对象
seurat_wnn <- readRDS("results/cDC1/Robjects/cDC1_final.rds")

# 检查是否包含作者注释
"author_annotation" %in% colnames(seurat_wnn@meta.data)
table(seurat_wnn$author_annotation)

# 可视化对比
DimPlot(seurat_wnn, reduction = "wnn_umap", group.by = "WNN_clusters")
DimPlot(seurat_wnn, reduction = "wnn_umap", group.by = "author_annotation")
```

---

## 📊 预期效果

添加作者注释后，您可以：

1. **对比 WNN 聚类 vs 作者注释**
   ```r
   # UMAP 可视化对比
   p1 <- DimPlot(seurat_wnn, reduction = "wnn_umap", group.by = "WNN_clusters")
   p2 <- DimPlot(seurat_wnn, reduction = "wnn_umap", group.by = "author_annotation")
   p1 | p2
   ```

2. **分析对应关系**
   ```r
   # 交叉表
   table(seurat_wnn$WNN_clusters, seurat_wnn$author_annotation)

   # 热图
   library(pheatmap)
   correspondence <- table(seurat_wnn$WNN_clusters, seurat_wnn$author_annotation)
   pheatmap(correspondence, scale = "row")
   ```

3. **验证 marker 基因**
   ```r
   # 基于作者注释找 marker
   Idents(seurat_wnn) <- seurat_wnn$author_annotation
   markers <- FindAllMarkers(seurat_wnn, only.pos = TRUE)
   ```

---

## 🎉 结论

**完全可以将 WNN 分析后的聚类结果加上作者的原始注释！**

- ✅ 数据已经在 `_01_subset_extracted.rds` 中准备好
- ✅ 只需要修改代码，将 `metadata_original` 中的注释字段添加到 Seurat 对象
- ✅ 三种方案任选其一，推荐使用方案一（修改 Section 5）

---

**需要我帮您实施哪个方案吗？**
