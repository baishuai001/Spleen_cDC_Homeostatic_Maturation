# CITE-seq 分析框架设计文档

## 📊 问题分析

### 需要支持的分析场景

1. **cDC1 单独分析** (RNA + ADT)
2. **cDC2 单独分析** (RNA + ADT)
3. **cDC1 + cDC2 整合分析** (RNA + ADT)

### 数据来源

根据 [GSE228544](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228544):

| 文件 | 包含细胞类型 |
|------|-------------|
| `GSE228544_seuratObj_paper_dietSAM2and3_WT_subset_2023.rds.gz` | cDC1 各亚型 |
| `GSE228544_seuratObj_paper_dietSAM2and3_WT_subset2_v2_2023.rds.gz` | cDC2 各亚型 |

---

## 🎯 解决方案对比

### 方案A: 分离脚本方案

```
scripts/
├── 01_cDC1_analysis.Rmd
├── 02_cDC2_analysis.Rmd
└── 03_integrated_analysis.Rmd
```

**优点:**
- 每个脚本独立运行,互不干扰
- 适合不同分析者分工

**缺点:**
- ❌ 代码高度重复 (~80% 相同代码)
- ❌ 质控参数难以统一
- ❌ 维护成本高 (改一处要改三处)
- ❌ 缓存无法共享
- ❌ 图表风格容易不一致

---

### 方案B: 单一参数化脚本 (推荐 ⭐)

```
scripts/
└── 01_CITE-seq_unified_analysis.Rmd
    ├── 参数配置区 (analysis_mode 切换)
    ├── 通用质控模块
    ├── 通用归一化模块
    ├── 通用聚类模块
    └── 模式特异性分析
```

**优点:**
- ✅ 代码复用率 95%+
- ✅ 统一的质控标准
- ✅ 缓存智能共享
- ✅ 易于维护和更新
- ✅ 图表风格一致
- ✅ 支持快速切换模式

**缺点:**
- 需要更好的代码组织
- 初期设计稍复杂 (但我已经帮你设计好了!)

---

## 🏗️ 推荐架构: 参数化统一框架

### 核心设计思想

```r
# 一个参数控制全局行为
analysis_mode <- "cDC1"  # 可选: "cDC1", "cDC2", "integrated"

# 所有下游代码根据 mode 自动调整
if (analysis_mode == "cDC1") {
    # 加载 cDC1 特异性参数
} else if (analysis_mode == "cDC2") {
    # 加载 cDC2 特异性参数
} else {
    # 整合分析参数
}
```

---

## 📐 详细实现方案

### 1. 参数配置系统

```r
################################################################################
##################### 分析模式配置 #############################################
################################################################################

# ==================== 核心参数: 选择分析模式 ====================
# 可选值:
#   "cDC1"       - 仅分析 cDC1 细胞
#   "cDC2"       - 仅分析 cDC2 细胞
#   "integrated" - cDC1 + cDC2 整合分析
analysis_mode <- "cDC1"

# ==================== 自动配置参数 ====================
config <- get_analysis_config(analysis_mode)

# config 包含:
# - source_file: 源数据文件名
# - cell_types: 要提取的细胞类型列表
# - dims_to_use: PCA 维度
# - resolution: 聚类分辨率
# - project_name: 项目名称
# - cache_prefix: 缓存文件前缀
```

### 2. 配置函数实现

```r
get_analysis_config <- function(mode) {
  configs <- list(
    cDC1 = list(
      source_file = "GSE228544_seuratObj_paper_dietSAM2and3_WT_subset_2023.rds.gz",
      cell_types = c("pre-cDC1s", "Proliferating cDC1s", "Early immature cDC1s",
                     "Late immature cDC1s", "Early mature cDC1s", "Late mature cDC1s"),
      dims_rna = 25,
      dims_adt = 15,
      resolution = 0.8,
      project_name = "GSE228544_cDC1",
      cache_prefix = "cDC1",
      colors = RColorBrewer::brewer.pal(6, "Set2")
    ),

    cDC2 = list(
      source_file = "GSE228544_seuratObj_paper_dietSAM2and3_WT_subset2_v2_2023.rds.gz",
      cell_types = c("pre-cDC2s", "Proliferating cDC2s", "Early immature cDC2s",
                     "Late immature cDC2s", "Early mature cDC2s", "Late mature cDC2s"),
      dims_rna = 30,
      dims_adt = 20,
      resolution = 0.8,
      project_name = "GSE228544_cDC2",
      cache_prefix = "cDC2",
      colors = RColorBrewer::brewer.pal(6, "Set1")
    ),

    integrated = list(
      source_files = c(
        "GSE228544_seuratObj_paper_dietSAM2and3_WT_subset_2023.rds.gz",
        "GSE228544_seuratObj_paper_dietSAM2and3_WT_subset2_v2_2023.rds.gz"
      ),
      cell_types_list = list(
        cDC1 = c("pre-cDC1s", "Proliferating cDC1s", "Early immature cDC1s",
                 "Late immature cDC1s", "Early mature cDC1s", "Late mature cDC1s"),
        cDC2 = c("pre-cDC2s", "Proliferating cDC2s", "Early immature cDC2s",
                 "Late immature cDC2s", "Early mature cDC2s", "Late mature cDC2s")
      ),
      dims_rna = 25,
      dims_adt = 20,
      resolution = 0.8,
      project_name = "GSE228544_integrated",
      cache_prefix = "integrated",
      batch_correction = TRUE,  # 整合分析需要批次校正
      integration_method = "harmony"  # 或 "CCA", "RPCA"
    )
  )

  if (!mode %in% names(configs)) {
    stop("Invalid analysis_mode. Choose from: ", paste(names(configs), collapse = ", "))
  }

  return(configs[[mode]])
}
```

### 3. 智能缓存系统升级

```r
# 缓存文件名包含模式前缀,避免冲突
cache_file <- paste0(config$cache_prefix, "_01_subset_extracted.rds")

# 示例:
# cDC1 模式: "cDC1_01_subset_extracted.rds"
# cDC2 模式: "cDC2_01_subset_extracted.rds"
# 整合模式: "integrated_01_subset_extracted.rds"
```

### 4. 数据提取模块 (模式自适应)

```r
# cDC1/cDC2 单独分析
if (config$cache_prefix %in% c("cDC1", "cDC2")) {

  subset_data <- load_or_compute(
    cache_file = paste0(config$cache_prefix, "_01_subset_extracted.rds"),
    description = paste0(config$cache_prefix, " 子集提取"),
    compute_fn = function() {
      # 解压和读取
      full_path <- file.path(data_dir, config$source_file)
      try(gunzip(full_path, remove = FALSE, overwrite = TRUE))

      rds_path <- gsub("\\.gz$", "", full_path)
      seuratObjAll <- readRDS(rds_path)

      # 提取特定细胞类型
      existing_types <- intersect(config$cell_types, levels(Idents(seuratObjAll)))
      cells_to_keep <- WhichCells(seuratObjAll, idents = existing_types)

      # 过滤数据
      rawDataRNA_filtered <- rawDataRNA[, cells_to_keep]
      rawDataADT_filtered <- rawDataADT[, cells_to_keep]

      return(list(
        rawDataRNA = rawDataRNA_filtered,
        rawDataADT = rawDataADT_filtered,
        cells = cells_to_keep,
        metadata = seuratObjAll@meta.data[cells_to_keep, ]
      ))
    }
  )

} else if (config$cache_prefix == "integrated") {
  # 整合分析: 提取 cDC1 + cDC2

  subset_data <- load_or_compute(
    cache_file = "integrated_01_subset_extracted.rds",
    description = "cDC1 + cDC2 整合数据提取",
    compute_fn = function() {
      all_cells <- list()
      all_metadata <- list()

      # 循环提取 cDC1 和 cDC2
      for (i in seq_along(config$source_files)) {
        subset_name <- names(config$cell_types_list)[i]
        file_path <- file.path(data_dir, config$source_files[i])

        # 解压读取
        try(gunzip(file_path, remove = FALSE, overwrite = TRUE))
        rds_path <- gsub("\\.gz$", "", file_path)
        seurat_obj <- readRDS(rds_path)

        # 提取细胞
        cell_types <- config$cell_types_list[[i]]
        existing_types <- intersect(cell_types, levels(Idents(seurat_obj)))
        cells <- WhichCells(seurat_obj, idents = existing_types)

        all_cells[[subset_name]] <- cells

        # 添加批次信息
        metadata <- seurat_obj@meta.data[cells, ]
        metadata$subset <- subset_name
        all_metadata[[subset_name]] <- metadata
      }

      # 合并细胞
      combined_cells <- unlist(all_cells)
      combined_metadata <- do.call(rbind, all_metadata)

      # 过滤 RNA 和 ADT
      rawDataRNA_filtered <- rawDataRNA[, combined_cells]
      rawDataADT_filtered <- rawDataADT[, combined_cells]

      return(list(
        rawDataRNA = rawDataRNA_filtered,
        rawDataADT = rawDataADT_filtered,
        cells = combined_cells,
        metadata = combined_metadata,
        batch_info = combined_metadata$subset
      ))
    }
  )
}
```

### 5. 整合分析特殊处理

```r
# 仅在整合模式下执行批次校正
if (config$cache_prefix == "integrated" && config$batch_correction) {

  # 使用 Harmony 进行批次校正
  if (config$integration_method == "harmony") {
    library(harmony)

    seurat_obj <- load_or_compute(
      cache_file = "integrated_04_harmony_integrated.rds",
      description = "Harmony 批次校正",
      compute_fn = function() {
        # 在 PCA 后运行 Harmony
        seurat_obj <- RunHarmony(
          seurat_obj,
          group.by.vars = "subset",  # 以 cDC1/cDC2 为批次
          reduction = "pca",
          dims.use = 1:config$dims_rna
        )

        # 使用 harmony reduction 进行下游分析
        seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:config$dims_rna)
        seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:config$dims_rna)
        seurat_obj <- FindClusters(seurat_obj, resolution = config$resolution)

        return(seurat_obj)
      }
    )
  }
}
```

---

## 📊 ADT 整合策略

### 策略 1: WNN (Weighted Nearest Neighbor) - 推荐

**适用场景:** 所有模式 (cDC1, cDC2, integrated)

```r
# WNN 是模式无关的,对所有分析方案都适用
seurat_obj <- FindMultiModalNeighbors(
  seurat_obj,
  reduction.list = list("pca", "adt_pca"),  # 或 "harmony" (整合模式)
  dims.list = list(1:config$dims_rna, 1:config$dims_adt)
)

seurat_obj <- RunUMAP(
  seurat_obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn_umap"
)
```

### 策略 2: 分层整合

**适用场景:** 整合分析

```r
# 先整合 RNA (批次校正)
seurat_obj <- RunHarmony(..., reduction = "pca")

# 再整合 ADT (批次校正)
seurat_obj <- RunHarmony(..., reduction = "adt_pca")

# 最后用 WNN 融合
seurat_obj <- FindMultiModalNeighbors(
  reduction.list = list("harmony", "adt_harmony")
)
```

---

## 🔄 模式切换方法

### 方法 1: 修改参数重新运行 (最简单)

```r
# 1. 修改 analysis_mode
analysis_mode <- "cDC2"  # 从 "cDC1" 改为 "cDC2"

# 2. 重新运行整个脚本
# 系统会自动:
# - 加载 cDC2 配置
# - 使用 cDC2 缓存文件
# - 应用 cDC2 参数
```

### 方法 2: 命令行参数 (高级)

```r
# 在 Rmd 开头添加
args <- commandArgs(trailingOnly = TRUE)
analysis_mode <- ifelse(length(args) > 0, args[1], "cDC1")
```

```bash
# 运行不同模式
Rscript -e "rmarkdown::render('analysis.Rmd', params = list(mode = 'cDC1'))"
Rscript -e "rmarkdown::render('analysis.Rmd', params = list(mode = 'cDC2'))"
Rscript -e "rmarkdown::render('analysis.Rmd', params = list(mode = 'integrated'))"
```

### 方法 3: RMarkdown 参数化 (推荐用于批量生成报告)

```yaml
---
title: "CITE-seq Analysis"
params:
  mode: "cDC1"
---
```

```r
# 在代码块中使用
analysis_mode <- params$mode
```

---

## 📁 推荐的项目结构

```
Spleen_cDC_Homeostatic_Maturation/
├── data/
│   └── GSE228544/
│       ├── GSE228544_seuratObj_..._subset_2023.rds.gz    # cDC1
│       └── GSE228544_seuratObj_..._subset2_v2_2023.rds.gz # cDC2
├── scripts/
│   └── 01_CITE-seq_unified_analysis.Rmd  # 统一脚本
├── results/
│   ├── Cache/
│   │   ├── cDC1_01_subset_extracted.rds
│   │   ├── cDC1_04_seurat_clustered.rds
│   │   ├── cDC2_01_subset_extracted.rds
│   │   ├── cDC2_04_seurat_clustered.rds
│   │   ├── integrated_01_subset_extracted.rds
│   │   └── integrated_04_harmony_integrated.rds
│   ├── Reports/
│   │   ├── cDC1_analysis.html
│   │   ├── cDC2_analysis.html
│   │   └── integrated_analysis.html
│   └── Plots/
│       ├── cDC1/
│       ├── cDC2/
│       └── integrated/
└── docs/
    ├── ANALYSIS_FRAMEWORK_DESIGN.md  # 本文档
    └── OPTIMIZATION_README.md
```

---

## 🎯 最佳实践建议

### 1. 开发流程

```
第一步: 用 cDC1 模式开发和调试完整流程
  ↓
第二步: 切换到 cDC2 模式验证代码通用性
  ↓
第三步: 切换到 integrated 模式测试整合分析
  ↓
第四步: 批量生成所有模式的报告
```

### 2. 参数调优

| 参数 | cDC1 | cDC2 | Integrated |
|------|------|------|------------|
| dims_rna | 20-25 | 25-30 | 25-30 |
| dims_adt | 10-15 | 15-20 | 15-20 |
| resolution | 0.6-0.8 | 0.6-0.8 | 0.4-0.6 |
| integration_method | - | - | harmony/CCA |

### 3. 质控标准统一

```r
# 所有模式使用相同的质控参数
qc_params <- list(
  nmad_low_feature = 5,
  nmad_high_feature = 5,
  nmad_low_UMI = 5,
  nmad_high_UMI = 5,
  nmad_high_mito = 10
)

# 确保 cDC1, cDC2, integrated 使用相同标准
```

### 4. 缓存管理策略

```r
# 查看特定模式的缓存
list_mode_cache <- function(mode) {
  list.files("results/Cache/", pattern = paste0("^", mode, "_"), full.names = TRUE)
}

# 清除特定模式的缓存
clear_mode_cache <- function(mode) {
  files <- list_mode_cache(mode)
  file.remove(files)
  cat("Removed", length(files), "cache files for mode:", mode, "\n")
}

# 用法
clear_mode_cache("cDC2")  # 仅清除 cDC2 的缓存
```

---

## 🔬 科学分析建议

### 单独分析的价值

1. **cDC1 单独分析**
   - 深入理解 cDC1 发育轨迹
   - 识别 cDC1 特异性标志基因
   - 细粒度的亚群分类

2. **cDC2 单独分析**
   - 探索 cDC2 成熟过程
   - 发现 cDC2 功能特征
   - 与 cDC1 对比分析

### 整合分析的价值

1. **共同发育特征**
   - 识别 cDC1 和 cDC2 的共同祖先标志
   - 比较成熟过程的异同

2. **差异基因分析**
   - cDC1 vs cDC2 的直接对比
   - 亚型特异性通路分析

3. **空间关系**
   - 在统一的 UMAP 空间中可视化两种细胞
   - 理解它们的发育分叉点

---

## 📊 结论

**推荐方案: 单一参数化脚本**

| 方面 | 评分 | 说明 |
|------|------|------|
| 代码复用 | ⭐⭐⭐⭐⭐ | 95%+ 代码共享 |
| 维护性 | ⭐⭐⭐⭐⭐ | 改一处,全部更新 |
| 可扩展性 | ⭐⭐⭐⭐⭐ | 轻松添加新模式 |
| 学习曲线 | ⭐⭐⭐⭐ | 初期稍复杂,但值得 |
| 性能 | ⭐⭐⭐⭐⭐ | 智能缓存,模式间复用 |

**下一步行动:**
1. 实现参数化统一脚本
2. 测试三种模式
3. 生成对比报告

---

**文档版本:** v1.0
**创建日期:** 2026-01-06
**作者:** Claude Code
