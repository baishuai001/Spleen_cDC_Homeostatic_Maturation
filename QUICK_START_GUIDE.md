# 🚀 快速使用指南

## 三种分析模式快速切换

### 方式 1: 修改参数 (最简单⭐)

打开 `01_CITE-seq_unified_analysis.Rmd`,找到第 30 行左右:

```r
# ==================== 在这里修改分析模式 ====================
analysis_mode <- "cDC1"  # 可选: "cDC1", "cDC2", "integrated"
```

**修改这一行即可:**

```r
# 分析 cDC1
analysis_mode <- "cDC1"

# 分析 cDC2
analysis_mode <- "cDC2"

# 整合分析 cDC1 + cDC2
analysis_mode <- "integrated"
```

然后运行所有代码块,系统会自动:
- ✅ 加载对应的配置参数
- ✅ 使用对应的缓存文件
- ✅ 生成对应的图表和报告

---

## 📊 三种模式对比

| 特性 | cDC1 | cDC2 | Integrated |
|------|------|------|------------|
| **数据源** | subset_2023.rds.gz | subset2_v2_2023.rds.gz | 两者合并 |
| **细胞类型** | 6 种 cDC1 亚型 | 6 种 cDC2 亚型 | 12 种 (6+6) |
| **RNA 维度** | 25 | 30 | 25 |
| **ADT 维度** | 15 | 20 | 20 |
| **批次校正** | 无 | 无 | ✅ Harmony |
| **分析重点** | cDC1 发育轨迹 | cDC2 成熟过程 | cDC1 vs cDC2 对比 |
| **关键标志物** | Xcr1, Clec9a, Itgae | Sirpa, Cd24a, Notch2 | 两者混合 |

---

## 🎯 典型使用场景

### 场景 1: 首次分析三种模式

```r
# 第 1 天: 分析 cDC1
analysis_mode <- "cDC1"
# 运行所有代码块 → 生成 cDC1_*.rds 缓存 → 得到 cDC1 报告

# 第 2 天: 分析 cDC2
analysis_mode <- "cDC2"
# 运行所有代码块 → 生成 cDC2_*.rds 缓存 → 得到 cDC2 报告

# 第 3 天: 整合分析
analysis_mode <- "integrated"
# 运行所有代码块 → 生成 integrated_*.rds 缓存 → 得到整合报告
```

**总耗时:** 每个模式首次运行约 90-120 分钟

---

### 场景 2: 修改可视化参数

```r
# 已经运行过 cDC1 分析,现在想调整图表
analysis_mode <- "cDC1"  # 保持模式不变

# 修改可视化代码块 (如 DimPlot 参数)
# 重新运行可视化部分即可

# ✅ 优势: 自动加载缓存,只需 2-3 分钟
```

---

### 场景 3: 修改聚类分辨率

假设你想让 cDC2 的聚类更细:

**步骤 1:** 打开 `01_CITE-seq_unified_analysis.Rmd`

**步骤 2:** 找到配置函数中的 cDC2 部分 (约 80 行):

```r
cDC2 = list(
  ...
  resolution = 0.8,  # 改为 1.2
  ...
)
```

**步骤 3:** 删除 cDC2 聚类缓存:

```r
file.remove("results/Cache/cDC2_04_seurat_clustered.rds")
file.remove("results/Cache/cDC2_07_wnn_integrated.rds")
```

**步骤 4:** 重新运行脚本

✅ **优势:** 只重新计算聚类步骤,之前的质控/归一化直接加载缓存

---

## 📁 缓存文件结构

```
results/Cache/
├── cDC1_01_subset_extracted.rds      # cDC1 子集
├── cDC1_04_seurat_clustered.rds      # cDC1 聚类
├── cDC1_07_wnn_integrated.rds        # cDC1 WNN
├── cDC2_01_subset_extracted.rds      # cDC2 子集
├── cDC2_04_seurat_clustered.rds      # cDC2 聚类
├── cDC2_07_wnn_integrated.rds        # cDC2 WNN
├── integrated_01_subset_extracted.rds # 整合子集
├── integrated_04_seurat_clustered.rds # 整合聚类 (含 Harmony)
└── integrated_07_wnn_integrated.rds   # 整合 WNN
```

**特点:** 不同模式的缓存独立,互不干扰

---

## 🛠️ 缓存管理

### 查看特定模式的缓存

```r
# 在 R 控制台运行
source("01_CITE-seq_unified_analysis.Rmd")  # 加载函数

list_mode_cache("cDC1")
#   file                             size_mb mtime
#   cDC1_07_wnn_integrated.rds        412.5  2026-01-06 10:30:00
#   cDC1_04_seurat_clustered.rds      305.2  2026-01-06 10:15:00
#   cDC1_01_subset_extracted.rds      198.3  2026-01-06 10:00:00
```

### 清除特定模式的所有缓存

```r
# ⚠️ 谨慎使用! 会删除该模式的所有缓存
clear_mode_cache("cDC2")
# ✓ Removed 3 cache files for mode: cDC2
```

### 清除特定步骤的缓存

```r
# 仅清除 cDC1 的 WNN 缓存,保留上游计算
file.remove("results/Cache/cDC1_07_wnn_integrated.rds")
```

---

## 🔬 科学分析建议

### 推荐的分析顺序

```
第一步: cDC1 单独分析
  ↓  理解 cDC1 的发育特征和亚群
  ↓
第二步: cDC2 单独分析
  ↓  理解 cDC2 的成熟过程
  ↓
第三步: 整合分析
  ↓  对比 cDC1 vs cDC2
  ↓  发现共同特征和差异
  ↓
第四步: 撰写论文
     使用三种模式的结果互相验证
```

### 整合分析的优势

1. **批次校正验证**
   - 查看 cDC1 和 cDC2 在 UMAP 上是否按生物学聚类
   - 如果按批次聚类 → 批次效应未消除
   - 如果按细胞类型聚类 → 批次校正成功

2. **直接对比**
   ```r
   # 在整合模式下,可以直接比较 cDC1 vs cDC2
   Idents(seurat_obj) <- seurat_obj$batch
   markers_cdc1_vs_cdc2 <- FindMarkers(seurat_obj, ident.1 = "cDC1", ident.2 = "cDC2")
   ```

3. **发育轨迹分析**
   - 使用 Monocle3 或 Slingshot
   - 在整合空间中追踪从前体到成熟的轨迹

---

## ⚡ 性能对比

| 操作 | 传统方案 (3个脚本) | 参数化方案 (本脚本) |
|------|------------------|-------------------|
| 首次运行 3 种模式 | ~300 分钟 | ~300 分钟 |
| 修改参数重跑 | ~300 分钟 | **~30 分钟** |
| 调整可视化 | ~300 分钟 | **~5 分钟** |
| 代码维护 | 修改 3 处 | 修改 1 处 |
| 质控标准统一性 | 容易不一致 | ✅ 完全一致 |

---

## 🎨 自定义配置

### 修改细胞类型

假设你只想分析部分亚型:

```r
# 在配置函数中修改
cDC1 = list(
  cell_types = c("Early mature cDC1s", "Late mature cDC1s"),  # 只分析成熟期
  ...
)
```

### 修改聚类参数

```r
cDC2 = list(
  dims_rna = 35,        # 增加 RNA 维度
  dims_adt = 25,        # 增加 ADT 维度
  resolution = 1.2,     # 更细的聚类
  ...
)
```

### 更换批次校正方法

```r
integrated = list(
  integration_method = "CCA",  # 改用 Seurat CCA
  # 或 "RPCA"
  ...
)
```

---

## ❓ 常见问题

### Q1: 切换模式后缓存冲突?

**A:** 不会! 每个模式有独立的缓存前缀:
- cDC1: `cDC1_*.rds`
- cDC2: `cDC2_*.rds`
- integrated: `integrated_*.rds`

### Q2: 如何确认正在使用哪个模式?

**A:** 运行脚本时会显示:

```
================================================================================
🎯 当前分析模式: cDC1
================================================================================
```

### Q3: 整合分析需要先运行单独分析吗?

**A:** 不需要! 三种模式完全独立,可以任意顺序运行。

### Q4: 缓存文件太大怎么办?

**A:** 定期清理不需要的模式:

```r
# 如果不需要 cDC2,清除其缓存
clear_mode_cache("cDC2")

# 或者只保留最终结果,删除中间缓存
file.remove(list.files("results/Cache/", pattern = "_0[1-6]", full.names = TRUE))
```

### Q5: 如何批量生成三种模式的报告?

```r
# 方法 1: 手动运行三次
# (每次修改 analysis_mode 并渲染)

# 方法 2: 使用脚本批量渲染
modes <- c("cDC1", "cDC2", "integrated")
for (mode in modes) {
  rmarkdown::render(
    "01_CITE-seq_unified_analysis.Rmd",
    output_file = paste0("results/Reports/", mode, "_analysis.html"),
    params = list(mode = mode)
  )
}
```

---

## 📞 获取帮助

- **设计文档:** `ANALYSIS_FRAMEWORK_DESIGN.md` (详细架构说明)
- **优化文档:** `OPTIMIZATION_README.md` (缓存系统说明)
- **本文档:** `QUICK_START_GUIDE.md`

---

**祝分析愉快! 🎉**

**版本:** v1.0
**日期:** 2026-01-06
**作者:** Claude Code
