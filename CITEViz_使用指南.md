# CITEViz ADT数据分析指南

## 概述

CITEViz是一个交互式的R-Shiny应用程序,可以像FlowJo一样对CITE-seq的ADT(抗体衍生标签)数据进行gating分析。它提供了流式细胞术风格的工作流程,让您能够通过交互式gating来分类和可视化细胞群。

## 前置要求

### 输入数据
- 已处理的Seurat对象(来自`3.RNA-ADT_HPCscript.R`和`4.script_CITEseq_SAM_WT_aggr.R`)
- 包含归一化的ADT assay(CLR转换)
- 包含降维结果(UMAP/tSNE)
- 包含细胞注释信息

### 软件要求
- R (>= 4.0)
- 已安装的R包:
  - Seurat
  - devtools
  - CITEViz
  - dplyr

## 安装CITEViz

```r
# 首次使用需要安装CITEViz
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# 从GitHub安装CITEViz
devtools::install_github("maxsonBraunLab/CITE-Viz")
```

**注意**: GitHub仓库名为"CITE-Viz"(带连字符),但R包名为"CITEViz"(无连字符)。

## 使用步骤

### 1. 准备Seurat对象

确保您的Seurat对象包含:
- ✅ ADT assay with normalized data (CLR transformation)
- ✅ 降维结果 (ADT_umap, ADT_tsne, SCT_umap, SCT_tsne)
- ✅ 细胞聚类和注释信息

### 2. 运行分析脚本

```bash
# 在R或RStudio中运行
Rscript 10.script_CITEViz_ADT_analysis.R
```

或者在R/RStudio交互界面中:
```r
source("10.script_CITEViz_ADT_analysis.R")
```

### 3. 修改输入路径

在运行脚本前,需要修改第35行的输入文件路径:

```r
# 选择以下选项之一:

# 选项1: 完整的WT_aggr对象
input_seurat <- "SAM2and3_WT/results/Robjects/seuratObj_SAM2and3_WT_ADT.rds"

# 选项2: cDC1子集
input_seurat <- "SAM2and3_WT_subset1/results/Robjects/seuratObj_SAM2and3_WT_subset1_ADT.rds"

# 选项3: cDC2子集
input_seurat <- "SAM2and3_WT_subset2/results/Robjects/seuratObj_SAM2and3_WT_subset2_ADT.rds"
```

## CITEViz交互界面使用

### 启动应用
运行脚本后,CITEViz会在您的默认浏览器中打开一个交互式界面。

### 主要功能模块

#### 1. **Gating标签页** (核心功能)
类似FlowJo的gating工作流:

- **选择标记物**:
  - X轴和Y轴选择不同的ADT标记物
  - 例如: X轴=CD11c, Y轴=CD8a

- **绘制门控(Gate)**:
  - 使用鼠标在散点图上圈选细胞群
  - 支持多边形、矩形等多种gate形状
  - 可以进行序列gating(在已gated的细胞中继续gate)

- **实时可视化**:
  - Gated的细胞会立即在UMAP/tSNE图上高亮显示
  - 可以看到选中细胞的空间分布

- **保存Gate**:
  - 为每个gate命名(例如: "Resident_cDC1s", "Mig_cDC1s")
  - Gate信息会保存到Seurat对象的metadata中

#### 2. **QC标签页**
查看质量控制指标:
- RNA counts分布
- Gene counts分布
- ADT counts分布
- 可以按不同metadata分组查看

#### 3. **Feature Plot标签页**
可视化标记物表达:
- 在降维图上显示特定ADT标记物的表达
- 支持多个标记物同时显示
- 可以调整颜色范围和cutoff值

### 推荐的Gating策略

基于您的研究,以下是一些有用的ADT标记物组合:

#### 识别主要DC群体:
```
1. 所有DCs: CD11c+ MHC-II+
2. cDC1s: CD11c+ CD8a+ CD11b-
3. cDC2s: CD11c+ CD11b+ CD172a+
4. 成熟DCs: CD86+ MHC-II(high)
5. 迁移性DCs: CCR7+ CD86+
6. 驻留cDC2s: ESAM+ CD11b+
```

#### Gating顺序示例:
```
Gate 1: CD11c+ cells (所有DCs)
  └─ Gate 2a: CD11c+ CD8a+ (cDC1s)
      └─ Gate 3a: CD11c+ CD8a+ CCR7+ (迁移性cDC1s)
      └─ Gate 3b: CD11c+ CD8a+ CCR7- (驻留cDC1s)
  └─ Gate 2b: CD11c+ CD11b+ (cDC2s)
      └─ Gate 3c: CD11c+ CD11b+ ESAM+ (ESAM+ 驻留cDC2s)
      └─ Gate 3d: CD11c+ CD11b+ ESAM- (ESAM- 驻留cDC2s)
```

## 保存和导出结果

### 1. 保存带有Gate信息的Seurat对象

```r
# CITEViz会在metadata中添加gate信息
# 关闭CITEViz后保存更新的对象
saveRDS(seuratObj, file = "results/seuratObj_with_CITEViz_gates.rds")
```

### 2. 提取特定gated群体

```r
# 提取特定gate的细胞
gated_cells <- WhichCells(seuratObj, expression = CITEViz_gate == "Resident_cDC1s")

# 创建子集对象进行下游分析
subset_obj <- subset(seuratObj, cells = gated_cells)

# 可视化
DimPlot(seuratObj, reduction = "ADT_umap", group.by = "CITEViz_gate")
```

### 3. 导出Gate统计信息

```r
# 统计每个gate的细胞数
table(seuratObj$CITEViz_gate)

# 导出为CSV
gate_stats <- as.data.frame(table(seuratObj$CITEViz_gate))
write.csv(gate_stats, "results/CITEViz_gate_statistics.csv")
```

## 手动Gating(可选)

如果您不想使用交互界面,也可以使用Seurat函数进行手动gating:

```r
# 设置ADT为默认assay
DefaultAssay(seuratObj) <- "ADT"

# 示例: Gate CD11c+ CD8a+ cells (cDC1s)
cd11c_high <- WhichCells(seuratObj, expression = CD11c > 2)
cd8a_high <- WhichCells(seuratObj, expression = CD8a > 1.5)
cDC1_gated <- intersect(cd11c_high, cd8a_high)

# 添加到metadata
seuratObj$manual_gate_cDC1 <- ifelse(
  colnames(seuratObj) %in% cDC1_gated,
  "cDC1",
  "Other"
)

# 可视化
DimPlot(seuratObj, reduction = "ADT_umap", group.by = "manual_gate_cDC1")
```

## CITEViz的优势

✅ **交互式gating**: 类似FlowJo的直观界面
✅ **实时可视化**: 立即在UMAP/tSNE上看到gated细胞
✅ **序列gating**: 支持在已gated的细胞中继续细分
✅ **无需手动设定阈值**: 通过可视化直接圈选
✅ **可重复性**: 导出gate参数用于后续分析
✅ **整合分析**: 同时查看ADT和RNA数据

## 故障排除

### 问题1: Matrix包版本错误

**错误信息**:
```
Error: package or namespace load failed for 'SeuratObject' in loadNamespace...
载入了名字空间'Matrix' 1.6-2，但需要的是>= 1.6.3
```

**解决方案**:

#### 方法1: 使用快速修复脚本 (推荐)
```r
source("quick_fix_matrix.R")
```

#### 方法2: 使用完整修复脚本
```r
source("fix_dependencies.R")
```

#### 方法3: 手动修复
```r
# 1. 卸载旧版本
remove.packages("Matrix")

# 2. 重启R会话
# RStudio: Session > Restart R
# 或命令行: q() 然后重新启动

# 3. 安装最新版本
install.packages("Matrix")

# 4. 验证版本
packageVersion("Matrix")  # 应该 >= 1.6.3

# 5. 重新加载Seurat
library(Seurat)
```

### 问题2: CITEViz安装失败

**可能原因**: 缺少依赖包或网络问题

**解决方案**:
```r
# 确保安装所有依赖
install.packages(c("shiny", "plotly", "DT", "ggplot2", "dplyr"))

# 使用devtools安装CITEViz
install.packages("devtools")
devtools::install_github("maxsonBraunLab/CITE-Viz", dependencies = TRUE)

# 如果GitHub访问有问题，可以尝试镜像或下载zip包手动安装
```

### 问题3: Seurat对象缺少ADT assay

**错误信息**: `Error: Seurat object does not contain ADT assay`

**解决方案**:
确保使用的是从`3.RNA-ADT_HPCscript.R`生成的包含ADT数据的Seurat对象。检查:
```r
# 检查可用的assays
names(seuratObj@assays)
# 应该包含: "RNA", "ADT", "SCT"

# 如果缺少ADT assay，可能需要重新运行script 3
```

### 问题4: ADT数据未归一化

**警告信息**: `Warning: ADT assay may not be normalized`

**解决方案**:
脚本会自动进行CLR归一化，但你也可以手动执行:
```r
seuratObj <- NormalizeData(
  seuratObj,
  assay = "ADT",
  normalization.method = "CLR"
)
```

### 问题5: CITEViz应用无法启动

**可能原因**: 端口占用或Shiny配置问题

**解决方案**:
```r
# 尝试指定不同的端口
options(shiny.port = 8888)
run_app(seurat_object = seuratObj)

# 或者在浏览器中手动打开
options(shiny.launch.browser = TRUE)
run_app(seurat_object = seuratObj)
```

### 问题6: 内存不足

**症状**: R崩溃或响应缓慢

**解决方案**:
```r
# 1. 增加R的内存限制 (Windows)
memory.limit(size = 16000)  # 16GB

# 2. 使用子集进行分析
# 先分析cDC1或cDC2子集，而不是完整对象
subset_obj <- subset(seuratObj, subset = SCT_clusters %in% c("0", "1", "2"))

# 3. 降采样
subset_obj <- subset(seuratObj, downsample = 5000)
```

### 问题7: R包版本冲突

**解决方案**:
```r
# 更新所有包到最新版本
update.packages(ask = FALSE)

# 检查关键包版本
packageVersion("Seurat")        # 推荐 >= 4.0
packageVersion("SeuratObject")  # 推荐 >= 4.0
packageVersion("Matrix")        # 必须 >= 1.6.3
packageVersion("dplyr")         # 推荐 >= 1.0

# 如果问题持续，考虑重新安装R (>= 4.2)
```

### 问题8: 找不到特定ADT标记物

**检查可用标记物**:
```r
# 查看所有可用的ADT标记物
rownames(seuratObj@assays$ADT)

# 在CITEViz中，ADT标记物名称应该与此列表匹配
```

## 常见问题

### Q1: CITEViz需要什么格式的输入?
A: 需要预处理的Seurat对象,包含归一化的ADT数据(推荐CLR转换)。

### Q2: 可以同时分析RNA和ADT数据吗?
A: 是的,CITEViz可以在ADT数据上gating,同时在RNA的UMAP上可视化。

### Q3: Gate信息如何保存?
A: Gate信息会作为新的metadata列添加到Seurat对象中。

### Q4: 可以导出gate坐标吗?
A: 可以,CITEViz允许导出gate的坐标和参数,用于重现分析。

### Q5: 支持批次效应校正吗?
A: CITEViz本身不做预处理,但可以使用已经校正过批次效应的Seurat对象。

## 参考资源

- **GitHub仓库**: https://github.com/maxsonBraunLab/CITEViz
- **发表文章**: CITEViz: interactively classify cell populations in CITE-Seq via a flow cytometry-like gating workflow using R-Shiny. BMC Bioinformatics (2024)
- **在线文档**: https://maxsonbraunlab.github.io/CITEViz/
- **Workshop教程**: https://github.com/gartician/CITEVizWorkshop

## 技术支持

如有问题,可以:
1. 查看GitHub Issues: https://github.com/maxsonBraunLab/CITEViz/issues
2. 参考BMC Bioinformatics论文中的方法学描述
3. 查看Bioconductor workshop材料

## 与FlowJo的对比

| 功能 | FlowJo | CITEViz |
|------|--------|---------|
| 交互式gating | ✅ | ✅ |
| 2D散点图 | ✅ | ✅ |
| 序列gating | ✅ | ✅ |
| 多种gate形状 | ✅ | ✅ |
| 整合降维可视化 | ❌ | ✅ |
| 整合RNA数据 | ❌ | ✅ |
| 单细胞分辨率 | ✅ | ✅ |
| 处理大规模数据 | 有限 | 更好 |

## 下游分析建议

使用CITEViz gating后,可以进行:

1. **差异表达分析**: 比较不同gated群体的基因表达
2. **轨迹推断**: 在gated的细胞上进行拟时序分析
3. **功能富集**: 对gated群体的marker基因进行GO/KEGG分析
4. **细胞通讯**: 分析不同gated群体之间的细胞间通讯
5. **亚群分析**: 对特定gated群体进行更精细的聚类

## 脚本更新日志

- **v1.0** (2024): 初始版本,支持基本的CITEViz工作流

---

**作者**: Claude
**日期**: 2024
**用途**: 分析Spleen cDC Homeostatic Maturation项目中的CITE-seq ADT数据
