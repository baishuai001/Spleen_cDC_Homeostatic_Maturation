# CITEViz 完整解决方案

## 问题诊断

### 根本原因
CITEViz 是基于 **Seurat v3/v4** 开发的，无法识别 **Seurat v5** 的新数据结构（`Assay5` 对象）。

### 错误症状
- `Error: [object Object]` - CITEViz 无法读取 Assay5 对象
- `Error: 缺少参数"..4"` - ggplot2 与旧代码冲突
- 空白图表或崩溃

---

## 解决方案概览

我们提供了 **3 种解决方案**，按推荐程度排序：

| 方案 | 难度 | 成功率 | 说明 |
|------|------|--------|------|
| **方案 1: 一键运行** | ⭐ 最简单 | 中等 | 自动转换 + 启动 CITEViz |
| **方案 2: 分步执行** | ⭐⭐ 简单 | 中等 | 手动控制每个步骤 |
| **方案 3: 使用替代工具** | ⭐ 最简单 | ⭐⭐⭐ 最高 | 放弃 CITEViz，使用 Seurat 原生分析 |

---

## 方案 1: 一键运行（推荐首选）

### 使用方法
```bash
Rscript run_citeviz_complete.R
```

### 工作流程
1. 自动将 Seurat v5 对象转换为 v4 格式
2. 降采样到 10K 细胞（加快加载速度）
3. 启动 CITEViz 并预加载数据

### 优点
- 无需手动干预
- 一条命令完成所有步骤

### 缺点
- 如果某个步骤失败，需要手动排查

---

## 方案 2: 分步执行（推荐用于排查问题）

### 步骤 1: 检查 CITEViz 安装
```bash
Rscript check_citeviz_installation.R
```

**预期输出：**
- ✓ CITEViz 已安装
- ✓ 包可以加载
- ⚠ Seurat v5 可能不兼容

### 步骤 2: 转换数据格式
```bash
Rscript convert_seurat_v5_to_v4_for_citeviz.R
```

**关键参数（可在脚本中修改）：**
```r
INPUT_FILE <- "results/integrated/Robjects/seurat_obj_annotated.rds"
OUTPUT_FILE <- "results/integrated/Robjects/seurat_obj_citeviz_v4_compatible.rds"
DOWNSAMPLE <- TRUE
MAX_CELLS <- 10000  # 降采样细胞数
```

**预期输出：**
- ✓ 所有 Assay 已转换为 v3/v4 格式
- ✓ 保存成功

**转换原理：**
```r
# Seurat v5 使用 Assay5 类
class(seurat_obj@assays[["ADT"]])  # [1] "Assay5"

# 转换为 v3/v4 的 Assay 类
seurat_obj@assays[["ADT"]] <- as(seurat_obj@assays[["ADT"]], "Assay")

# 转换后
class(seurat_obj@assays[["ADT"]])  # [1] "Assay"
```

### 步骤 3: 启动 CITEViz

#### 方式 A: 预加载数据（推荐）
```bash
Rscript launch_citeviz.R --preload
```

**优点：**
- 直接加载转换后的 v4 兼容文件
- 无需在网页中手动上传
- 减少兼容性问题

#### 方式 B: 空启动（需手动上传）
```bash
Rscript launch_citeviz.R
```

然后在 CITEViz 网页中：
1. 点击 "Browse..."
2. 选择 `results/integrated/Robjects/seurat_obj_citeviz_v4_compatible.rds`
3. 等待加载

---

## 方案 3: 使用替代工具（最稳定）

如果 CITEViz 持续出错，**强烈建议使用这个方案**。

### 使用 Seurat 原生 ADT 分析
```bash
Rscript adt_analysis_seurat.R
```

### 功能对比

| 功能 | CITEViz | adt_analysis_seurat.R |
|------|---------|----------------------|
| 交互式圈门 | ✓ | ✗ |
| UMAP 可视化 | ✓ | ✓ |
| ADT 表达热图 | ✓ | ✓ |
| 相关性分析 | 部分 | ✓ |
| 统计图表 | 基础 | ✓ 详细 |
| 稳定性 | ⚠ v5 不兼容 | ✓ 完全兼容 |
| 输出 | 网页交互 | PDF 报告 |

### 输出内容
- `results/integrated/adt_analysis/`
  - `adt_umap_overview.pdf` - UMAP 总览
  - `adt_expression_heatmap.pdf` - 表达热图
  - `adt_dotplot.pdf` - 点图
  - `adt_ridgeplot.pdf` - 山脊图
  - `key_markers_comparison.pdf` - 关键标记比较
  - `adt_correlation_matrix.pdf` - 相关性矩阵

### 优点
- ✓ 100% Seurat v5 兼容
- ✓ 无需转换数据
- ✓ 生成出版级图表
- ✓ 可编程扩展

### 缺点
- ✗ 无交互式圈门功能（需要的话可用 FlowJo 或 Cytobank）

---

## 故障排查

### 问题 1: 转换后仍报错 `[object Object]`

**可能原因：**
- 转换不完整
- CITEViz 版本过旧
- 依赖包冲突

**解决方法：**
```r
# 检查转换结果
seurat_obj <- readRDS("results/integrated/Robjects/seurat_obj_citeviz_v4_compatible.rds")

# 检查 Assay 类型
sapply(names(seurat_obj@assays), function(x) {
  class(seurat_obj@assays[[x]])[1]
})

# 预期输出（全部应该是 "Assay"，不应该是 "Assay5"）
#    RNA    ADT
# "Assay" "Assay"
```

如果仍有 `Assay5`，说明转换失败。

### 问题 2: 转换脚本报错

**常见错误：**
```
Error: Can't find input file
```

**解决方法：**
检查文件路径是否正确：
```bash
ls -lh results/integrated/Robjects/seurat_obj_annotated.rds
```

### 问题 3: CITEViz 启动但无法加载数据

**可能原因：**
- 文件太大（>2GB）
- 内存不足
- 浏览器缓存问题

**解决方法：**
1. 增加降采样力度（改为 5000 细胞）
2. 清除浏览器缓存
3. 使用 Chrome 或 Firefox 浏览器

### 问题 4: `Error: 缺少参数"..4"`

**原因：**
- ggplot2 版本与 CITEViz 不兼容

**临时解决（不推荐）：**
```r
# 降级 ggplot2（可能影响其他脚本）
remotes::install_version("ggplot2", version = "3.3.6")
```

**推荐解决：**
使用方案 3（`adt_analysis_seurat.R`）

---

## 技术细节

### Seurat v5 vs v4 数据结构对比

#### Seurat v4/v3（CITEViz 兼容）
```r
# 数据存储在 slots 中
seurat@assays$ADT@data          # 归一化数据
seurat@assays$ADT@counts        # 原始计数
seurat@assays$ADT@scale.data    # 标准化数据

# 使用旧 API
GetAssayData(seurat, assay = "ADT", slot = "data")
```

#### Seurat v5（CITEViz 不兼容）
```r
# 数据存储在 layers 中
seurat@assays$ADT@layers$data          # 归一化数据
seurat@assays$ADT@layers$counts        # 原始计数
seurat@assays$ADT@layers$scale.data    # 标准化数据

# 使用新 API
GetAssayData(seurat, assay = "ADT", layer = "data")
```

### 转换关键代码
```r
# 方法 1: 使用 as() 强制转换
seurat@assays[["ADT"]] <- as(seurat@assays[["ADT"]], "Assay")

# 方法 2: 手动重建（如果方法 1 失败）
counts <- LayerData(assay_v5, layer = "counts")
data <- LayerData(assay_v5, layer = "data")
new_assay <- CreateAssayObject(counts = counts)
new_assay <- SetAssayData(new_assay, slot = "data", new.data = data)
```

---

## 最终建议

### 如果你的目标是快速分析 ADT 数据
→ **直接使用方案 3**（`adt_analysis_seurat.R`）

### 如果你必须使用 CITEViz 的交互式圈门功能
→ **使用方案 1**（`run_citeviz_complete.R`），如果失败则改用方案 2 逐步排查

### 如果 CITEViz 持续失败
→ **考虑降级 Seurat 到 v4**（不推荐，可能影响其他分析）

```r
# 降级到 Seurat v4（慎用！）
remotes::install_version("Seurat", version = "4.4.0")
remotes::install_version("SeuratObject", version = "4.1.4")
```

---

## 文件清单

| 文件 | 用途 | 何时使用 |
|------|------|---------|
| `check_citeviz_installation.R` | 检查安装 | 首次使用 |
| `convert_seurat_v5_to_v4_for_citeviz.R` | 格式转换 | 方案 1/2 |
| `launch_citeviz.R` | 启动 CITEViz | 方案 2 |
| `run_citeviz_complete.R` | 一键运行 | 方案 1 |
| `adt_analysis_seurat.R` | 替代方案 | 方案 3 |

---

## 联系方式

如果所有方案都失败，可能的原因：
1. CITEViz 已不再维护，与新版 R 包不兼容
2. 特定系统环境问题

建议：
- 查看 CITEViz GitHub Issues
- 联系软件作者
- 使用其他 CITE-seq 分析工具（Seurat, Scanpy, AnnData）
