# CITEViz 工具快速开始指南

## 前置要求

### 必需文件
在使用 CITEViz 工具之前，您需要一个 Seurat 对象文件，通常名为：
- `seurat_obj_annotated.rds`

### 如何生成
运行下游分析脚本：

```r
# 在 R 中运行
rmarkdown::render('4.downstream_analysis.Rmd')
```

或

```r
# 使用 knitr
knitr::knit('4.downstream_analysis.Rmd')
```

这将生成文件：`results/integrated/Robjects/seurat_obj_annotated.rds`

---

## 使用指南

### 选项 1: 一键运行（推荐）

最简单的方式，自动完成所有步骤：

```bash
Rscript run_citeviz_complete.R
```

**优点：**
- 无需手动干预
- 自动转换格式
- 自动启动 CITEViz

**缺点：**
- 如果某步失败，需要手动排查

---

### 选项 2: 分步执行（推荐用于排查问题）

#### 步骤 1: 检查文件是否存在

```bash
Rscript find_seurat_object.R
```

这将：
- 搜索所有可能的 Seurat 对象文件
- 显示文件大小和位置
- 给出下一步建议

#### 步骤 2: 转换格式（Seurat v5 → v4）

```bash
Rscript convert_seurat_v5_to_v4_for_citeviz.R
```

这将：
- 自动检测输入文件
- 转换 Assay5 → Assay 格式
- 降采样到 10,000 细胞（可配置）
- 生成 v4 兼容文件

#### 步骤 3: 启动 CITEViz

```bash
Rscript launch_citeviz.R --preload
```

或（空启动，需手动上传文件）：

```bash
Rscript launch_citeviz.R
```

---

### 选项 3: 使用替代工具（最稳定，**强烈推荐**）

如果 CITEViz 持续出现问题，使用 Seurat 原生分析：

```bash
Rscript adt_analysis_seurat.R
```

**优点：**
- ✓ 100% Seurat v5 兼容
- ✓ 自动检测输入文件
- ✓ 生成出版级 PDF 图表
- ✓ 无需格式转换

**输出：**
- `results/integrated/adt_analysis/*.pdf`
- 包括 UMAP、热图、点图、山脊图等

---

## 常见问题

### Q1: 运行时提示"未找到 Seurat 对象文件"

**原因：** 您还没有运行前面的分析步骤

**解决：**
1. 运行 `4.downstream_analysis.Rmd` 生成 Seurat 对象
2. 或运行 `find_seurat_object.R` 查找现有文件

### Q2: CITEViz 显示 `[object Object]` 错误

**原因：** Seurat v5 兼容性问题

**解决：**
1. 确保已运行格式转换：`Rscript convert_seurat_v5_to_v4_for_citeviz.R`
2. 使用预加载模式：`Rscript launch_citeviz.R --preload`
3. 或改用替代工具：`Rscript adt_analysis_seurat.R`

### Q3: 转换脚本报错"cannot open connection"

**原因：** 输入文件不存在或路径错误

**解决：**
1. 检查文件是否存在：`ls -lh results/integrated/Robjects/seurat_obj_annotated.rds`
2. 运行检测脚本：`Rscript find_seurat_object.R`
3. 手动指定路径：修改转换脚本中的 `INPUT_FILE` 变量

### Q4: Rscript 命令不可用

**原因：** R 未安装或不在 PATH 中

**解决：**
1. 安装 R：`conda install r-base` 或从官网下载
2. 或在 R 中直接运行：`source('run_citeviz_complete.R')`

---

## 脚本清单

| 脚本 | 用途 | 何时使用 |
|------|------|---------|
| `find_seurat_object.R` | 查找可用文件 | 首次使用 / 找不到文件时 |
| `convert_seurat_v5_to_v4_for_citeviz.R` | 格式转换 | 使用 CITEViz 前 |
| `launch_citeviz.R` | 启动 CITEViz | 格式转换后 |
| `run_citeviz_complete.R` | 一键运行 | 快速使用 |
| `adt_analysis_seurat.R` | 替代方案 | CITEViz 失败时 |
| `check_citeviz_installation.R` | 检查安装 | 排查安装问题 |

---

## 配置选项

### 转换脚本配置

编辑 `convert_seurat_v5_to_v4_for_citeviz.R`：

```r
# 降采样设置
DOWNSAMPLE <- TRUE    # 是否降采样
MAX_CELLS <- 10000    # 最大细胞数（推荐 5000-10000）

# 如果文件太大，可以降低细胞数
MAX_CELLS <- 5000     # 更快加载
```

### ADT 分析配置

编辑 `adt_analysis_seurat.R`：

```r
# 如果需要分析特定的 ADT 标记
key_adts <- c(
  "adt-CD11c",
  "adt-CD317",
  "adt-Sirpa",
  # 添加您感兴趣的标记
)
```

---

## 完整工作流示例

### 场景 1: 首次使用

```bash
# 1. 检查是否有数据文件
Rscript find_seurat_object.R

# 如果没有，先运行下游分析（在 R 中）
R -e "rmarkdown::render('4.downstream_analysis.Rmd')"

# 2. 一键运行 CITEViz
Rscript run_citeviz_complete.R
```

### 场景 2: CITEViz 出现问题

```bash
# 1. 检查安装
Rscript check_citeviz_installation.R

# 2. 尝试分步执行
Rscript convert_seurat_v5_to_v4_for_citeviz.R
Rscript launch_citeviz.R --preload

# 3. 如果仍失败，使用替代方案
Rscript adt_analysis_seurat.R
```

### 场景 3: 只想快速分析 ADT 数据

```bash
# 直接使用 Seurat 分析（跳过 CITEViz）
Rscript adt_analysis_seurat.R
```

---

## 技术支持

### 查看详细文档

- 完整解决方案：`CITEViz_完整解决方案.md`
- 包含故障排查、技术细节、方案对比

### 检查日志

所有脚本都会输出详细的运行日志，包括：
- 文件检测结果
- 转换进度
- 错误信息

### 常用调试命令

```bash
# 检查 R 版本
R --version

# 检查 Seurat 版本
R -e "packageVersion('Seurat')"

# 检查文件大小
ls -lh results/integrated/Robjects/*.rds

# 在 R 中检查对象
R -e "obj <- readRDS('path/to/file.rds'); print(class(obj)); print(names(obj@assays))"
```

---

## 推荐流程

基于我们的经验，推荐以下流程：

1. **如果您是新用户** → 使用选项 3（`adt_analysis_seurat.R`）
   - 最简单、最稳定
   - 生成的图表质量高
   - 避免 CITEViz 兼容性问题

2. **如果您需要交互式圈门** → 使用选项 2（分步执行）
   - 可以控制每个步骤
   - 便于排查问题
   - 灵活性高

3. **如果您熟悉工具** → 使用选项 1（一键运行）
   - 快速启动
   - 适合重复使用

---

## 更新日志

- **2026-01-28**: 添加自动文件检测功能
- **2026-01-28**: 修复 Seurat v5 兼容性问题
- **2026-01-28**: 创建 find_seurat_object.R
- **2026-01-23**: 初始版本

---

## 相关文件

- `CITEViz_完整解决方案.md` - 详细技术文档
- `CITEViz_使用教程.Rmd` - R Markdown 教程
- `CITEViz_使用说明.md` - 使用说明
