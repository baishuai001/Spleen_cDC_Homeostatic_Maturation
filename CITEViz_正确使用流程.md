# CITEViz 正确使用流程

## 您的疑问是对的！

您提出的问题非常关键：**在启动脚本中下采样没有意义，因为 CITEViz 启动后还需要上传文件。**

这是正确的观察！我已经修正了流程。

---

## 正确的两步流程

### 第 1 步：准备数据文件

**运行一次**（准备不同大小的文件供上传）：

```r
# 在 R 或 RStudio 中运行
source("prepare_citeviz_data.R")
```

**输出文件**（3个）：
```
citeviz_fast_5k.rds       ← 5000 细胞（推荐首次使用）
citeviz_medium_8k.rds     ← 8000 细胞（平衡版）
citeviz_full_compatible.rds ← 所有细胞（完整版）
```

**为什么需要这一步**：
- ✅ 创建兼容 CITEViz 的对象格式
- ✅ 预先下采样，生成不同大小的文件
- ✅ 您可以选择合适大小的文件上传
- ✅ 下采样在文件中，上传后就是采样后的数据

---

### 第 2 步：启动 CITEViz 并上传文件

**每次使用时运行**：

```r
# 在 R 或 RStudio 中运行
source("launch_citeviz.R")
```

**操作流程**：
1. 启动后，CITEViz 在浏览器中打开
2. 看到 "Browse..." 按钮
3. 点击上传文件：
   - 首次尝试：`citeviz_fast_5k.rds`（最快）
   - 如果不卡：`citeviz_medium_8k.rds`
   - 最终分析：`citeviz_full_compatible.rds`

---

## 对比：错误 vs 正确的理解

### ❌ 错误理解（我之前的脚本）

```r
# start_citeviz.R 的问题：
seurat_obj <- 加载数据
seurat_obj <- 下采样(seurat_obj)  # ← 在内存中下采样
run_app(seurat_object = seurat_obj) # ← 启动应用

# 问题：如果 run_app() 不接受对象参数，
# 启动后仍需上传文件，内存中的下采样就浪费了！
```

### ✅ 正确理解（新的流程）

```r
# prepare_citeviz_data.R：
seurat_obj <- 加载数据
seurat_sampled <- 下采样(seurat_obj)
保存文件(seurat_sampled, "citeviz_fast_5k.rds")  # ← 保存到文件

# launch_citeviz.R：
run_app()  # ← 只启动应用

# 然后在网页上上传 citeviz_fast_5k.rds
```

**关键差异**：
- ❌ 旧方式：下采样在内存中，可能浪费
- ✅ 新方式：下采样在文件中，上传后生效

---

## 为什么要准备 3 个不同大小的文件？

| 文件 | 细胞数 | 用途 | 响应速度 |
|------|--------|------|----------|
| **fast_5k** | 5,000 | 快速探索、门控测试 | ⚡️⚡️⚡️⚡️⚡️ |
| **medium_8k** | 8,000 | 常规分析、平衡版 | ⚡️⚡️⚡️⚡️ |
| **full** | 18,193 | 最终分析、出图 | ⚡️⚡️ |

**使用策略**：
1. 先用 `fast_5k` 快速浏览和测试门控策略
2. 确定分析思路后，用 `medium_8k` 或 `full` 做最终分析

---

## 完整使用示例

### 首次使用（完整流程）

```r
# 1. 设置工作目录
setwd("/home/user/Spleen_cDC_Homeostatic_Maturation")

# 2. 准备数据文件（只需运行一次）
source("prepare_citeviz_data.R")

# 输出：
# ✓ 保存: citeviz_fast_5k.rds
# ✓ 保存: citeviz_medium_8k.rds
# ✓ 保存: citeviz_full_compatible.rds

# 3. 启动 CITEViz
source("launch_citeviz.R")

# 4. 在浏览器中：
#    - 点击 "Browse..."
#    - 选择 citeviz_fast_5k.rds
#    - 开始分析
```

### 后续使用（快速启动）

```r
# 直接启动（文件已准备好）
source("launch_citeviz.R")

# 然后在网页上传文件
```

---

## 脚本对比表

| 脚本 | 功能 | 运行频率 | 输出 |
|------|------|----------|------|
| **prepare_citeviz_data.R** | 准备数据文件 | 运行一次 | 3 个 RDS 文件 |
| **launch_citeviz.R** | 启动 CITEViz | 每次使用 | 浏览器应用 |
| ~~start_citeviz.R~~ | （旧脚本，已弃用） | - | - |

---

## 常见问题

### Q1: 为什么不能在启动时直接传递对象？

**A**: CITEViz 的设计可能是：
- 通过文件上传界面加载数据（更灵活）
- 而不是通过函数参数传递对象

即使 `run_app(seurat_object = obj)` 语法存在，也可能：
- 不稳定
- 或仍然会显示上传界面

### Q2: 我可以只准备一个文件吗？

**A**: 可以！修改 `prepare_citeviz_data.R`：

```r
# 只保留快速版
seurat_fast <- smart_sample(seurat_obj, 5000)
seurat_fast <- create_compatible(seurat_fast)
saveRDS(seurat_fast, "citeviz_data.rds")
```

### Q3: 文件在哪里？

**A**: 在项目根目录：
```bash
ls -lh citeviz*.rds
```

---

## 总结

**您的观察完全正确**！下采样应该：
1. ✅ **在文件准备阶段** - 保存为不同大小的文件
2. ❌ **不是在启动阶段** - 启动只负责打开应用

**正确流程**：
```
准备数据（一次） → 启动应用（每次） → 上传文件（选择大小） → 开始分析
```

感谢您的仔细观察！这确保了流程的正确性。
