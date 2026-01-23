# CITEViz 使用说明

## 两种使用方式对比

### 方式 1：直接加载服务器数据（推荐）✅

```bash
cd /home/user/Spleen_cDC_Homeostatic_Maturation
Rscript run_citeviz.R
```

**优点**：
- ✅ 无需上传，直接读取服务器文件
- ✅ 启动快速
- ✅ 节省网络带宽
- ✅ 内存使用高效

**资源占用**：
- 📍 位置：云服务器
- 💾 内存：服务器内存
- 🔢 计算：服务器 CPU
- 📊 数据：服务器磁盘

---

### 方式 2：网页上传 RDS 文件

```bash
# 先启动 CITEViz（无参数）
R
> library(CITEViz)
> run_app()

# 然后在网页界面点击 "上传文件"
```

**缺点**：
- ❌ 需要先下载到本地，再上传到服务器（多此一举）
- ❌ 大文件上传慢（占用网络带宽）
- ❌ 临时文件占用磁盘空间

**资源占用**：
- 📍 位置：仍然是云服务器（不是本地！）
- 💾 内存：服务器内存
- 🔢 计算：服务器 CPU
- 📊 数据：服务器磁盘 + 浏览器临时缓存

**重要**：即使是"上传"方式，**所有计算仍在服务器进行**！

---

## 使用流程

### 1. 直接运行脚本（推荐）

```bash
# 在项目目录下
cd /home/user/Spleen_cDC_Homeostatic_Maturation

# 运行脚本
Rscript run_citeviz.R
```

输出：
```
=================================================
CITEViz 启动脚本
=================================================

正在加载 Seurat 对象...
路径: results/cDC1/Robjects/seurat_obj_annotated.rds
✓ 加载成功
  细胞数: 18193
  基因数: 20000
  Assays: RNA, ADT
  降维: umap, wnn.umap

启动 CITEViz...
Listening on http://127.0.0.1:XXXX
```

### 2. 访问应用

**如果在本地浏览器访问云服务器**：

```bash
# 需要端口转发（SSH tunnel）
ssh -L 8888:localhost:XXXX user@your-server-ip

# 然后在本地浏览器访问：
http://localhost:8888
```

**如果云服务器有公网 IP**：

```bash
# 修改 run_citeviz.R，添加 host 参数
run_app(seurat_object = seurat_obj,
        host = "0.0.0.0",  # 允许外部访问
        port = 3838)

# 访问
http://your-server-ip:3838
```

⚠️ **安全提醒**：如果开放公网访问，建议配置防火墙或使用反向代理（nginx）。

---

## 资源监控

### 查看 R 进程资源占用

```bash
# 查看 R 进程
ps aux | grep R

# 实时监控
top -u your-username

# 查看内存使用
free -h
```

### 典型资源占用（18K 细胞）

- **内存**：~2-4 GB
- **CPU**：门控操作时占用较高
- **磁盘**：RDS 文件大小（通常几百 MB 到几 GB）

---

## 故障排除

### 问题 1：端口被占用

```
Error: Port XXXX is already in use
```

**解决**：
```r
# 指定其他端口
run_app(seurat_object = seurat_obj, port = 8888)
```

### 问题 2：找不到 ADT assay

```
Warning: ADT assay not found
```

**解决**：
```r
# 检查 assay 名称
Assays(seurat_obj)

# 如果是其他名称，重命名
seurat_obj[["ADT"]] <- seurat_obj[["Protein"]]
```

### 问题 3：内存不足

```
Error: cannot allocate vector of size XX GB
```

**解决**：
```r
# 只加载部分细胞
seurat_obj_subset <- seurat_obj[, sample(1:ncol(seurat_obj), 5000)]
run_app(seurat_object = seurat_obj_subset)
```

---

## 数据准备检查

在运行 CITEViz 前，确保数据满足要求：

```r
library(Seurat)

# 加载数据
seurat_obj <- readRDS("results/cDC1/Robjects/seurat_obj_annotated.rds")

# 1. 检查 ADT assay
if (!"ADT" %in% Assays(seurat_obj)) {
  stop("需要 ADT assay")
}

# 2. 检查 ADT 归一化
DefaultAssay(seurat_obj) <- "ADT"
if (!"data" %in% slotNames(seurat_obj@assays$ADT)) {
  # 进行 CLR 归一化
  seurat_obj <- NormalizeData(seurat_obj,
                               normalization.method = 'CLR',
                               margin = 2,
                               assay = "ADT")
}

# 3. 检查降维
if (length(seurat_obj@reductions) == 0) {
  # 运行 UMAP
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
}

# 4. 保存准备好的对象
saveRDS(seurat_obj, "seurat_obj_for_citeviz.rds")
```

---

## 总结

**关键点**：
1. ✅ 无论哪种方式，**所有计算都在服务器上**
2. ✅ 浏览器只是显示界面，不占用本地计算资源
3. ✅ 推荐直接加载服务器文件，避免不必要的上传
4. ✅ 数据和应用在同一台服务器，最高效

**命令速查**：
```bash
# 启动 CITEViz（直接加载数据）
Rscript run_citeviz.R

# SSH 端口转发（本地访问）
ssh -L 8888:localhost:XXXX user@server

# 监控资源
top
htop
```
