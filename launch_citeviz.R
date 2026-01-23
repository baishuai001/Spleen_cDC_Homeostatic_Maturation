#!/usr/bin/env Rscript
# 简单的 CITEViz 启动脚本

library(CITEViz)

cat("=================================================\n")
cat("启动 CITEViz\n")
cat("=================================================\n\n")

cat("提示:\n")
cat("1. CITEViz 将在浏览器中打开\n")
cat("2. 点击 'Browse...' 上传数据文件\n")
cat("3. 推荐使用预处理的文件:\n")
cat("   - citeviz_fast_5k.rds (快速，推荐)\n")
cat("   - citeviz_medium_8k.rds (中等)\n")
cat("   - citeviz_full_compatible.rds (完整)\n\n")

cat("如果还没有准备数据文件，先运行:\n")
cat("  source('prepare_citeviz_data.R')\n\n")

cat("启动 CITEViz...\n")
cat("=================================================\n\n")

# 尝试不同的启动函数
if (exists("run_app", where = "package:CITEViz")) {
  run_app()
} else if (exists("launchApp", where = "package:CITEViz")) {
  launchApp()
} else if (exists("CITEViz", where = "package:CITEViz")) {
  CITEViz()
} else {
  cat("错误：未找到 CITEViz 启动函数\n")
  cat("可用函数:\n")
  print(ls("package:CITEViz"))
  stop("请检查 CITEViz 包是否正确安装")
}
