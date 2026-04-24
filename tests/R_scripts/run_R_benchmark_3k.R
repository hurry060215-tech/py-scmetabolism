# ---------------------------------------------------
# 1. 环境初始化（解决警告报错与版本冲突）
# ---------------------------------------------------
.libPaths(unique(c("C:/Users/17904/Documents/R/win-library/4.5", .libPaths())))

# 关键：将警告调回默认，防止 Matrix 包的过时警告中断运行
options(warn = 0)
# 关键：强制 Seurat 以后只创建老版本 Assay (V3 风格)
options(Seurat.object.assay.version = "v3")

library(Seurat)
library(scMetabolism)
library(GSVA)
library(AUCell)

set.seed(42)

# ---------------------------------------------------
# 2. 读取数据 (表达矩阵 + GMT)
# ---------------------------------------------------
cat("\n--- 正在加载数据 ---\n")
# 读取表达矩阵
exprs <- as.matrix(read.csv("D:/test/py-scMetabolism/py-scmetabolism/data/adipocyte/exprs_3k.csv", row.names=1))
# 读取基因集
gmt_file <- "C:/Users/17904/Documents/R/win-library/4.5/scMetabolism/data/KEGG_metabolism_nc.gmt"
pathway_lines <- readLines(gmt_file)
pathway_list <- list()
for (line in pathway_lines) {
  parts <- strsplit(line, "\t")[[1]]
  pathway_list[[parts[1]]] <- parts[3:length(parts)][parts[3:length(parts)] != "" & parts[3:length(parts)] != "NA"]
}

# 创建 Seurat 对象 (不进行 NormalizeData，与源码一致)
countexp.Seurat <- CreateSeuratObject(counts = exprs)

time_results <- list()

# ---------------------------------------------------
# 3. 官方函数运行与"双中括号"提取 (避开 GetAssayData)
# ---------------------------------------------------

# --- [1] 官方 VISION ---
cat("\nRunning Official VISION ...\n")
t0 <- Sys.time()
res_vision_obj <- suppressWarnings(
  sc.metabolism.Seurat(obj = countexp.Seurat, method = "VISION",
                       imputation = FALSE, ncores = 1, metabolism.type = "KEGG")
)
time_results$VISION <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# 【关键修复】：从 score 槽中提取结果，这是 scMetabolism 包的标准存储位置
if ("METABOLISM" %in% names(res_vision_obj@assays)) {
  # 与源码一致：从 score 槽中提取结果
  res_mat_v <- as.matrix(res_vision_obj[["METABOLISM"]]$score)
  write.csv(res_mat_v, "metabolism_res_VISION.csv", quote = TRUE)
  cat("VISION 保存成功，维度:", nrow(res_mat_v), "x", ncol(res_mat_v), "\n")
} else {
  cat("错误：未在对象中找到 METABOLISM Assay。\n")
}

# --- [2] 官方 AUCell ---
cat("\nRunning Official AUCell ...\n")
t0 <- Sys.time()
res_aucell_obj <- sc.metabolism.Seurat(obj = countexp.Seurat, method = "AUCell",
                                       imputation = FALSE, ncores = 1, metabolism.type = "KEGG")
time_results$AUCell <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

if ("METABOLISM" %in% names(res_aucell_obj@assays)) {
  # 同理从 score 槽中提取 AUCell 结果
  res_mat_a <- as.matrix(res_aucell_obj[["METABOLISM"]]$score)
  write.csv(res_mat_a, "metabolism_res_AUCell.csv", quote = TRUE)
  cat("AUCell 保存成功。\n")
}

# --- [3] 官方 ssGSEA ---
cat("\nRunning Official ssGSEA ...\n")
t0 <- Sys.time()
res_ssgsea_obj <- sc.metabolism.Seurat(obj = countexp.Seurat, method = "ssGSEA",
                                       imputation = FALSE, ncores = 1, metabolism.type = "KEGG")
time_results$ssGSEA <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

if ("METABOLISM" %in% names(res_ssgsea_obj@assays)) {
  # 从 score 槽中提取 ssGSEA 结果
  res_mat_s <- as.matrix(res_ssgsea_obj[["METABOLISM"]]$score)
  write.csv(res_mat_s, "metabolism_res_ssGSEA.csv", quote = TRUE)
  cat("ssGSEA 保存成功。\n")
}

# --- [4] 官方 GSVA ---
cat("\nRunning Official GSVA (Poisson kcdf) ...\n")
t0 <- Sys.time()
res_gsva_obj <- sc.metabolism.Seurat(obj = countexp.Seurat, method = "GSVA",
                                     imputation = FALSE, ncores = 1, metabolism.type = "KEGG")
time_results$GSVA <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

if ("METABOLISM" %in% names(res_gsva_obj@assays)) {
  # 从 score 槽中提取 GSVA 结果
  res_mat_g <- as.matrix(res_gsva_obj[["METABOLISM"]]$score)
  write.csv(res_mat_g, "metabolism_res_GSVA.csv", quote = TRUE)
  cat("GSVA 保存成功。\n")
}

# ---------------------------------------------------
# 4. 官方函数运行完成
# ---------------------------------------------------
# 已通过官方函数完成了所有方法的运行：VISION、AUCell、ssGSEA 和 GSVA
# 结果已保存到相应的 CSV 文件中

# ---------------------------------------------------
# 5. 最终统计汇总
# ---------------------------------------------------
cat("\n" , rep("=", 40), "\n")
cat("R 执行时间汇总 (单核):\n")
cat(sprintf("  VISION:  %.1fs\n", time_results$VISION))
cat(sprintf("  AUCell:  %.1fs\n", time_results$AUCell))
cat(sprintf("  ssGSEA:  %.1fs\n", time_results$ssGSEA))
cat(sprintf("  GSVA:    %.1fs\n", time_results$GSVA))
cat(rep("=", 40), "\n")
