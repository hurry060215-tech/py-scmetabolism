# ---------------------------------------------------
# 重新生成 VISION/AUCell 结果 (使用 quote=TRUE)
# ---------------------------------------------------
.libPaths(unique(c("C:/Users/17904/Documents/R/win-library/4.5", .libPaths())))
options(warn = 0)
options(Seurat.object.assay.version = "v3")

library(Seurat)
library(scMetabolism)
library(AUCell)

set.seed(42)

# 读取数据
cat("\n--- 正在加载数据 ---\n")
exprs <- as.matrix(read.csv("D:/test/py-scMetabolism/py-scmetabolism/data/adipocyte/exprs_3k.csv", row.names=1))

# 创建 Seurat 对象
countexp.Seurat <- CreateSeuratObject(counts = exprs)
countexp.Seurat <- NormalizeData(countexp.Seurat, verbose = FALSE)

time_results <- list()

# --- VISION ---
cat("\nRunning Official VISION ...\n")
t0 <- Sys.time()
res_vision_obj <- suppressWarnings(
  sc.metabolism.Seurat(obj = countexp.Seurat, method = "VISION",
                       imputation = FALSE, ncores = 1, metabolism.type = "KEGG")
)
time_results$VISION <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

if ("METABOLISM" %in% names(res_vision_obj@assays)) {
  res_mat_v <- as.matrix(res_vision_obj[["METABOLISM"]][["data"]])
  write.csv(res_mat_v, "D:/test/py-scMetabolism/py-scmetabolism/examples/r/metabolism_res_VISION.csv", quote = TRUE)
  cat("VISION 保存成功，维度:", nrow(res_mat_v), "x", ncol(res_mat_v), "\n")
} else {
  cat("错误：未在对象中找到 METABOLISM Assay。\n")
}

# --- AUCell ---
cat("\nRunning Official AUCell ...\n")
t0 <- Sys.time()
res_aucell_obj <- sc.metabolism.Seurat(obj = countexp.Seurat, method = "AUCell",
                                       imputation = FALSE, ncores = 1, metabolism.type = "KEGG")
time_results$AUCell <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

if ("METABOLISM" %in% names(res_aucell_obj@assays)) {
  res_mat_a <- as.matrix(res_aucell_obj[["METABOLISM"]][["data"]])
  write.csv(res_mat_a, "D:/test/py-scMetabolism/py-scmetabolism/examples/r/metabolism_res_AUCell.csv", quote = TRUE)
  cat("AUCell 保存成功，维度:", nrow(res_mat_a), "x", ncol(res_mat_a), "\n")
}

# 汇总
cat("\n", rep("=", 40), "\n")
cat("R 执行时间汇总 (单核):\n")
cat(sprintf("  VISION:  %.1fs\n", time_results$VISION))
cat(sprintf("  AUCell:  %.1fs\n", time_results$AUCell))
cat(rep("=", 40), "\n")
