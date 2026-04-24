# ---------------------------------------------------
# 重新生成 ssGSEA/GSVA 结果 (使用标准流程)
# ---------------------------------------------------
.libPaths(unique(c("C:/Users/17904/Documents/R/win-library/4.5", .libPaths())))
options(warn = 0)

library(GSVA)

set.seed(42)

# 读取数据
cat("\n--- 正在加载数据 ---\n")
exprs <- as.matrix(read.csv("D:/test/py-scMetabolism/py-scmetabolism/data/adipocyte/exprs_3k.csv", row.names=1))

# 读取基因集
gmt_file <- "C:/Users/17904/Documents/R/win-library/4.5/scMetabolism/data/KEGG_metabolism_nc.gmt"
pathway_lines <- readLines(gmt_file)
pathway_list <- list()
for (line in pathway_lines) {
  parts <- strsplit(line, "\t")[[1]]
  pathway_list[[parts[1]]] <- parts[3:length(parts)][parts[3:length(parts)] != "" & parts[3:length(parts)] != "NA"]
}

time_results <- list()

# --- ssGSEA (ssgseaParam: 纯秩方法，无kcdf参数) ---
cat("\nRunning ssGSEA ...\n")
t0 <- Sys.time()
param_ss <- ssgseaParam(exprs, pathway_list, normalize=TRUE)
res_ssgsea <- gsva(param_ss)
time_results$ssGSEA <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
write.csv(as.matrix(res_ssgsea), "D:/test/py-scMetabolism/py-scmetabolism/examples/r/metabolism_res_ssGSEA.csv", quote = TRUE)
cat("ssGSEA done:", nrow(res_ssgsea), "x", ncol(res_ssgsea), "\n")

# --- GSVA (gsvaParam: kcdf="Poisson" 与 scMetabolism 原始包一致) ---
cat("\nRunning GSVA (Poisson kcdf) ...\n")
t0 <- Sys.time()
param_gv <- gsvaParam(exprs, pathway_list, kcdf="Poisson")
res_gsva <- gsva(param_gv)
time_results$GSVA <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
write.csv(as.matrix(res_gsva), "D:/test/py-scMetabolism/py-scmetabolism/examples/r/metabolism_res_GSVA.csv", quote = TRUE)
cat("GSVA done:", nrow(res_gsva), "x", ncol(res_gsva), "\n")

# 汇总
cat("\n", rep("=", 40), "\n")
cat("R 执行时间汇总 (单核):\n")
cat(sprintf("  ssGSEA:  %.1fs\n", time_results$ssGSEA))
cat(sprintf("  GSVA:    %.1fs\n", time_results$GSVA))
cat(rep("=", 40), "\n")
