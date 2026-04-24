#!/usr/bin/env Rscript
# Install scMetabolism from local directory
cat("Installing scMetabolism from local directory...\n")

# First check if devtools is available
if (!requireNamespace("devtools", quietly = TRUE)) {
  cat("Installing devtools...\n")
  install.packages("devtools", repos = "https://cloud.r-project.org", lib = Sys.getenv("R_LIBS_USER"))
}

library(devtools)

# Install from local directory
scmetabolism_path <- "D:/test/py-scMetabolism/scMetabolism-main"
cat(paste("Installing from:", scmetabolism_path, "\n"))

# Install dependencies first
cat("Installing dependencies...\n")
deps <- c("AUCell", "GSVA", "GSEABase", "VISION")
for (dep in deps) {
  if (!requireNamespace(dep, quietly = TRUE)) {
    cat(paste("Installing", dep, "...\n"))
    install.packages(dep, repos = "https://cloud.r-project.org", lib = Sys.getenv("R_LIBS_USER"))
  }
}

# Now install scMetabolism
install_local(scmetabolism_path)

# Check installation
if (requireNamespace("scMetabolism", quietly = TRUE)) {
  cat("scMetabolism installed successfully!\n")
  library(scMetabolism)
  cat(paste("Version:", packageVersion("scMetabolism"), "\n"))
} else {
  cat("Failed to install scMetabolism.\n")
}