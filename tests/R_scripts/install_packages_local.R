#!/usr/bin/env Rscript
# Install required packages to local directory

# Set local library path
local_lib <- "R_libs"
if (!dir.exists(local_lib)) {
  dir.create(local_lib, recursive = TRUE)
}

# Add to library paths
.libPaths(c(local_lib, .libPaths()))
cat("Library paths:\n")
print(.libPaths())

# CRAN packages to install
cran_packages <- c("AUCell", "GSVA", "GSEABase", "devtools", "remotes")

cat("\nInstalling CRAN packages...\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE, lib.loc = local_lib)) {
    cat(paste("Installing", pkg, "...\n"))
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org", lib = local_lib)
      cat(paste("  Successfully installed", pkg, "\n"))
    }, error = function(e) {
      cat(paste("  Failed to install", pkg, ":", conditionMessage(e), "\n"))
    })
  } else {
    cat(paste(pkg, "already installed.\n"))
  }
}

# Check if VISION is available (might be from Bioconductor)
cat("\nChecking for VISION package...\n")
if (!requireNamespace("VISION", quietly = TRUE, lib.loc = local_lib)) {
  cat("VISION not found. Trying to install...\n")
  # VISION might be on Bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = local_lib)
  }
  library(BiocManager, lib.loc = local_lib)
  BiocManager::install("VISION", lib = local_lib, update = FALSE, ask = FALSE)
  cat("VISION installed from Bioconductor.\n")
} else {
  cat("VISION already installed.\n")
}

# Install scMetabolism from GitHub
cat("\nInstalling scMetabolism from GitHub...\n")
if (!requireNamespace("scMetabolism", quietly = TRUE, lib.loc = local_lib)) {
  # Try with devtools
  if (requireNamespace("devtools", quietly = TRUE, lib.loc = local_lib)) {
    library(devtools, lib.loc = local_lib)
    tryCatch({
      devtools::install_github("YingchengWu/scMetabolism", lib = local_lib)
      cat("scMetabolism installed from GitHub.\n")
    }, error = function(e) {
      cat(paste("Failed to install scMetabolism from GitHub:", conditionMessage(e), "\n"))
      # Try from local directory
      cat("Trying to install from local directory...\n")
      scmetabolism_dir <- "D:/test/py-scMetabolism/scMetabolism-main"
      if (dir.exists(scmetabolism_dir)) {
        devtools::install_local(scmetabolism_dir, lib = local_lib)
        cat("scMetabolism installed from local directory.\n")
      } else {
        cat("Local scMetabolism directory not found.\n")
      }
    })
  } else {
    cat("devtools not available for installing scMetabolism.\n")
  }
} else {
  cat("scMetabolism already installed.\n")
}

# Verify installation
cat("\nVerifying installation...\n")
required <- c("scMetabolism", "AUCell", "GSVA", "GSEABase", "VISION")
for (pkg in required) {
  if (requireNamespace(pkg, quietly = TRUE, lib.loc = local_lib)) {
    cat(paste("✓", pkg, "is available.\n"))
  } else {
    cat(paste("✗", pkg, "is NOT available.\n"))
  }
}