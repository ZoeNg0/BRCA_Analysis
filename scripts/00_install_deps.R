# Install dependencies

repos <- "https://cloud.r-project.org"

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = repos)
}

cran_pkgs <- c("RobustRankAggreg", "data.table", "R.utils", "matrixStats", "dplyr",
               "ggplot2", "jsonlite", "survival", "survminer", "igraph", "msigdbr", "png")
bioc_pkgs <- c("GEOquery", "limma", "edgeR", "sva", "DESeq2", "TCGAbiolinks",
               "SummarizedExperiment", "maftools", "EnhancedVolcano", "AnnotationDbi",
               "org.Hs.eg.db", "clusterProfiler", "ReactomePA", "DOSE", "enrichplot")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message(">>> installing (CRAN): ", p)
    install.packages(p, repos = repos)
  } else message("--- already present: ", p)
}

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message(">>> installing (Bioc): ", p)
    BiocManager::install(p, ask = FALSE, update = FALSE)
  } else message("--- already present: ", p)
}

message("\n=== final status ===")
for (p in c(cran_pkgs, bioc_pkgs)) {
  message(sprintf("%-18s %s", p,
    if (requireNamespace(p, quietly = TRUE)) "OK" else "MISSING"))
}
