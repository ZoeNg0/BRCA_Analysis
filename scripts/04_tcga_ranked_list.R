# TCGA-BRCA DESeq2 results mapped to gene symbols

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
source(file.path(getwd(), "R", "utils.R"))

dea_file <- file.path(PROJ_ROOT, "exports", "tables", "DEA.csv")
dea <- fread(dea_file)
setnames(dea, 1, "gene_id")
msg("TCGA: ", nrow(dea), " genes from ", basename(dea_file))

dea[, ensembl := strip_ens_version(gene_id)]

sym <- suppressMessages(
  AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(dea$ensembl),
                        column = "SYMBOL", keytype = "ENSEMBL",
                        multiVals = "first"))
dea[, symbol := sym[ensembl]]
msg("mapped via org.Hs.eg.db: ", sum(!is.na(dea$symbol)), "/", nrow(dea))

fb_file <- list.files(file.path(DIR_DATA, "GSE58135_fpkm"),
                      pattern = "fpkm_tracking", full.names = TRUE)[1]
if (!is.na(fb_file) && any(is.na(dea$symbol))) {
  fb <- unique(fread(fb_file, select = c("tracking_id", "gene_short_name"),
                     showProgress = FALSE))
  fb[, ensembl := strip_ens_version(tracking_id)]
  fb <- fb[!duplicated(ensembl) & gene_short_name != "-"]
  hit <- is.na(dea$symbol)
  dea[hit, symbol := fb$gene_short_name[match(ensembl, fb$ensembl)]]
  msg("fallback recovered ", sum(hit & !is.na(dea$symbol)), " more symbols")
}

res <- data.frame(symbol = dea$symbol,
                  logFC  = dea$log2FoldChange,
                  pvalue = dea$pvalue,
                  padj   = dea$padj,
                  ensembl = dea$ensembl,
                  stringsAsFactors = FALSE)

res <- res[!is.na(res$pvalue) & !is.na(res$logFC), ]
n_before <- nrow(res)
res <- collapse_by_symbol(res, "symbol", "logFC")
msg("collapsed ", n_before, " -> ", nrow(res), " unique gene symbols")

msg("significant (padj<0.05, |logFC|>1): up=",
    sum(res$padj < 0.05 & res$logFC >  1, na.rm = TRUE), " down=",
    sum(res$padj < 0.05 & res$logFC < -1, na.rm = TRUE))

write_de(res, "TCGA")
