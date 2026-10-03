# Generic Affymetrix cohort DE - tumor vs normal, unpaired.

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(Biobase)
})
source(file.path(getwd(), "R", "utils.R"))

GSE <- Sys.getenv("GSE", "")
if (!nzchar(GSE)) stop("set GSE, e.g. GSE=GSE42568 Rscript scripts/06_de_affy.R")

NORMAL_PAT <- c(GSE42568 = "normal",
                GSE86374 = "normal tissue")
if (!GSE %in% names(NORMAL_PAT))
  stop("no normal-sample pattern defined for ", GSE)

options(timeout = 3600)
cache <- file.path(DIR_DATA, paste0(GSE, "_eset.rds"))

if (file.exists(cache)) {
  msg(GSE, ": cached eset")
  es <- readRDS(cache)
} else {
  msg(GSE, ": loading series matrix ...")
  es <- getGEO(GSE, GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = DIR_DATA)
  if (is.list(es)) es <- es[[1]]
  saveRDS(es, cache)
}
msg(GSE, ": ", nrow(es), " features x ", ncol(es), " samples")

pd <- pData(es)
ex <- exprs(es)

src <- as.character(pd$source_name_ch1)
grp <- ifelse(grepl(NORMAL_PAT[[GSE]], src, ignore.case = TRUE), "normal", "tumor")
group <- factor(grp, levels = c("normal", "tumor"))
msg(GSE, ": ", sum(group == "tumor"), " tumor vs ", sum(group == "normal"), " normal")
if (any(table(group) < 3)) stop("a group has <3 samples — check NORMAL_PAT")

if (max(ex, na.rm = TRUE) > 100) {
  msg(GSE, ": values look linear (max=", round(max(ex, na.rm = TRUE)), ") — log2 transforming")
  ex <- log2(ex + 1)
}
ex <- normalizeBetweenArrays(ex, method = "quantile")

keep_probe <- rowSums(is.na(ex)) < ncol(ex) * 0.5
msg(GSE, ": dropping ", sum(!keep_probe), " probes that are >50% NA")
ex <- ex[keep_probe, , drop = FALSE]

design <- model.matrix(~ group)
fit <- lmFit(ex, design)
fit <- eBayes(fit, trend = TRUE)

tt <- topTable(fit, coef = "grouptumor", number = Inf, sort.by = "P")
tt$probe <- rownames(tt)

fd <- fData(es)
cand <- c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Gene.symbol",
          "GeneSymbol", "gene_symbol", "SYMBOL", "Gene_Symbol", "GeneName")
sym_col <- intersect(cand, colnames(fd))

if (length(sym_col)) {
  sym <- as.character(fd[tt$probe, sym_col[1]])
  sym <- sub(" ?///.*$", "", sym)
  msg(GSE, ": symbols from fData column '", sym_col[1], "'")
} else if ("gene_assignment" %in% colnames(fd)) {
  ga <- as.character(fd[tt$probe, "gene_assignment"])
  sym <- vapply(strsplit(ga, " /// ", fixed = TRUE), function(blocks) {
    if (!length(blocks) || is.na(blocks[1])) return(NA_character_)
    parts <- trimws(strsplit(blocks[1], " // ", fixed = TRUE)[[1]])
    if (length(parts) >= 2 && nzchar(parts[2])) parts[2] else NA_character_
  }, character(1))
  sym[sym %in% c("---", "")] <- NA_character_
  msg(GSE, ": symbols parsed from 'gene_assignment' (GPL6244 format)")
} else {
  stop("no usable annotation column in fData; found: ",
       paste(colnames(fd), collapse = ", "))
}

sym[sym %in% c("---", "")] <- NA_character_
tt$symbol <- sym

res <- data.frame(symbol = tt$symbol,
                  logFC  = tt$logFC,
                  pvalue = tt$P.Value,
                  padj   = tt$adj.P.Val,
                  probe  = tt$probe,
                  stringsAsFactors = FALSE)

n_before <- nrow(res)
res <- collapse_by_symbol(res, "symbol", "logFC")
msg(GSE, ": collapsed ", n_before, " probes -> ", nrow(res), " unique gene symbols")

obj <- list(expr    = ex,
            design  = design,
            fit     = fit,
            results = tt,
            collapsed = res,
            pheno   = data.frame(sample = colnames(ex), group = group,
                                 source = src, stringsAsFactors = FALSE))
f_obj <- file.path(PROJ_ROOT, "exports", "r_objects", paste0(GSE, "_de.rds"))
saveRDS(obj, f_obj); msg("saved R objects -> ", f_obj)
f_mat <- file.path(PROJ_ROOT, "exports", "r_objects", paste0(GSE, "_expr.rds"))
saveRDS(ex, f_mat);  msg("saved expression matrix -> ", f_mat)

msg(GSE, ": significant (padj<0.05, |logFC|>1): up=",
    sum(res$padj < 0.05 & res$logFC >  1, na.rm = TRUE), " down=",
    sum(res$padj < 0.05 & res$logFC < -1, na.rm = TRUE))

write_de(res, GSE)
