# GSE70947 - paired differential expression, tumor vs matched adjacent normal.

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(Biobase)
})
source(file.path(getwd(), "R", "utils.R"))

es <- readRDS(file.path(DIR_DATA, "GSE70947_eset.rds"))
msg("GSE70947: ", nrow(es), " features x ", ncol(es), " samples")

pd <- pData(es)
ex <- exprs(es)

titles  <- as.character(pd$title)
patient <- sub("-.*$", "", titles)
tissue  <- sub("^[^-]*-", "", titles)
tissue  <- ifelse(grepl("normal", tissue, ignore.case = TRUE), "normal", "tumor")

stopifnot(length(unique(patient)) == 148)
tab <- table(patient, tissue)
if (!all(tab == 1)) stop("pairing is not 1:1 — inspect sample titles")
msg("pairing verified: ", length(unique(patient)), " patients x 2 tissues")

tissue  <- factor(tissue, levels = c("normal", "tumor"))
patient <- factor(patient)

keep_probe <- rowSums(is.na(ex)) < ncol(ex) * 0.5
msg("dropping ", sum(!keep_probe), " probes that are >50% NA")
ex <- ex[keep_probe, , drop = FALSE]

design <- model.matrix(~ patient + tissue)
msg("design: ", nrow(design), " samples x ", ncol(design), " coefficients")

fit <- lmFit(ex, design)
fit <- eBayes(fit, trend = TRUE)

tt <- topTable(fit, coef = "tissuetumor", number = Inf, sort.by = "P")
tt$probe <- rownames(tt)

fd <- fData(es)
cand <- c("GENE_SYMBOL", "Gene symbol", "Gene.symbol", "GeneSymbol",
          "gene_symbol", "SYMBOL", "GENE_NAME", "Gene_Symbol",
          "GeneName")
sym_col <- intersect(cand, colnames(fd))
if (!length(sym_col)) {
  score <- sapply(fd, function(x) mean(grepl("^[A-Z][A-Z0-9-]{1,14}$",
                                             as.character(x)), na.rm = TRUE))
  sym_col <- names(which.max(score))
  msg("no standard symbol column; inferred '", sym_col,
      "' (", round(100 * max(score)), "% HGNC-like)")
} else {
  sym_col <- sym_col[1]
  msg("gene symbols from fData column '", sym_col, "'")
}

tt$symbol <- as.character(fd[tt$probe, sym_col])

res <- data.frame(symbol = tt$symbol,
                  logFC  = tt$logFC,
                  pvalue = tt$P.Value,
                  padj   = tt$adj.P.Val,
                  probe  = tt$probe,
                  stringsAsFactors = FALSE)

n_before <- nrow(res)
res <- collapse_by_symbol(res, "symbol", "logFC")
msg("collapsed ", n_before, " probes -> ", nrow(res), " unique gene symbols")

obj <- list(expr    = ex,
            design  = design,
            fit     = fit,
            results = tt,
            collapsed = res,
            pheno   = data.frame(sample = colnames(ex), patient = patient,
                                 tissue = tissue, stringsAsFactors = FALSE))
f_obj <- file.path(PROJ_ROOT, "exports", "r_objects", "GSE70947_de.rds")
saveRDS(obj, f_obj)
msg("saved R objects -> ", f_obj)

f_mat <- file.path(PROJ_ROOT, "exports", "r_objects", "GSE70947_expr.rds")
saveRDS(ex, f_mat)
msg("saved expression matrix -> ", f_mat)

msg("significant (padj<0.05, |logFC|>1): up=",
    sum(res$padj < 0.05 & res$logFC >  1, na.rm = TRUE), " down=",
    sum(res$padj < 0.05 & res$logFC < -1, na.rm = TRUE))

write_de(res, "GSE70947")
