# GSE58135 - differential expression, tumor vs normal breast tissue.

suppressPackageStartupMessages({
  library(limma)
  library(data.table)
})
source(file.path(getwd(), "R", "utils.R"))

NORMAL_SET <- Sys.getenv("NORMAL_SET", "all")

meta <- fread(file.path(DIR_DATA, "GSE58135_sample_metadata.csv"))
supp_dir <- file.path(DIR_DATA, "GSE58135_fpkm")

meta[, file := file.path(supp_dir, basename(fpkm))]
missing <- meta[!file.exists(file)]
if (nrow(missing)) {
  warning(nrow(missing), " FPKM files missing; dropping them")
  meta <- meta[file.exists(file)]
}
msg("GSE58135: reading ", nrow(meta), " FPKM files ...")

read_fpkm <- function(f) {
  d <- fread(f, select = c("tracking_id", "gene_short_name", "FPKM"),
             showProgress = FALSE)
  setorder(d, -FPKM)
  d[!duplicated(tracking_id)]
}

first <- read_fpkm(meta$file[1])
genes <- data.table(tracking_id = first$tracking_id,
                    symbol      = first$gene_short_name)

mat <- matrix(NA_real_, nrow = nrow(genes), ncol = nrow(meta),
              dimnames = list(genes$tracking_id, meta$gsm))

for (i in seq_len(nrow(meta))) {
  d <- read_fpkm(meta$file[i])
  mat[, i] <- d$FPKM[match(rownames(mat), d$tracking_id)]
  if (i %% 25 == 0) msg("  ", i, "/", nrow(meta))
}
msg("matrix: ", nrow(mat), " genes x ", ncol(mat), " samples")

grp <- meta$group
if (NORMAL_SET == "healthy") {
  keep_s <- grp %in% c("tumor", "healthy_normal")
} else {
  keep_s <- grp %in% c("tumor", "adjacent_normal", "healthy_normal")
}
mat  <- mat[, keep_s, drop = FALSE]
meta <- meta[keep_s]
group <- factor(ifelse(meta$group == "tumor", "tumor", "normal"),
                levels = c("normal", "tumor"))
msg("contrast: ", sum(group == "tumor"), " tumor vs ",
    sum(group == "normal"), " normal  (NORMAL_SET=", NORMAL_SET, ")")

mat[is.na(mat)] <- 0
min_grp <- min(table(group))
keep_g  <- rowSums(mat > 1) >= min_grp
msg("expression filter (FPKM>1 in >=", min_grp, " samples): ",
    sum(keep_g), "/", length(keep_g), " genes retained")
mat <- mat[keep_g, , drop = FALSE]

lg <- log2(mat + 1)

lg <- normalizeBetweenArrays(lg, method = "quantile")

design <- model.matrix(~ group)
fit <- lmFit(lg, design)
fit <- eBayes(fit, trend = TRUE)

tt <- topTable(fit, coef = "grouptumor", number = Inf, sort.by = "P")
tt$tracking_id <- rownames(tt)
tt$symbol <- genes$symbol[match(strip_ens_version(tt$tracking_id),
                                strip_ens_version(genes$tracking_id))]

res <- data.frame(symbol = tt$symbol,
                  logFC  = tt$logFC,
                  pvalue = tt$P.Value,
                  padj   = tt$adj.P.Val,
                  ensembl = strip_ens_version(tt$tracking_id),
                  stringsAsFactors = FALSE)

n_before <- nrow(res)
res <- collapse_by_symbol(res, "symbol", "logFC")
msg("collapsed ", n_before, " -> ", nrow(res), " unique gene symbols")

obj <- list(fpkm    = mat,
            expr    = lg,
            design  = design,
            fit     = fit,
            results = tt,
            collapsed = res,
            pheno   = data.frame(sample = colnames(lg), group = meta$group,
                                 subtype = meta$subtype, contrast = group,
                                 stringsAsFactors = FALSE))
suffix <- if (NORMAL_SET == "healthy") "_healthy" else ""
f_obj <- file.path(PROJ_ROOT, "exports", "r_objects",
                   paste0("GSE58135_de", suffix, ".rds"))
saveRDS(obj, f_obj)
msg("saved R objects -> ", f_obj)

f_mat <- file.path(PROJ_ROOT, "exports", "r_objects",
                   paste0("GSE58135_expr", suffix, ".rds"))
saveRDS(lg, f_mat)
msg("saved expression matrix -> ", f_mat)

msg("significant (padj<0.05, |logFC|>1): up=",
    sum(res$padj < 0.05 & res$logFC >  1, na.rm = TRUE), " down=",
    sum(res$padj < 0.05 & res$logFC < -1, na.rm = TRUE))

cohort <- if (NORMAL_SET == "healthy") "GSE58135_healthy" else "GSE58135"
write_de(res, cohort)
