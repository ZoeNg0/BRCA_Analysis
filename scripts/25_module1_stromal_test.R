# Module 1: tumour EMT vs stromal content

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi)
})
source(file.path(getwd(), "R", "utils.R"))

DIR_NET <- file.path(DIR_OUT, "network")
DIR_EN  <- file.path(DIR_OUT, "hub_enrichment")
net <- readRDS(file.path(DIR_NET, "tcga_network.rds"))
mod <- as.data.table(net$modules)

counts <- readRDS(file.path(PROJ_ROOT, "exports", "r_objects", "rna.rds"))
counts <- counts[, substr(colnames(counts), 14, 15) == "01", drop = FALSE]
counts <- counts[, !duplicated(substr(colnames(counts), 1, 12)), drop = FALSE]
dge  <- calcNormFactors(DGEList(counts = counts), method = "TMM")
lcpm <- cpm(dge, log = TRUE, prior.count = 1)
rownames(lcpm) <- strip_ens_version(rownames(lcpm))
rm(counts, dge); invisible(gc())

s2e <- as.data.table(suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = rownames(lcpm), keytype = "ENSEMBL", columns = "SYMBOL")))
s2e <- s2e[!is.na(SYMBOL)][!duplicated(SYMBOL)]
get_expr <- function(sym) {
  e <- s2e[SYMBOL %in% sym]
  m <- lcpm[e$ENSEMBL, , drop = FALSE]; rownames(m) <- e$SYMBOL; m
}

MARK <- list(
  fibroblast  = c("PDGFRB","THY1","FAP","COL1A1","LUM","DCN"),
  epithelial  = c("EPCAM","KRT8","KRT18","KRT19","CDH1","ESR1"),
  endothelial = c("PECAM1","VWF","CDH5"),
  immune      = c("PTPRC","CD3E","CD68"))

score_module <- function(m) {
  g <- mod[module == m]$symbol
  X <- get_expr(g)
  X <- X[rowSums(is.na(X)) == 0, , drop = FALSE]
  X <- t(scale(t(X)))
  p <- prcomp(t(X[complete.cases(X), , drop = FALSE]), center = TRUE, scale. = FALSE)
  s <- p$x[, 1]
  if (cor(s, colMeans(X, na.rm = TRUE)) < 0) s <- -s
  list(score = s, var = summary(p)$importance[2, 1], n = nrow(X))
}

rows <- list()
for (m in c(1, 2, 3, 4)) {
  sc <- score_module(m)
  msg("module ", m, ": ", sc$n, " genes, PC1 explains ",
      round(100 * sc$var), "% of variance")
  for (ct in names(MARK)) {
    M <- get_expr(MARK[[ct]])
    r <- apply(M, 1, function(v) cor(v, sc$score, method = "spearman"))
    rows[[paste(m, ct)]] <- data.table(module = m, compartment = ct,
      marker = names(r), rho = as.numeric(r))
  }
}
R <- rbindlist(rows)
fwrite(R, file.path(DIR_EN, "module_compartment_correlation.csv"))

msg("\n=== Spearman rho: module score vs cell-type markers (1,095 TCGA tumours) ===")
tab <- dcast(R[, .(rho = round(mean(rho), 2)), by = .(module, compartment)],
             module ~ compartment, value.var = "rho")
print(as.data.frame(tab))

msg("\n=== module 1, marker by marker ===")
print(as.data.frame(R[module == 1][order(compartment, -rho),
      .(compartment, marker, rho = round(rho, 2))]))
