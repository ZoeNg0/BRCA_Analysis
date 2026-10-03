# TCGA co-expression network and hub genes

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(igraph)
})
source(file.path(getwd(), "R", "utils.R"))

PADJ_MAX <- 0.05
LFC_MIN  <- 1.0
RHO_MIN  <- 0.70
N_HUB    <- 200
BLOCK    <- 500
EXPR_CPM <- 0.5
EXPR_PCT <- 0.10

DIR_NET <- file.path(DIR_OUT, "network")
dir.create(DIR_NET, recursive = TRUE, showWarnings = FALSE)

de <- as.data.table(read_de("TCGA"))
deg <- de[!is.na(padj) & padj < PADJ_MAX & abs(logFC) > LFC_MIN]
msg("TCGA DEGs (padj<", PADJ_MAX, ", |logFC|>", LFC_MIN, "): ", nrow(deg),
    "  (up ", sum(deg$logFC > 0), " / down ", sum(deg$logFC < 0), ")")

counts <- readRDS(file.path(PROJ_ROOT, "exports", "r_objects", "rna.rds"))
counts <- counts[, substr(colnames(counts), 14, 15) == "01", drop = FALSE]
counts <- counts[, !duplicated(substr(colnames(counts), 1, 12)), drop = FALSE]
dge  <- calcNormFactors(DGEList(counts = counts), method = "TMM")
cpm_lin <- cpm(dge)
lcpm <- cpm(dge, log = TRUE, prior.count = 1)
rownames(lcpm) <- strip_ens_version(rownames(lcpm))
rownames(cpm_lin) <- rownames(lcpm)
keep <- rowMeans(cpm_lin >= EXPR_CPM) >= EXPR_PCT
lcpm <- lcpm[keep, , drop = FALSE]
rm(counts, dge, cpm_lin); invisible(gc())
msg("TCGA tumours: ", ncol(lcpm), " patients x ", nrow(lcpm), " expressed genes")

suppressPackageStartupMessages({library(org.Hs.eg.db); library(AnnotationDbi)})
map <- as.data.table(suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = deg$symbol, keytype = "SYMBOL",
  columns = c("ENSEMBL", "ENTREZID"))))
map <- map[!is.na(ENSEMBL) & ENSEMBL %in% rownames(lcpm)][!duplicated(SYMBOL)]
map <- map[!duplicated(ENSEMBL)]
msg("DEGs with an expressed ENSEMBL id: ", nrow(map))

expr_map <- as.data.table(suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = rownames(lcpm), keytype = "ENSEMBL", columns = "SYMBOL")))
expressed_symbols <- unique(expr_map[!is.na(SYMBOL)]$SYMBOL)
saveRDS(expressed_symbols, file.path(DIR_NET, "expressed_symbols.rds"))
msg("expressed-gene background: ", length(expressed_symbols), " symbols")

X <- lcpm[map$ENSEMBL, , drop = FALSE]
rownames(X) <- map$SYMBOL
rm(lcpm); invisible(gc())

msg("ranking ", nrow(X), " genes across ", ncol(X), " tumours ...")
R <- t(apply(X, 1, rank))
R <- t(scale(t(R)))
R[is.na(R)] <- 0
denom <- ncol(R) - 1
n <- nrow(R)

msg("building co-expression edges at |rho| >= ", RHO_MIN, " ...")
edges <- vector("list", ceiling(n / BLOCK))
k <- 0L
for (start in seq(1, n, by = BLOCK)) {
  idx <- start:min(start + BLOCK - 1L, n)
  cb  <- tcrossprod(R[idx, , drop = FALSE], R) / denom
  for (ii in seq_along(idx)) {
    gi <- idx[ii]
    if (gi == n) next
    v  <- cb[ii, (gi + 1L):n]
    hit <- which(abs(v) >= RHO_MIN)
    if (length(hit)) {
      k <- k + 1L
      edges[[k]] <- data.table(from = gi, to = gi + hit, rho = v[hit])
    }
  }
  rm(cb)
  if (start %% (BLOCK * 4) == 1) msg("  ", min(idx), "/", n, " genes scanned")
}
E <- rbindlist(edges[seq_len(k)])
rm(edges); invisible(gc())
msg("edges retained: ", nrow(E), "  (of ", format(choose(n, 2), big.mark = ","),
    " possible pairs)")

E[, `:=`(gene1 = rownames(R)[from], gene2 = rownames(R)[to])]

g <- graph_from_data_frame(E[, .(gene1, gene2)], directed = FALSE,
                           vertices = data.frame(name = rownames(R)))
E(g)$weight <- abs(E$rho)
deg_v <- degree(g)
msg("graph: ", vcount(g), " nodes, ", ecount(g), " edges; ",
    sum(deg_v == 0), " isolated nodes")

gc_sub <- induced_subgraph(g, which(deg_v > 0))
set.seed(1)
cl <- cluster_louvain(gc_sub, weights = E(gc_sub)$weight)
mod <- data.table(symbol = V(gc_sub)$name, module = as.integer(membership(cl)))
msg("modules: ", length(unique(mod$module)),
    "  (largest ", max(table(mod$module)), " genes)")

hub <- data.table(symbol = names(deg_v), degree = as.integer(deg_v))
hub <- merge(hub, mod, by = "symbol", all.x = TRUE)
hub <- merge(hub, deg[, .(symbol, logFC, padj)], by = "symbol", all.x = TRUE)
setorder(hub, -degree)
hub[, rank := seq_len(.N)]
hub[, is_hub := rank <= N_HUB]

fwrite(hub, file.path(DIR_NET, "tcga_deg_network_degree.csv"))
fwrite(E[, .(gene1, gene2, rho)], file.path(DIR_NET, "tcga_deg_network_edges.csv.gz"))
fwrite(hub[is_hub == TRUE], file.path(DIR_NET, "tcga_hub_genes.csv"))
saveRDS(list(hub = hub, modules = mod, params =
             list(PADJ_MAX = PADJ_MAX, LFC_MIN = LFC_MIN, RHO_MIN = RHO_MIN,
                  N_HUB = N_HUB, n_tumours = ncol(X))),
        file.path(DIR_NET, "tcga_network.rds"))

msg("\n=== top 30 TCGA hub genes by degree ===")
print(as.data.frame(head(hub[, .(rank, symbol, degree, module, logFC, padj)], 30)))

msg("\n=== module sizes among the ", N_HUB, " hubs ===")
print(hub[is_hub == TRUE, .N, by = module][order(-N)])
