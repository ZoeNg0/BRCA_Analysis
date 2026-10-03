# EMT module hub genes

suppressPackageStartupMessages({
  library(data.table); library(igraph); library(clusterProfiler); library(ReactomePA)
  library(DOSE); library(msigdbr); library(org.Hs.eg.db); library(AnnotationDbi)
})
source(file.path(getwd(), "R", "utils.R"))

DIR_NET <- file.path(DIR_OUT, "network")
DIR_EN  <- file.path(DIR_OUT, "hub_enrichment")
QCUT    <- 0.05
N_HUB_M <- 30

net <- readRDS(file.path(DIR_NET, "tcga_network.rds"))
hub <- as.data.table(net$hub); mod <- as.data.table(net$modules)
expressed <- readRDS(file.path(DIR_NET, "expressed_symbols.rds"))
E <- fread(file.path(DIR_NET, "tcga_deg_network_edges.csv.gz"))

ent <- function(s) {
  m <- as.data.table(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(s), keytype = "SYMBOL", columns = "ENTREZID")))
  m[!is.na(ENTREZID)]
}
uni_e <- unique(ent(expressed)$ENTREZID)

intra_degree <- function(m) {
  members <- mod[module == m]$symbol
  sub <- E[gene1 %in% members & gene2 %in% members]
  g <- graph_from_data_frame(sub[, .(gene1, gene2)], directed = FALSE,
                             vertices = data.frame(name = members))
  d <- degree(g)
  data.table(symbol = names(d), intra_degree = as.integer(d), module = m)
}

msig <- function(coll, sub = NULL) {
  d <- tryCatch(if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll)
                else msigdbr(species = "Homo sapiens", collection = coll,
                             subcollection = sub), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  ec <- intersect(c("ncbi_gene", "entrez_gene"), names(d))[1]
  x <- unique(data.frame(term = d$gs_name, gene = as.character(d[[ec]])))
  x[!is.na(x$gene), ]
}
MSIG <- list(MSigDB_Hallmark = list("H", NULL), MSigDB_C2_PID = list("C2", "CP:PID"),
             MSigDB_C2_BioCarta = list("C2", "CP:BIOCARTA"),
             MSigDB_C2_WikiPathways = list("C2", "CP:WIKIPATHWAYS"),
             MSigDB_C6_Oncogenic = list("C6", NULL), MSigDB_C8_CellType = list("C8", NULL),
             MSigDB_C2_CGP = list("C2", "CGP"))

enrich_set <- function(genes_sym, label) {
  g <- unique(ent(genes_sym)$ENTREZID)
  out <- list()
  add <- function(nm, x) {
    d <- tryCatch(as.data.table(x), error = function(e) data.table())
    if (nrow(d)) { d[, `:=`(database = nm, gene_set = label)]; out[[nm]] <<- d }
    msg(sprintf("   %-24s %s", nm, if (nrow(d)) paste0(nrow(d), " terms") else "none"))
  }
  for (o in c("BP", "MF", "CC"))
    add(paste0("GO_", o), setReadable(enrichGO(g, org.Hs.eg.db, keyType = "ENTREZID",
        ont = o, universe = uni_e, pAdjustMethod = "BH", qvalueCutoff = QCUT),
        org.Hs.eg.db, "ENTREZID"))
  add("KEGG", setReadable(enrichKEGG(g, "hsa", universe = uni_e,
      pAdjustMethod = "BH", qvalueCutoff = QCUT), org.Hs.eg.db, "ENTREZID"))
  add("Reactome", enrichPathway(g, universe = uni_e, pAdjustMethod = "BH",
      qvalueCutoff = QCUT, readable = TRUE))
  add("WikiPathways", setReadable(enrichWP(g, organism = "Homo sapiens",
      universe = uni_e, pAdjustMethod = "BH", qvalueCutoff = QCUT),
      org.Hs.eg.db, "ENTREZID"))
  add("DisGeNET", setReadable(enrichDGN(g, universe = uni_e, pAdjustMethod = "BH",
      qvalueCutoff = QCUT), org.Hs.eg.db, "ENTREZID"))
  for (nm in names(MSIG)) {
    t2g <- msig(MSIG[[nm]][[1]], MSIG[[nm]][[2]]); if (is.null(t2g)) next
    add(nm, setReadable(enricher(g, TERM2GENE = t2g, universe = uni_e,
        pAdjustMethod = "BH", qvalueCutoff = QCUT), org.Hs.eg.db, "ENTREZID"))
  }
  rbindlist(out, fill = TRUE)
}

big <- mod[, .N, by = module][N >= 100][order(-N)]$module
res <- list(); hubs <- list()
for (m in big) {
  idg <- intra_degree(m)
  setorder(idg, -intra_degree)
  idg[, module_rank := seq_len(.N)]
  hubs[[as.character(m)]] <- idg
  msg("\n=== module ", m, " (", nrow(idg), " genes) — enrichment ===")
  res[[as.character(m)]] <- enrich_set(idg$symbol, paste0("module_", m))
}
H <- rbindlist(hubs, fill = TRUE)
H <- merge(H, hub[, .(symbol, degree, logFC, padj)], by = "symbol", all.x = TRUE)
setorder(H, module, -intra_degree)
fwrite(H, file.path(DIR_NET, "module_intra_degree.csv"))

R <- rbindlist(res, fill = TRUE)
fwrite(R, file.path(DIR_EN, "module_enrichment_all_databases.csv"))

msg("\n\n########  MODULE 1 — the EMT / stromal-invasive programme  ########")
msg("\n--- top ", N_HUB_M, " module-1 hub genes (by intra-module degree) ---")
m1 <- H[module == 1][order(-intra_degree)]
print(as.data.frame(head(m1[, .(module_rank, symbol, intra_degree,
      whole_network_degree = degree, logFC, padj = signif(padj, 3))], N_HUB_M)))
fwrite(head(m1, N_HUB_M), file.path(DIR_EN, "module1_EMT_hub_genes.csv"))

msg("\n--- module 1: motility-relevant enriched terms, every database ---")
MOT <- "mesenchym|migrat|motil|invas|matrix|adhes|wound|focal|integrin|collagen|angiogen|tgf|metasta|actin|locomot"
d1 <- R[gene_set == "module_1" & tolower(Description) %like% MOT][order(p.adjust)]
print(as.data.frame(head(d1[, .(database, term = substr(Description, 1, 54),
      padj = signif(p.adjust, 3), Count)], 35)))

msg("\n--- module 1: top term per database ---")
print(as.data.frame(R[gene_set == "module_1"][order(p.adjust)][
  , .SD[1], by = database][, .(database, term = substr(Description, 1, 50),
  padj = signif(p.adjust, 3), Count)]))

msg("\n--- module 2 (contrast): top term per database ---")
print(as.data.frame(R[gene_set == "module_2"][order(p.adjust)][
  , .SD[1], by = database][, .(database, term = substr(Description, 1, 50),
  padj = signif(p.adjust, 3), Count)]))
