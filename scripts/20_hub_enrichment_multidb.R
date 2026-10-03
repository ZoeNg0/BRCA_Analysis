# Over-representation analysis of hub genes

suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(ReactomePA)
  library(DOSE); library(msigdbr); library(org.Hs.eg.db); library(AnnotationDbi)
})
source(file.path(getwd(), "R", "utils.R"))

DIR_NET <- file.path(DIR_OUT, "network")
DIR_EN  <- file.path(DIR_OUT, "hub_enrichment")
dir.create(DIR_EN, recursive = TRUE, showWarnings = FALSE)
QCUT <- 0.05

net  <- readRDS(file.path(DIR_NET, "tcga_network.rds"))
hub  <- as.data.table(net$hub)
msg("hub genes: ", sum(hub$is_hub), " of ", nrow(hub), " network nodes")

expressed <- readRDS(file.path(DIR_NET, "expressed_symbols.rds"))
sym2ent <- function(s) {
  m <- as.data.table(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(s), keytype = "SYMBOL", columns = "ENTREZID")))
  unique(m[!is.na(ENTREZID)]$ENTREZID)
}
hub_sym <- hub[is_hub == TRUE]$symbol
hub_ent <- sym2ent(hub_sym)
bg <- list(expressed = sym2ent(expressed), deg = sym2ent(hub$symbol))
msg("hub ENTREZ: ", length(hub_ent),
    " | background expressed: ", length(bg$expressed),
    " | background DEG: ", length(bg$deg))

msig_get <- function(coll, sub = NULL) {
  d <- tryCatch({
    if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll)
    else msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
  }, error = function(e) NULL)
  if (is.null(d)) return(NULL)
  ec <- intersect(c("ncbi_gene", "entrez_gene"), names(d))[1]
  x <- unique(data.frame(term = d$gs_name, gene = as.character(d[[ec]])))
  x[!is.na(x$gene), ]
}

MSIG <- list(
  MSigDB_Hallmark        = list("H",  NULL),
  MSigDB_C2_CGP          = list("C2", "CGP"),
  MSigDB_C2_BioCarta     = list("C2", "CP:BIOCARTA"),
  MSigDB_C2_PID          = list("C2", "CP:PID"),
  MSigDB_C2_WikiPathways = list("C2", "CP:WIKIPATHWAYS"),
  MSigDB_C2_KEGG_MEDICUS = list("C2", "CP:KEGG_MEDICUS"),
  MSigDB_C3_TFT_GTRD     = list("C3", "TFT:GTRD"),
  MSigDB_C3_MIR_MIRDB    = list("C3", "MIR:MIRDB"),
  MSigDB_C4_CancerModule = list("C4", "CM"),
  MSigDB_C5_HPO          = list("C5", "HPO"),
  MSigDB_C6_Oncogenic    = list("C6", NULL),
  MSigDB_C8_CellType     = list("C8", NULL)
)

tidy <- function(x, db, bgname) {
  if (is.null(x)) return(data.table())
  d <- as.data.table(x)
  if (!nrow(d)) return(data.table())
  d[, `:=`(database = db, background = bgname)]
  d[]
}
pretty_ms <- function(x)
  tools::toTitleCase(tolower(gsub("_", " ",
    sub("^(HALLMARK|BIOCARTA|PID|WP|REACTOME|KEGG_MEDICUS|GOBP|HP)_", "", x))))

run_all <- function(genes, universe, bgname) {
  res <- list()
  add <- function(nm, expr) {
    r <- tryCatch(expr, error = function(e) { msg("   ! ", nm, ": ", conditionMessage(e)); NULL })
    d <- tidy(r, nm, bgname)
    msg(sprintf("   %-24s %s", nm, if (nrow(d)) paste0(nrow(d), " terms") else "none"))
    if (nrow(d)) res[[nm]] <<- d
  }
  for (ont in c("BP", "MF", "CC"))
    add(paste0("GO_", ont), setReadable(enrichGO(genes, OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID", ont = ont, universe = universe,
        pAdjustMethod = "BH", qvalueCutoff = QCUT), org.Hs.eg.db, "ENTREZID"))
  add("KEGG", setReadable(enrichKEGG(genes, organism = "hsa", universe = universe,
      pAdjustMethod = "BH", qvalueCutoff = QCUT), org.Hs.eg.db, "ENTREZID"))
  add("KEGG_Module", setReadable(enrichMKEGG(genes, organism = "hsa",
      universe = universe, pAdjustMethod = "BH", qvalueCutoff = QCUT),
      org.Hs.eg.db, "ENTREZID"))
  add("Reactome", enrichPathway(genes, universe = universe, pAdjustMethod = "BH",
      qvalueCutoff = QCUT, readable = TRUE))
  add("WikiPathways", setReadable(enrichWP(genes, organism = "Homo sapiens",
      universe = universe, pAdjustMethod = "BH", qvalueCutoff = QCUT),
      org.Hs.eg.db, "ENTREZID"))
  add("DiseaseOntology", setReadable(enrichDO(genes, universe = universe,
      pAdjustMethod = "BH", qvalueCutoff = QCUT), org.Hs.eg.db, "ENTREZID"))
  add("DisGeNET", setReadable(enrichDGN(genes, universe = universe,
      pAdjustMethod = "BH", qvalueCutoff = QCUT), org.Hs.eg.db, "ENTREZID"))
  add("NCG_CancerGenes", setReadable(enrichNCG(genes, universe = universe,
      pAdjustMethod = "BH", qvalueCutoff = QCUT), org.Hs.eg.db, "ENTREZID"))
  for (nm in names(MSIG)) {
    t2g <- msig_get(MSIG[[nm]][[1]], MSIG[[nm]][[2]])
    if (is.null(t2g)) { msg("   ! ", nm, ": unavailable"); next }
    add(nm, setReadable(enricher(genes, TERM2GENE = t2g, universe = universe,
        pAdjustMethod = "BH", qvalueCutoff = QCUT), org.Hs.eg.db, "ENTREZID"))
  }
  rbindlist(res, fill = TRUE)
}

all_res <- list()
for (bn in names(bg)) {
  msg("\n=== background: ", bn, " (", length(bg[[bn]]), " genes) ===")
  all_res[[bn]] <- run_all(hub_ent, bg[[bn]], bn)
}
res <- rbindlist(all_res, fill = TRUE)
res[, Description_pretty := pretty_ms(Description)]
fwrite(res, file.path(DIR_EN, "hub_enrichment_all_databases.csv"))
msg("\ntotal enriched terms: ", nrow(res), " across ",
    uniqueN(res$database), " databases")

PROLIF <- "prolifer|cell cycle|mitotic|mitosis|dna replication|g2m|e2f|chromosome segregation|spindle|kinetochore|nuclear division|cytokinesis|dna repair|myc target"
EMT    <- "epithelial.mesenchymal|\\bemt\\b|cell migration|cell motility|locomotion|invasion|invasive|metastasis|metastatic|extracellular matrix|focal adhesion|cell adhesion|wound healing|actin cytoskeleton|tgf.?beta|wnt|cadherin|integrin|matrix metallo|angiogenesis"

flag <- res[, .(database, background, ID, term = Description_pretty,
                padj = p.adjust, qvalue, GeneRatio, BgRatio, Count, geneID)]
flag[, theme := fifelse(grepl(PROLIF, tolower(term)), "proliferation",
              fifelse(grepl(EMT, tolower(term)), "EMT / migration / invasion", "other"))]
fwrite(flag, file.path(DIR_EN, "hub_terms_themed.csv"))

msg("\n=== themed term counts (FDR < ", QCUT, ") ===")
print(dcast(flag[, .N, by = .(background, theme)], background ~ theme,
            value.var = "N", fill = 0))

for (bn in names(bg)) {
  msg("\n--- background ", bn, ": EMT / migration / invasion terms ---")
  d <- flag[background == bn & theme == "EMT / migration / invasion"][order(padj)]
  if (!nrow(d)) cat("   NONE enriched at FDR < ", QCUT, "\n", sep = "")
  else print(as.data.frame(head(d[, .(database, term, padj, Count)], 25)))
}

msg("\n--- top proliferation terms (expressed background) ---")
print(as.data.frame(head(flag[background == "expressed" & theme == "proliferation"][
  order(padj), .(database, term, padj, Count)], 20)))
