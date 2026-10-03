# Hub modules: proliferation vs EMT scores

suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(msigdbr)
  library(org.Hs.eg.db); library(AnnotationDbi)
})
source(file.path(getwd(), "R", "utils.R"))

DIR_NET <- file.path(DIR_OUT, "network")
DIR_EN  <- file.path(DIR_OUT, "hub_enrichment")

net <- readRDS(file.path(DIR_NET, "tcga_network.rds"))
hub <- as.data.table(net$hub)
hub_sym <- hub[is_hub == TRUE]$symbol
expressed <- readRDS(file.path(DIR_NET, "expressed_symbols.rds"))

ent <- function(s) {
  m <- as.data.table(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(s), keytype = "SYMBOL", columns = "ENTREZID")))
  m[!is.na(ENTREZID)]
}
hub_e <- unique(ent(hub_sym)$ENTREZID)
uni_e <- unique(ent(expressed)$ENTREZID)

msig_all <- rbindlist(lapply(
  list(c("H", NA), c("C2", "CP:REACTOME"), c("C2", "CP:WIKIPATHWAYS"),
       c("C2", "CP:BIOCARTA"), c("C2", "CP:PID"), c("C5", "GO:BP")),
  function(k) {
    d <- tryCatch(if (is.na(k[2])) msigdbr(species = "Homo sapiens", collection = k[1])
                  else msigdbr(species = "Homo sapiens", collection = k[1],
                               subcollection = k[2]), error = function(e) NULL)
    if (is.null(d)) return(data.table())
    ec <- intersect(c("ncbi_gene", "entrez_gene"), names(d))[1]
    data.table(term = d$gs_name, gene = as.character(d[[ec]]))
  }), fill = TRUE)
msig_all <- unique(msig_all[!is.na(gene)])

TARGETS <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_ANGIOGENESIS", "HALLMARK_APICAL_JUNCTION",
  "GOBP_CELL_MIGRATION", "GOBP_CELL_MOTILITY", "GOBP_AMEBOIDAL_TYPE_CELL_MIGRATION",
  "GOBP_REGULATION_OF_EPITHELIAL_CELL_MIGRATION",
  "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION",
  "GOBP_WOUND_HEALING", "GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION", "GOBP_CELL_CELL_ADHESION",
  "GOBP_CELL_SUBSTRATE_ADHESION", "GOBP_ACTIN_FILAMENT_ORGANIZATION",
  "GOBP_MITOTIC_NUCLEAR_DIVISION", "GOBP_DNA_REPLICATION",
  "REACTOME_CELL_CYCLE", "REACTOME_CELL_CYCLE_MITOTIC",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
  "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX",
  "REACTOME_INTEGRIN_SIGNALING", "REACTOME_SMOOTH_MUSCLE_CONTRACTION",
  "WP_FOCAL_ADHESION", "WP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION_IN_COLORECTAL_CANCER",
  "WP_CELL_CYCLE", "PID_INTEGRIN1_PATHWAY", "BIOCARTA_CELLCYCLE_PATHWAY"
)

N <- length(uni_e); K <- length(intersect(hub_e, uni_e))
rows <- list()
for (tm in TARGETS) {
  sg <- unique(msig_all[term == tm]$gene)
  sg <- intersect(sg, uni_e)
  if (!length(sg)) { rows[[tm]] <- data.table(term = tm, set_size = 0L,
      overlap = NA_integer_, expected = NA_real_, fold = NA_real_,
      pvalue = NA_real_, hits = NA_character_); next }
  ov <- intersect(hub_e, sg)
  p  <- phyper(length(ov) - 1L, length(sg), N - length(sg), K, lower.tail = FALSE)
  exp_n <- K * length(sg) / N
  sy <- unique(ent(hub_sym)[ENTREZID %in% ov]$SYMBOL)
  rows[[tm]] <- data.table(term = tm, set_size = length(sg), overlap = length(ov),
    expected = round(exp_n, 2), fold = round(length(ov) / exp_n, 2),
    pvalue = p, hits = paste(sort(sy), collapse = ";"))
}
dt <- rbindlist(rows)
dt[, padj := p.adjust(pvalue, "BH")]
dt[, theme := fifelse(grepl("MITOTIC|CELL_?CYCLE|E2F|G2M|MYC|DNA_REPL", term), "proliferation",
             fifelse(grepl("MESENCHYMAL|MIGRATION|MOTILITY|ADHESION|MATRIX|WOUND|FOCAL|TGF|ANGIO|INTEGRIN|ACTIN|JUNCTION|MUSCLE", term), "EMT / motility", "other"))]
setorder(dt, theme, pvalue)
fwrite(dt, file.path(DIR_EN, "directed_emt_vs_proliferation.csv"))

msg("\n=== DIRECTED TEST: ", K, " hub genes vs ", N, " expressed genes ===")
for (th in c("proliferation", "EMT / motility")) {
  cat("\n---- ", toupper(th), " ----\n", sep = "")
  d <- dt[theme == th]
  print(as.data.frame(d[, .(term = substr(term, 1, 52), set_size, overlap,
                            expected, fold, padj = signif(padj, 3))]))
}

cgp <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
ec  <- intersect(c("ncbi_gene", "entrez_gene"), names(cgp))[1]
cgp <- data.table(term = cgp$gs_name, gene = as.character(cgp[[ec]]))

prolif_ref <- unique(msig_all[term %in% c("HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS", "GOBP_MITOTIC_NUCLEAR_DIVISION",
  "GOBP_DNA_REPLICATION", "GOBP_CHROMOSOME_SEGREGATION")]$gene)
emt_ref <- unique(msig_all[term %in% c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "GOBP_CELL_MIGRATION", "GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION")]$gene)

res20 <- fread(file.path(DIR_EN, "hub_enrichment_all_databases.csv"))
susp <- unique(res20[database == "MSigDB_C2_CGP" & background == "expressed" &
       tolower(Description) %like% "metasta|invasi|mesenchym"]$Description)

fr <- list()
for (s in susp) {
  ov <- intersect(hub_e, unique(cgp[term == s]$gene))
  if (!length(ov)) next
  fr[[s]] <- data.table(signature = s, overlap = length(ov),
    n_proliferation = length(intersect(ov, prolif_ref)),
    n_emt_motility  = length(intersect(ov, emt_ref)),
    n_neither = length(setdiff(ov, union(prolif_ref, emt_ref))))
}
fo <- rbindlist(fr)
if (nrow(fo)) {
  fo[, pct_proliferation := round(100 * n_proliferation / overlap)]
  fo[, pct_emt := round(100 * n_emt_motility / overlap)]
  setorder(fo, -overlap)
  fwrite(fo, file.path(DIR_EN, "metastasis_signature_forensics.csv"))
  msg("\n=== What is actually driving the 'metastasis/invasive' signature hits? ===")
  print(as.data.frame(fo[, .(signature = substr(signature, 1, 46), overlap,
        n_proliferation, n_emt_motility, pct_proliferation, pct_emt)]))
}

EMT_GENES <- c("VIM","CDH2","CDH1","FN1","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2",
               "MMP2","MMP9","MMP11","MMP14","INHBA","COL10A1","COL1A1","SPARC",
               "TGFB1","ACTA2","POSTN","SERPINE1","ITGB1","CTNNB1")
em <- hub[symbol %in% EMT_GENES][order(rank)]
msg("\n=== canonical EMT/invasion genes in the TCGA DEG network ===")
if (nrow(em)) print(as.data.frame(em[, .(symbol, degree, rank, module, logFC,
                                          is_hub)])) else cat("  none present\n")
cat("\n  hub degree: median ", median(hub[is_hub == TRUE]$degree),
    " | EMT-gene degree: median ", if (nrow(em)) median(em$degree) else NA, "\n", sep = "")
fwrite(em, file.path(DIR_EN, "emt_gene_network_position.csv"))

mod <- as.data.table(net$modules)
big <- mod[, .N, by = module][N >= 30][order(-N)]
msg("\n=== per-module Hallmark enrichment (modules with >= 30 genes) ===")
hm  <- msig_all[term %like% "^HALLMARK_"]
mrows <- list()
for (m in big$module) {
  g <- unique(ent(mod[module == m]$symbol)$ENTREZID)
  e <- tryCatch(as.data.table(enricher(g, TERM2GENE = hm, universe = uni_e,
        pAdjustMethod = "BH", qvalueCutoff = 0.05)), error = function(e) data.table())
  top <- if (nrow(e)) paste(sub("^HALLMARK_", "", head(e[order(p.adjust)]$Description, 3)),
                            collapse = "; ") else "none"
  emt_p <- if (nrow(e) && "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" %in% e$Description)
             signif(e[Description == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]$p.adjust, 3) else NA
  mrows[[as.character(m)]] <- data.table(module = m, n_genes = big[module == m]$N,
      top_hallmarks = top, EMT_padj = emt_p)
}
mo <- rbindlist(mrows)
setorder(mo, -n_genes)
fwrite(mo, file.path(DIR_EN, "module_hallmark_summary.csv"))
print(as.data.frame(mo))
