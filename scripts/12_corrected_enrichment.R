# Positionally corrected enrichment.

suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(msigdbr)
  library(org.Hs.eg.db); library(AnnotationDbi)
})
source(file.path(getwd(), "R", "utils.R"))

GENES  <- strsplit(Sys.getenv("GENES", "PPP4C,TOMM34"), ",")[[1]]
DIR_EN <- file.path(DIR_OUT, "enrichment")

ens_universe <- NULL

hm <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
               error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
ent_col <- intersect(c("ncbi_gene", "entrez_gene"), names(hm))[1]
hm_t2g  <- unique(data.frame(term = hm$gs_name, gene = as.character(hm[[ent_col]])))
pretty_hm <- function(x) tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", x))))

out <- list()
for (g in GENES) {
  nbf <- file.path(DIR_OUT, paste0(g, "_neighbourhoods.rds"))
  if (!file.exists(nbf)) { warning("missing ", nbf, " — run 11_target_dossier.R for ", g); next }
  nb <- readRDS(nbf)

  if (is.null(ens_universe)) {
    lc <- readRDS(file.path(PROJ_ROOT, "exports", "r_objects", "rna.rds"))
    ens_universe <- strip_ens_version(rownames(lc)); rm(lc)
  }
  u <- as.data.table(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = ens_universe, keytype = "ENSEMBL", columns = "ENTREZID")))
  u <- u[!is.na(ENTREZID)][!duplicated(ENSEMBL)]
  universe <- unique(u$ENTREZID)

  for (mode in c("raw", "offchrom")) {
    nbs <- if (mode == "raw") nb$nb_raw else nb$nb_off
    ent <- unique(u[ENSEMBL %in% nbs, ENTREZID])
    hl <- tryCatch(as.data.table(enricher(ent, TERM2GENE = hm_t2g, universe = universe,
                                          pAdjustMethod = "BH", qvalueCutoff = 0.05)),
                   error = function(e) data.table())
    if (nrow(hl))
      out[[paste(g, mode)]] <- data.table(gene = g, mode = mode,
        term = pretty_hm(hl$Description), padj = hl$p.adjust, count = hl$Count)
    msg(g, " [", mode, "] genes=", length(ent), "  hallmarks=", nrow(hl))
  }
}

res <- rbindlist(out, fill = TRUE)
fwrite(res, file.path(DIR_EN, "hallmark_positionally_corrected.csv"))

for (g in GENES) {
  msg("\n=== ", g, " ===")
  for (md in c("raw", "offchrom")) {
    d <- res[gene == g & mode == md][order(padj)]
    cat("  ", toupper(md), " (", nrow(d), " hallmarks)\n", sep = "")
    if (!nrow(d)) cat("     none\n")
    else for (k in seq_len(min(5, nrow(d))))
      cat(sprintf("     %-42s FDR=%.2g n=%d\n", substr(d$term[k],1,42), d$padj[k], d$count[k]))
  }
}
