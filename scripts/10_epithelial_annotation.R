# Cell-lineage annotation for the nine shortlisted candidates.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
source(file.path(getwd(), "R", "utils.R"))

GENES <- c("GGCT","DLG3","CCT6A","PPP4C","TMED3","TUBA1C","DCTPP1","TOMM34","DYNLT1")
PUBMED <- c(GGCT = 16, DLG3 = 10, CCT6A = 9, PPP4C = 5, TMED3 = 5,
            TUBA1C = 14, DCTPP1 = 15, TOMM34 = 6, DYNLT1 = 5)

MARKERS <- list(
  epithelial  = c("EPCAM","CDH1","KRT8","KRT18","KRT19","KRT5","KRT14","ELF3"),
  immune      = c("PTPRC","CD3E","CD68","CD19","LYZ"),
  fibroblast  = c("COL1A1","COL1A2","THY1","FAP","PDGFRB","VIM"),
  endothelial = c("PECAM1","VWF","CDH5","ENG"),
  adipocyte   = c("ADIPOQ","LEP","PLIN1","FABP4")
)

map <- as.data.table(AnnotationDbi::select(
  org.Hs.eg.db, keys = GENES, keytype = "SYMBOL",
  columns = c("ENSEMBL", "GENENAME")))
map <- map[!duplicated(SYMBOL)]

hpa_rows <- list()
for (i in seq_len(nrow(map))) {
  sym <- map$SYMBOL[i]; ens <- map$ENSEMBL[i]
  url <- sprintf("https://www.proteinatlas.org/%s.json", ens)
  j <- tryCatch(fromJSON(url), error = function(e) NULL)
  if (is.null(j)) { warning("HPA fetch failed: ", sym); next }

  flat <- function(x) if (is.null(x) || !length(x)) NA_character_
                      else paste(unlist(x), collapse = "; ")
  named <- function(x) {
    if (is.null(x) || !length(x)) return(NA_character_)
    paste(sprintf("%s (%s)", names(unlist(x)), unlist(x)), collapse = "; ")
  }
  prog <- j[["Cancer prognostics - Breast Invasive Carcinoma (TCGA)"]]

  hpa_rows[[sym]] <- data.table(
    gene              = sym,
    description       = map$GENENAME[i],
    sc_specificity    = flat(j[["RNA single cell type specificity"]]),
    sc_distribution   = flat(j[["RNA single cell type distribution"]]),
    sc_enhanced_types = named(j[["RNA single cell type specific nCPM"]]),
    sc_group_types    = named(j[["RNA single cell type group specific nCPM"]]),
    tissue_breast_nTPM = {
      u2 <- sprintf(paste0("https://www.proteinatlas.org/api/search_download.php",
                           "?search=%s&format=json&columns=g,t_RNA_breast&compress=no"), sym)
      v <- tryCatch(fromJSON(u2), error = function(e) NULL)
      if (is.null(v) || !nrow(v)) NA_character_
      else as.character(v[[grep("breast", names(v), value = TRUE)[1]]][1])
    },
    subcellular       = flat(j[["Subcellular main location"]]),
    protein_class     = flat(j[["Protein class"]]),
    molecular_function = flat(j[["Molecular function"]]),
    hpa_brca_prognostic = if (is.null(prog)) NA_character_
      else sprintf("%s (p=%s)",
                   if (isTRUE(prog$is_prognostic))
                     paste(prog$prognostic, "-", prog$`prognostic type`)
                   else "unprognostic",
                   prog$p_val)
  )
  msg("HPA: ", sym)
}
hpa <- rbindlist(hpa_rows, fill = TRUE)
fwrite(hpa, file.path(DIR_OUT, "candidates_HPA_annotation.csv"))

counts <- readRDS(file.path(PROJ_ROOT, "exports", "r_objects", "rna.rds"))
counts <- counts[, substr(colnames(counts), 14, 15) == "01", drop = FALSE]
pat    <- substr(colnames(counts), 1, 12)
counts <- counts[, !duplicated(pat), drop = FALSE]

dge  <- DGEList(counts = counts)
dge  <- calcNormFactors(dge, method = "TMM")
lcpm <- cpm(dge, log = TRUE, prior.count = 1)
rownames(lcpm) <- strip_ens_version(rownames(lcpm))
lcpm <- lcpm[rowMeans(lcpm) > 0, , drop = FALSE]
msg("TCGA tumours: ", ncol(lcpm), " patients, ", nrow(lcpm), " expressed genes")

sym2ens <- function(s) {
  m <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = s, keytype = "SYMBOL", columns = "ENSEMBL"))
  m <- m[!is.na(m$ENSEMBL) & m$ENSEMBL %in% rownames(lcpm), ]
  m[!duplicated(m$SYMBOL), ]
}

zrow <- function(m) t(scale(t(m)))
scores <- sapply(names(MARKERS), function(k) {
  mm <- sym2ens(MARKERS[[k]])
  colMeans(zrow(lcpm[mm$ENSEMBL, , drop = FALSE]), na.rm = TRUE)
})
msg("compartment scores built from ",
    paste(sprintf("%s=%d", names(MARKERS), sapply(MARKERS, length)), collapse = ", "))

gm <- sym2ens(GENES)
cor_tab <- rbindlist(lapply(seq_len(nrow(gm)), function(i) {
  e <- lcpm[gm$ENSEMBL[i], ]
  as.data.table(c(list(gene = gm$SYMBOL[i]),
                  setNames(as.list(round(cor(e, scores, method = "spearman"), 3)),
                           colnames(scores))))
}))
cor_tab[, top_compartment := names(MARKERS)[max.col(as.matrix(.SD))],
        .SDcols = names(MARKERS)]
cor_tab[, epithelial_margin := round(epithelial -
          pmax(immune, fibroblast, endothelial, adipocyte), 3)]
setorder(cor_tab, -epithelial)

fwrite(cor_tab, file.path(DIR_OUT, "candidates_lineage_correlation.csv"))

out <- merge(cor_tab, hpa[, .(gene, sc_specificity, sc_group_types,
                              tissue_breast_nTPM, hpa_brca_prognostic)],
             by = "gene", all.x = TRUE)
out[, pubmed_hits := PUBMED[gene]]
setorder(out, -epithelial)
fwrite(out, file.path(DIR_OUT, "candidates_epithelial_summary.csv"))

msg("\n=== lineage correlation (Spearman vs compartment score) ===")
print(as.data.frame(cor_tab))
msg("\n=== HPA ===")
print(as.data.frame(hpa[, .(gene, sc_specificity, tissue_breast_nTPM, hpa_brca_prognostic)]))
