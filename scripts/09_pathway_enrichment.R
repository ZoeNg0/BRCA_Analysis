# Pathway enrichment for each shortlisted candidate.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})
source(file.path(getwd(), "R", "utils.R"))

GENES  <- c("GGCT","DLG3","CCT6A","PPP4C","TMED3","TUBA1C","DCTPP1","TOMM34","DYNLT1")
TOP_N  <- 300
DIR_EN <- file.path(DIR_OUT, "enrichment")
DIR_EP <- file.path(DIR_PLOTS, "enrichment")
for (d in c(DIR_EN, DIR_EP)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(PROJ_ROOT, "exports", "r_objects", "rna.rds"))
counts <- counts[, substr(colnames(counts), 14, 15) == "01", drop = FALSE]
counts <- counts[, !duplicated(substr(colnames(counts), 1, 12)), drop = FALSE]
dge  <- calcNormFactors(DGEList(counts = counts), method = "TMM")
lcpm <- cpm(dge, log = TRUE, prior.count = 1)
rownames(lcpm) <- strip_ens_version(rownames(lcpm))
lcpm <- lcpm[rowMeans(lcpm) > 0, , drop = FALSE]
msg("TCGA tumours: ", ncol(lcpm), " patients x ", nrow(lcpm), " expressed genes")

uni <- as.data.table(suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = rownames(lcpm), keytype = "ENSEMBL", columns = "ENTREZID")))
uni <- uni[!is.na(ENTREZID)][!duplicated(ENSEMBL)]
universe <- unique(uni$ENTREZID)
msg("enrichment background: ", length(universe), " genes with ENTREZ ids")

gm <- as.data.table(suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = GENES, keytype = "SYMBOL", columns = "ENSEMBL")))
gm <- gm[!is.na(ENSEMBL) & ENSEMBL %in% rownames(lcpm)][!duplicated(SYMBOL)]

msg("pre-ranking expression matrix for Spearman ...")
rank_t <- t(t(apply(lcpm, 1, rank)))
msg("ranked matrix: ", nrow(rank_t), " samples x ", ncol(rank_t), " genes")

hm <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
               error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
ent_col <- intersect(c("ncbi_gene", "entrez_gene"), names(hm))[1]
hm_t2g  <- unique(data.frame(term = hm$gs_name, gene = as.character(hm[[ent_col]])))
hm_t2g  <- hm_t2g[!is.na(hm_t2g$gene), ]
msg("Hallmark sets: ", length(unique(hm_t2g$term)))

pretty_hm <- function(x) tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", x))))

top_terms <- list(); hm_all <- list(); membership <- list()

for (i in seq_len(nrow(gm))) {
  g <- gm$SYMBOL[i]; eid <- gm$ENSEMBL[i]
  msg("=== ", g, " ===")

  r  <- cor(rank_t[, eid], rank_t)[1, ]
  r  <- sort(r[names(r) != eid], decreasing = TRUE)
  nb <- names(head(r, TOP_N))
  msg("  neighbourhood rho: ", round(r[nb[1]], 3), " to ", round(r[nb[TOP_N]], 3))
  nb_entrez <- unique(uni[ENSEMBL %in% nb, ENTREZID])

  go <- enrichGO(nb_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                 universe = universe, pAdjustMethod = "BH",
                 qvalueCutoff = 0.05, readable = TRUE)
  kg <- tryCatch(setReadable(enrichKEGG(nb_entrez, organism = "hsa",
                                        universe = universe, pAdjustMethod = "BH",
                                        qvalueCutoff = 0.05),
                             org.Hs.eg.db, "ENTREZID"),
                 error = function(e) NULL)
  re <- tryCatch(enrichPathway(nb_entrez, universe = universe, pAdjustMethod = "BH",
                               qvalueCutoff = 0.05, readable = TRUE),
                 error = function(e) NULL)
  hl <- tryCatch(setReadable(enricher(nb_entrez, TERM2GENE = hm_t2g,
                                      universe = universe, pAdjustMethod = "BH",
                                      qvalueCutoff = 0.05),
                             org.Hs.eg.db, "ENTREZID"),
                 error = function(e) NULL)

  dt <- function(x) if (!is.null(x)) as.data.table(x) else data.table()
  gd <- dt(go); kd <- dt(kg); rd <- dt(re); hd <- dt(hl)
  for (nm in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
    d <- switch(nm, GO_BP = gd, KEGG = kd, Reactome = rd, Hallmark = hd)
    if (nrow(d)) fwrite(d, file.path(DIR_EN, sprintf("%s_%s.csv", g, nm)))
  }
  msg("  GO:BP ", nrow(gd), " | KEGG ", nrow(kd),
      " | Reactome ", nrow(rd), " | Hallmark ", nrow(hd))

  for (nm in c("GO_BP", "KEGG", "Reactome", "Hallmark")) {
    d <- switch(nm, GO_BP = gd, KEGG = kd, Reactome = rd, Hallmark = hd)
    if (nrow(d))
      top_terms[[paste(g, nm)]] <- data.table(
        gene = g, source = nm,
        term = if (nm == "Hallmark") pretty_hm(d$Description) else d$Description,
        padj = d$p.adjust, count = d$Count)[1:min(10, nrow(d))]
  }
  if (nrow(hd))
    hm_all[[g]] <- data.table(gene = g, term = pretty_hm(hd$Description),
                              padj = hd$p.adjust, count = hd$Count)

  ge <- uni[ENSEMBL == eid, ENTREZID][1]
  if (!is.na(ge)) {
    mrows <- list()
    hmm <- unique(hm_t2g$term[hm_t2g$gene == ge])
    if (length(hmm)) mrows[[1]] <- data.table(gene = g, db = "Hallmark",
                                              pathway = pretty_hm(hmm))
    kt <- tryCatch({
      kk <- clusterProfiler::download_KEGG("hsa")
      ids <- kk$KEGGPATHID2EXTID$to[kk$KEGGPATHID2EXTID$from %in%
               kk$KEGGPATHID2EXTID$from[kk$KEGGPATHID2EXTID$to == ge]]
      pid <- unique(kk$KEGGPATHID2EXTID$from[kk$KEGGPATHID2EXTID$to == ge])
      unique(kk$KEGGPATHID2NAME$to[kk$KEGGPATHID2NAME$from %in% pid])
    }, error = function(e) character(0))
    if (length(kt)) mrows[[2]] <- data.table(gene = g, db = "KEGG", pathway = kt)
    if (length(mrows)) membership[[g]] <- rbindlist(mrows, fill = TRUE)
  }

  if (nrow(hd)) {
    ph <- ggplot(head(hd[order(p.adjust)], 10),
                 aes(x = Count, y = reorder(pretty_hm(Description), Count))) +
      geom_col(fill = "#2a78d6", width = 0.62) +
      labs(title = paste0(g, " — Hallmarks of cancer"),
           subtitle = sprintf("Enrichment of the %d genes most co-expressed with %s across %d TCGA-BRCA tumours",
                              TOP_N, g, ncol(lcpm)),
           x = "Genes in the neighbourhood", y = NULL) +
      theme_minimal(base_size = 10) +
      theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA),
            panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(colour = "#e6e5e1", linewidth = 0.25),
            plot.title    = element_text(face = "bold", size = 12, colour = "#0b0b0b"),
            plot.subtitle = element_text(size = 8, colour = "#52514e", margin = margin(b = 6)),
            axis.text.y   = element_text(size = 9, colour = "#0b0b0b"),
            axis.text.x   = element_text(colour = "#52514e"))
    ggsave(file.path(DIR_EP, paste0("hallmark_", g, ".png")), ph,
           width = 7.6, height = 4.4, dpi = 260, bg = "#fcfcfb")
  }
  if (nrow(gd)) {
    pg <- dotplot(go, showCategory = 12) +
      labs(title = paste0(g, " — GO biological process")) +
      theme_minimal(base_size = 10) +
      theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA),
            panel.grid.minor = element_blank(),
            plot.title  = element_text(face = "bold", size = 12, colour = "#0b0b0b"),
            axis.text.y = element_text(size = 8, colour = "#0b0b0b"),
            axis.text.x = element_text(colour = "#52514e"))
    ggsave(file.path(DIR_EP, paste0("GO_", g, ".png")), pg,
           width = 8.2, height = 6, dpi = 260, bg = "#fcfcfb")
  }
}

terms <- rbindlist(top_terms, fill = TRUE)
fwrite(terms, file.path(DIR_EN, "top_terms_per_gene.csv"))
if (length(membership))
  fwrite(rbindlist(membership, fill = TRUE),
         file.path(DIR_EN, "pathway_membership_per_gene.csv"))

hmm <- rbindlist(hm_all, fill = TRUE)
if (nrow(hmm)) {
  fwrite(hmm, file.path(DIR_EN, "hallmark_all_candidates.csv"))
  keep <- hmm[, .N, by = term][order(-N)][1:min(16, .N), term]
  hp <- hmm[term %in% keep]
  hp[, gene := factor(gene, levels = GENES)]
  hp[, term := factor(term, levels = rev(sort(unique(term))))]

  p <- ggplot(hp, aes(x = gene, y = term, fill = -log10(padj))) +
    geom_tile(colour = "#fcfcfb", linewidth = 1.4) +
    scale_fill_gradient(low = "#cde2fb", high = "#184f95",
                        name = expression(-log[10]~"FDR")) +
    labs(title = "Hallmarks of cancer across the nine candidates",
         subtitle = paste0("Enrichment of each gene's ", TOP_N,
                           "-gene co-expression neighbourhood in TCGA-BRCA tumours. ",
                           "Blank = not enriched (FDR > 0.05)."),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA),
          panel.grid = element_blank(),
          plot.title    = element_text(face = "bold", size = 13, colour = "#0b0b0b"),
          plot.subtitle = element_text(size = 8.5, colour = "#52514e", margin = margin(b = 8)),
          axis.text.x   = element_text(face = "italic", colour = "#0b0b0b", angle = 40, hjust = 1),
          axis.text.y   = element_text(colour = "#0b0b0b", size = 9),
          legend.position = "right")
  f <- file.path(DIR_PLOTS, "hallmark_candidate_matrix.png")
  ggsave(f, p, width = 9.4, height = 6.4, dpi = 300, bg = "#fcfcfb")
  msg("wrote ", f)
}

msg("\n=== Hallmark enrichment per candidate (top 5) ===")
for (g in GENES) {
  d <- hmm[gene == g][order(padj)]
  cat("\n", g, "  (", nrow(d), " hallmarks)\n", sep = "")
  if (nrow(d)) for (k in seq_len(min(5, nrow(d))))
    cat(sprintf("   %-46s FDR=%.2g  n=%d\n", substr(d$term[k],1,46), d$padj[k], d$count[k]))
}
