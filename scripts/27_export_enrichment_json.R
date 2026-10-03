# Hub gene enrichment (JSON export)

suppressPackageStartupMessages({
  library(jsonlite)
})
source(file.path(getwd(), "R", "utils.R"))

DIR_EXP <- file.path(PROJ_ROOT, "exports")
DIR_EN  <- file.path(DIR_OUT, "hub_enrichment")
DIR_NET <- file.path(DIR_OUT, "network")
PADJ_CUT <- 0.05

SLOT_PROLIF <- 1L
SLOT_EMT    <- 2L

THEME_LABEL <- c(proliferation = "Proliferation",
                 `EMT / motility` = "EMT / motility")

f_in <- file.path(DIR_EN, "directed_emt_vs_proliferation.csv")
if (!file.exists(f_in))
  stop("missing ", f_in, " — run scripts/21_emt_vs_proliferation.R first")
d <- read.csv(f_in, stringsAsFactors = FALSE)
msg("read ", basename(f_in), ": ", nrow(d), " gene sets")

unresolved <- d[d$set_size == 0, "term"]
if (length(unresolved))
  msg("dropping ", length(unresolved), " unresolved set name(s): ",
      paste(unresolved, collapse = ", "))
d <- d[d$set_size > 0, , drop = FALSE]
d <- d[d$theme %in% names(THEME_LABEL), , drop = FALSE]
stopifnot(nrow(d) > 0)

src_of <- function(x) {
  p <- sub("_.*$", "", x)
  ifelse(p == "GOBP", "GO:BP",
  ifelse(p == "HALLMARK", "Hallmark",
  ifelse(p == "REACTOME", "Reactome",
  ifelse(p == "WP", "WikiPathways",
  ifelse(p == "PID", "PID",
  ifelse(p == "BIOCARTA", "BioCarta", p))))))
}

FIX <- c("Dna" = "DNA", "G2m" = "G2M", "Myc" = "MYC", "E2f" = "E2F",
         "Integrin1" = "Integrin-1", "Cellcycle" = "Cell Cycle")
pretty <- function(x) {
  y <- sub("^(HALLMARK|GOBP|REACTOME|WP|PID|BIOCARTA)_", "", x)
  y <- tools::toTitleCase(tolower(gsub("_", " ", y)))
  y <- sub("^Epithelial Mesenchymal Transition$", "EMT", y)
  y <- sub("Epithelial to Mesenchymal Transition", "Epithelial-to-Mesenchymal Transition", y)
  for (nm in names(FIX)) y <- gsub(paste0("\\b", nm, "\\b"), FIX[[nm]], y, perl = TRUE)
  y <- gsub("Tgf Beta", "TGF-beta", y)
  gsub("Cell Cell Adhesion", "Cell-Cell Adhesion", y)
}

genes_of <- function(h) {
  if (is.na(h) || !nzchar(h)) return(character(0))
  sort(unique(strsplit(h, ";", fixed = TRUE)[[1]]))
}

rows <- lapply(seq_len(nrow(d)), function(i) {
  r <- d[i, ]
  g <- genes_of(r$hits)
  stopifnot(length(g) == r$overlap)
  row <- list(
    set         = pretty(r$term),
    setId       = r$term,
    source      = src_of(r$term),
    theme       = THEME_LABEL[[r$theme]],
    slot        = if (r$theme == "proliferation") SLOT_PROLIF else SLOT_EMT,
    setSize     = as.integer(r$set_size),
    overlap     = as.integer(r$overlap),
    expected    = round(as.numeric(r$expected), 2),
    fold        = round(as.numeric(r$fold), 2),
    depleted    = isTRUE(r$overlap > 0 && as.numeric(r$fold) < 1),
    p           = as.numeric(r$pvalue),
    padj        = as.numeric(r$padj),
    significant = isTRUE(as.numeric(r$padj) < PADJ_CUT),
    genes       = I(g))
  if (r$overlap == 0) row$fold <- NULL
  row
})

ord <- order(factor(vapply(rows, `[[`, character(1), "theme"),
                    levels = unname(THEME_LABEL)),
             -vapply(rows, function(r) if (is.null(r$fold)) -Inf else r$fold, numeric(1)))
rows <- rows[ord]

net    <- readRDS(file.path(DIR_NET, "tcga_network.rds"))
params <- net$params
n_hub  <- sum(net$hub$is_hub)
universe <- length(readRDS(file.path(DIR_NET, "expressed_symbols.rds")))

out <- list(
  title    = "Hub gene enrichment: proliferation vs EMT",
  subtitle = sprintf(paste("Fold enrichment of each gene set among the %d hub genes of the",
                           "TCGA-BRCA DEG co-expression network (n = %s tumours), against all",
                           "%s expressed genes."),
                     n_hub, format(params$n_tumours, big.mark = ","),
                     format(universe, big.mark = ",")),
  unit     = "fold enrichment over expectation",
  axis     = list(scale = "log", reference = 1,
                  referenceLabel = "no enrichment"),
  themes   = list(
    list(label = THEME_LABEL[["proliferation"]],  slot = SLOT_PROLIF),
    list(label = THEME_LABEL[["EMT / motility"]], slot = SLOT_EMT)),
  test     = list(
    name      = "Hypergeometric (one-sided, over-representation)",
    adjust    = "Benjamini-Hochberg across the sets shown",
    query     = n_hub,
    universe  = universe,
    tumours   = params$n_tumours),
  thresholds = list(padj = PADJ_CUT),
  network  = list(
    degPadj  = params$PADJ_MAX,
    degLfc   = params$LFC_MIN,
    rhoMin   = params$RHO_MIN,
    hubs     = n_hub),
  summary  = list(
    sets              = length(rows),
    significant       = sum(vapply(rows, `[[`, logical(1), "significant")),
    proliferationSets = sum(vapply(rows, `[[`, character(1), "theme") == THEME_LABEL[["proliferation"]]),
    emtSets           = sum(vapply(rows, `[[`, character(1), "theme") == THEME_LABEL[["EMT / motility"]]),
    maxFold           = max(vapply(rows, function(r) if (is.null(r$fold)) 0 else r$fold, numeric(1))),
    zeroOverlapSets   = sum(vapply(rows, `[[`, integer(1), "overlap") == 0L)),
  rows = rows)

f_json <- file.path(DIR_EXP, "hub_enrichment.json")
write(toJSON(out, auto_unbox = TRUE, digits = NA, pretty = 2, na = "null"), f_json)
msg("wrote ", f_json, "  (", length(rows), " sets)")

vals <- do.call(rbind, lapply(rows, function(r) data.frame(
  set = r$set, set_id = r$setId, source = r$source, theme = r$theme,
  set_size = r$setSize, overlap = r$overlap, expected = r$expected,
  fold = if (is.null(r$fold)) NA_real_ else r$fold,
  pvalue = r$p, padj = r$padj, significant = r$significant,
  genes = paste(unclass(r$genes), collapse = ";"),
  stringsAsFactors = FALSE)))
f_csv <- file.path(DIR_EXP, "hub_enrichment_values.csv")
write.csv(vals, f_csv, row.names = FALSE)
msg("wrote ", f_csv)

chk <- fromJSON(f_json, simplifyVector = FALSE)
stopifnot(length(chk$rows) == nrow(d))
stopifnot(chk$test$query == n_hub, chk$test$universe == universe)
seen_theme <- character(0)
for (r in chk$rows) {
  stopifnot(nzchar(r$set), nchar(r$set) <= nchar(r$setId) + 8,
            nzchar(r$source), nzchar(r$theme),
            r$setSize > 0, r$overlap >= 0, r$overlap <= r$setSize,
            is.numeric(r$padj), r$padj >= 0, r$padj <= 1,
            is.numeric(r$p), r$p >= 0, r$p <= 1,
            r$significant == (r$padj < PADJ_CUT),
            r$slot %in% c(SLOT_PROLIF, SLOT_EMT),
            length(r$genes) == r$overlap)
  if (r$overlap == 0) stopifnot(is.null(r$fold))
  else stopifnot(is.numeric(r$fold), r$fold > 0,
                 r$depleted == (r$fold < 1))
  seen_theme <- union(seen_theme, r$theme)
}
stopifnot(setequal(seen_theme, unname(THEME_LABEL)))
stopifnot(!any(is.na(unlist(Filter(Negate(is.null), chk$summary)))))

by_id <- setNames(chk$rows, vapply(chk$rows, `[[`, character(1), "setId"))
stopifnot(by_id[["HALLMARK_G2M_CHECKPOINT"]]$fold > 20,
          by_id[["HALLMARK_G2M_CHECKPOINT"]]$significant,
          by_id[["HALLMARK_E2F_TARGETS"]]$fold > 20,
          by_id[["GOBP_WOUND_HEALING"]]$depleted,
          !by_id[["GOBP_WOUND_HEALING"]]$significant,
          by_id[["HALLMARK_TGF_BETA_SIGNALING"]]$overlap == 0L,
          by_id[["REACTOME_INTEGRIN_SIGNALING"]]$overlap == 0L)
msg("assertions passed")

pkg <- function(p) as.character(packageVersion(p))
sig <- vapply(chk$rows, `[[`, logical(1), "significant")
prolif <- vapply(chk$rows, `[[`, character(1), "theme") == THEME_LABEL[["proliferation"]]

row_line <- function(r) sprintf("| `%s` | %s | %d / %d | %s | %s |",
  r$setId, r$theme, r$overlap, chk$test$query,
  if (is.null(r$fold)) "—" else sprintf("%.2f×", r$fold),
  if (r$significant) sprintf("%.2g", r$padj) else sprintf("%.2g (ns)", r$padj))

prov <- c(
"# Provenance — hub gene enrichment figure data", "",
sprintf("Generated %s by `scripts/27_export_enrichment_json.R`.", format(Sys.Date(), "%Y-%m-%d")), "",
"Web version of `plots/hub_proliferation_vs_emt.png`.", "",
"## Regenerate", "", "```bash",
"Rscript scripts/19_tcga_hub_genes.R          # network + hub genes",
"Rscript scripts/21_emt_vs_proliferation.R    # directed test (writes the CSV read here)",
"Rscript scripts/27_export_enrichment_json.R  # this file",
"```", "",
"Nothing is recomputed in the export. Every value is read from",
"`exports/validation/hub_enrichment/directed_emt_vs_proliferation.csv`, which is the same",
"file `scripts/23_plot_prolif_vs_emt.R` draws, so the interactive panel and the publication",
"figure cannot drift apart.", "",
"## What the panel shows", "",
sprintf("Hub genes are the top %d nodes by degree in a co-expression network built from TCGA-BRCA", n_hub),
sprintf("alone: nodes are genes with padj < %.2g and |log2FC| > %.1f in tumour vs normal; edges are",
        params$PADJ_MAX, params$LFC_MIN),
sprintf("|Spearman rho| >= %.2f across %s primary tumours. No external interaction database",
        params$RHO_MIN, format(params$n_tumours, big.mark = ",")),
"(STRING, BioGRID) contributes — that is what restricts the evidence to TCGA. The gene set",
"definitions are of course external; they are the vocabulary of the answer, not a source of",
"expression evidence.", "",
sprintf("Test: hypergeometric over-representation of each set among the %d hubs against the %s",
        n_hub, format(universe, big.mark = ",")),
"genes expressed in these tumours, Benjamini-Hochberg corrected across the sets shown.", "",
"## Every tested set ships, passing or not", "",
sprintf("%d of %d sets reach FDR < %.2g. The remainder are **kept** and carry `significant: false`.",
        sum(sig), length(sig), PADJ_CUT),
"This is deliberate and load-bearing: the result of this analysis is that motility is absent",
"from the hub genes, and a renderer that filtered to significant terms would delete exactly",
"that finding. Two fields exist to let the front end render the negatives correctly:", "",
"- `significant` — `false` rows should be drawn faded, not dropped.",
"- `fold` is **absent** when `overlap` is 0. Those sets contain no hub genes at all and must be",
"  drawn with no bar. Clamping them to a small positive number renders a long leftward bar on",
"  a log axis, which reads as strong depletion when the truth is no overlap.",
"- `depleted` is `true` where fold < 1. GO wound healing — the direct analogue of a scratch",
"  assay — is depleted, and that is a result, not a rounding artifact.", "",
"## Fields", "",
"| Field | Meaning |", "|---|---|",
"| `set` | display label, namespace prefix stripped |",
"| `setId` | original MSigDB set name, the stable key |",
"| `source` | Hallmark, GO:BP, Reactome, WikiPathways, PID, BioCarta |",
"| `theme` | Proliferation or EMT / motility |",
"| `slot` | colour slot (see below) |",
"| `setSize` | genes in the set that are expressed in TCGA-BRCA |",
"| `overlap` | hub genes in the set |",
"| `expected` | overlap expected by chance |",
"| `fold` | overlap / expected; **absent** when overlap is 0 |",
"| `depleted` | fold < 1 |",
"| `p`, `padj` | hypergeometric p, BH-adjusted |",
"| `significant` | padj < 0.05 |",
"| `genes` | every overlapping hub gene symbol, for the tooltip |", "",
"## Colour slots", "",
sprintf("Proliferation = slot %d, EMT / motility = slot %d.", SLOT_PROLIF, SLOT_EMT), "",
"`cohort_boxes.json` fixes 3 = Normal and 8 = Tumour; `subtype_composition.json` uses 1-5 for",
"the PAM50 subtypes and already reuses slot 3 (Basal-like) against that Normal meaning. Slot",
"integers are therefore scoped to a figure family on this site, not global. Pathway enrichment",
"is a third, unrelated family — no cohort and no subtype on either axis — so slots 1 and 2 are",
"reused here on purpose, to keep the web panel the same blue/orange as the publication figure",
"rather than drifting to a different pair of colours for the same two themes.", "",
"## Sets exported", "",
"| Set | Theme | Overlap | Fold | FDR |", "|---|---|---|---|---|",
vapply(chk$rows, row_line, character(1)), "",
if (length(unresolved))
  c(sprintf("## Not exported (%d)", length(unresolved)), "",
    "These names did not resolve in the installed MSigDB release, so they carry no statistics.",
    "They are a lookup failure, not a measured negative, and are dropped rather than shipped as",
    "empty bars:", "",
    paste0("- `", unresolved, "`"), "")
else c("All tested sets resolved in the installed MSigDB release.", ""),
"## Outputs", "",
"| File | Contents |", "|---|---|",
"| `exports/hub_enrichment.json` | D — the panel |",
"| `exports/hub_enrichment_values.csv` | D — one row per set, genes as a `;` list |",
"| `exports/hub_enrichment_provenance.md` | This file |", "",
"## Software", "",
sprintf("- %s", R.version.string),
sprintf("- jsonlite %s", pkg("jsonlite")),
sprintf("- msigdbr %s (set definitions, upstream)", pkg("msigdbr")),
sprintf("- clusterProfiler %s (upstream)", pkg("clusterProfiler")),
sprintf("- Platform: %s", R.version$platform), "",
sprintf("Generated %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

f_prov <- file.path(DIR_EXP, "hub_enrichment_provenance.md")
writeLines(prov, f_prov)
msg("wrote ", f_prov)
msg("done.")
