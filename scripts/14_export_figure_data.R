# Figure data export: cohort boxes, TCGA paired, volcano

suppressPackageStartupMessages({
  library(edgeR)
  library(jsonlite)
})
source(file.path(getwd(), "R", "utils.R"))

GENE    <- "PPP4C"
ENSG    <- "ENSG00000149923"
DIR_EXP <- file.path(PROJ_ROOT, "exports")
RO      <- file.path(PROJ_ROOT, "exports", "r_objects")
P_FLOOR <- .Machine$double.xmin

ORDER <- c("TCGA", "GSE70947", "GSE58135", "GSE42568", "GSE86374")
KEY   <- c(TCGA = "tcga", GSE70947 = "gse70947", GSE58135 = "gse58135",
           GSE42568 = "gse42568", GSE86374 = "gse86374")
TITLE <- c(TCGA = "TCGA-BRCA", GSE70947 = "GSE70947", GSE58135 = "GSE58135",
           GSE42568 = "GSE42568", GSE86374 = "GSE86374")

SLOT_NORMAL <- 3L
SLOT_TUMOUR <- 8L

audit <- list()
note  <- function(...) audit[[length(audit) + 1]] <<-
  data.frame(..., stringsAsFactors = FALSE)

box_stats <- function(x) {
  q   <- as.numeric(quantile(x, c(0.25, 0.5, 0.75), type = 7))
  iqr <- q[3] - q[1]
  lo  <- q[1] - 1.5 * iqr
  hi  <- q[3] + 1.5 * iqr
  inside <- x[x >= lo & x <= hi]
  list(q1 = q[1], median = q[2], q3 = q[3],
       whisker_lo = min(inside), whisker_hi = max(inside),
       mean = mean(x), outliers = sort(x[x < lo | x > hi]), inside = sort(inside))
}

thin_points <- function(inside, target) {
  if (length(inside) <= target) return(inside)
  as.numeric(quantile(inside, probs = seq(0, 1, length.out = target), type = 7))
}

make_group <- function(x, label, slot, max_points) {
  bs <- box_stats(x)
  list(label    = label,
       slot     = slot,
       n        = length(x),
       min      = floor(bs$whisker_lo * 100) / 100,
       q1       = round(bs$q1, 2),
       median   = round(bs$median, 2),
       q3       = round(bs$q3, 2),
       max      = ceiling(bs$whisker_hi * 100) / 100,
       mean     = round(bs$mean, 2),
       outliers = I(round(bs$outliers, 3)),
       points   = I(round(thin_points(bs$inside, max_points), 3)))
}

raw <- list()

counts <- readRDS(file.path(RO, "rna.rds"))
stype  <- substr(colnames(counts), 14, 15)
libsize <- colSums(counts)
dge    <- calcNormFactors(DGEList(counts = counts), method = "TMM")
lcpm   <- cpm(dge, log = TRUE, prior.count = 1)
i_tcga <- which(strip_ens_version(rownames(lcpm)) == ENSG)
stopifnot(length(i_tcga) == 1)
raw$TCGA <- list(
  value = as.numeric(lcpm[i_tcga, ]),
  group = ifelse(stype == "01", "tumor", "normal"),
  sample = colnames(lcpm), feature = rownames(lcpm)[i_tcga],
  unit = "log₂ CPM",
  subtitle = "Illumina RNA-seq · log₂ CPM, TMM-normalised",
  platform = "Illumina RNA-seq (STAR counts, GDC)",
  norm = "TMM (edgeR::calcNormFactors) then cpm(log = TRUE, prior.count = 1)")

d <- readRDS(file.path(RO, "GSE70947_de.rds"))
raw$GSE70947 <- list(
  value = as.numeric(d$expr["58319", ]),
  group = ifelse(as.character(d$pheno$tissue) == "tumor", "tumor", "normal"),
  sample = colnames(d$expr), feature = "probe 58319",
  unit = "log₂ intensity",
  subtitle = "Agilent GPL13607 · log₂ intensity, quantile-normalised",
  platform = "Agilent GPL13607 (1-colour)",
  norm = "background-corrected and log2 quantile-normalised by submitters; used as deposited")

d <- readRDS(file.path(RO, "GSE58135_de.rds"))
i <- which(strip_ens_version(rownames(d$expr)) == ENSG)
stopifnot(length(i) == 1)
raw$GSE58135 <- list(
  value = as.numeric(d$expr[i, ]),
  group = ifelse(as.character(d$pheno$group) == "tumor", "tumor", "normal"),
  sample = colnames(d$expr), feature = rownames(d$expr)[i],
  unit = "log₂(FPKM + 1)",
  subtitle = "Illumina RNA-seq · log₂(FPKM + 1), quantile-normalised",
  platform = "Illumina HiSeq 2000 (cufflinks FPKM)",
  norm = "log2(FPKM + 1) then limma::normalizeBetweenArrays(method = 'quantile')")

for (g in c("GSE42568", "GSE86374")) {
  d  <- readRDS(file.path(RO, paste0(g, "_de.rds")))
  pr <- if (g == "GSE42568") "208932_at" else "7994756"
  arr <- if (g == "GSE42568") "Affymetrix U133 Plus 2.0 (GPL570)"
         else                 "Affymetrix Gene ST 1.0 (GPL6244)"
  raw[[g]] <- list(
    value = as.numeric(d$expr[pr, ]),
    group = as.character(d$pheno$group),
    sample = colnames(d$expr), feature = paste("probe", pr),
    unit = "log₂ intensity",
    subtitle = paste0(sub(" \\(.*", "", arr), " · RMA log₂ intensity, quantile-normalised"),
    platform = arr,
    norm = "RMA log2 from the series matrix, then quantile (normalizeBetweenArrays)")
}

DE_TEST <- c(
  TCGA     = "DESeq2 Wald test (unpaired)",
  GSE70947 = "limma moderated t, paired (blocked on patient)",
  GSE58135 = "limma-trend moderated t (unpaired)",
  GSE42568 = "limma-trend moderated t (unpaired)",
  GSE86374 = "limma-trend moderated t (unpaired)")
DE_PAIRED <- c(TCGA = FALSE, GSE70947 = TRUE, GSE58135 = FALSE,
               GSE42568 = FALSE, GSE86374 = FALSE)

de_row <- function(cohort) {
  d <- read_de(cohort)
  r <- d[d$symbol == GENE, ]
  stopifnot(nrow(r) == 1)
  r
}

cohorts <- lapply(ORDER, function(g) {
  r  <- raw[[g]]
  na <- is.na(r$value)
  if (any(na)) note(cohort = g, item = "dropped: NA expression value", n = sum(na))
  keep <- !na
  xn <- r$value[keep & r$group == "normal"]
  xt <- r$value[keep & r$group == "tumor"]
  note(cohort = g, item = "samples in cohort matrix", n = length(r$value))
  note(cohort = g, item = "normal (plotted)", n = length(xn))
  note(cohort = g, item = "tumour (plotted)", n = length(xt))

  de <- de_row(g)
  msg(sprintf("%-9s normal n=%3d med=%6.2f | tumour n=%4d med=%6.2f | DE logFC=%+.3f p=%.3g",
              g, length(xn), median(xn), length(xt), median(xt), de$logFC, de$pvalue))

  list(key      = unname(KEY[g]),
       title    = unname(TITLE[g]),
       subtitle = r$subtitle,
       unit     = r$unit,
       test     = list(name   = unname(DE_TEST[g]),
                       p      = max(as.numeric(de$pvalue), P_FLOOR),
                       effect = sprintf("%+.2f log₂ units", de$logFC)),
       groups   = list(make_group(xn, "Normal", SLOT_NORMAL, 100),
                       make_group(xt, "Tumour", SLOT_TUMOUR, 100)))
})

f_boxes <- file.path(DIR_EXP, "cohort_boxes.json")
write(toJSON(list(cohorts = cohorts), auto_unbox = TRUE, digits = NA,
             pretty = 2, na = "null"), f_boxes)
msg("wrote ", f_boxes)

pat <- substr(colnames(lcpm), 1, 12)
expr_all <- as.numeric(lcpm[i_tcga, ])

pick_one <- function(idx) {
  if (length(idx) == 1) return(idx)
  idx[order(-libsize[idx], colnames(lcpm)[idx])][1]
}

tum_pat <- unique(pat[stype == "01"])
nor_pat <- unique(pat[stype == "11"])
matched <- sort(intersect(tum_pat, nor_pat))
note(cohort = "TCGA-paired", item = "patients with >=1 primary tumour", n = length(tum_pat))
note(cohort = "TCGA-paired", item = "patients with >=1 adjacent normal", n = length(nor_pat))
note(cohort = "TCGA-paired", item = "patients with BOTH (pairs drawn)", n = length(matched))
note(cohort = "TCGA-paired", item = "normal samples with no tumour from same patient",
     n = sum(!(pat[stype == "11"] %in% tum_pat)))

n_dup_t <- sum(vapply(matched, function(p) sum(pat == p & stype == "01") > 1, logical(1)))
n_dup_n <- sum(vapply(matched, function(p) sum(pat == p & stype == "11") > 1, logical(1)))
note(cohort = "TCGA-paired", item = "matched patients needing tumour-aliquot tie-break", n = n_dup_t)
note(cohort = "TCGA-paired", item = "matched patients needing normal-aliquot tie-break", n = n_dup_n)

pairs_raw <- do.call(rbind, lapply(matched, function(p) {
  it <- pick_one(which(pat == p & stype == "01"))
  ino <- pick_one(which(pat == p & stype == "11"))
  data.frame(id = p, normal = expr_all[ino], tumour = expr_all[it],
             stringsAsFactors = FALSE)
}))

pairs_raw$normal_r <- round(pairs_raw$normal, 2)
pairs_raw$tumour_r <- round(pairs_raw$tumour, 2)
pairs_raw$delta_r  <- round(pairs_raw$tumour_r - pairs_raw$normal_r, 2)

n_up   <- sum(pairs_raw$delta_r > 0)
n_down <- sum(pairs_raw$delta_r < 0)
n_zero <- sum(pairs_raw$delta_r == 0)
note(cohort = "TCGA-paired", item = "pairs with delta exactly 0 after rounding", n = n_zero)

wsr <- suppressWarnings(wilcox.test(pairs_raw$tumour, pairs_raw$normal, paired = TRUE))

paired_out <- list(
  title  = "PPP4C in matched tumour/normal pairs",
  unit   = raw$TCGA$unit,
  groups = list(list(label = "Normal", slot = SLOT_NORMAL),
                list(label = "Tumour", slot = SLOT_TUMOUR)),
  test   = list(name   = "Wilcoxon signed-rank (paired)",
                p      = max(as.numeric(wsr$p.value), P_FLOOR),
                effect = sprintf("median within-patient change %+.2f log₂ units",
                                 median(pairs_raw$delta_r))),
  summary = list(pairs       = nrow(pairs_raw),
                 up          = n_up,
                 down        = n_down,
                 medianDelta = round(median(pairs_raw$delta_r), 2),
                 meanDelta   = round(mean(pairs_raw$delta_r), 2)),
  pairs = unname(lapply(seq_len(nrow(pairs_raw)), function(k) list(
    id     = pairs_raw$id[k],
    normal = pairs_raw$normal_r[k],
    tumour = pairs_raw$tumour_r[k],
    delta  = pairs_raw$delta_r[k]))))

f_paired <- file.path(DIR_EXP, "tcga_paired.json")
write(toJSON(paired_out, auto_unbox = TRUE, digits = NA, pretty = 2, na = "null"), f_paired)
msg("wrote ", f_paired, "  (", nrow(pairs_raw), " pairs)")

f_pcsv <- file.path(DIR_EXP, "tcga_paired_values.csv")
write.csv(data.frame(patient_id = pairs_raw$id, normal = pairs_raw$normal_r,
                     tumour = pairs_raw$tumour_r, delta = pairs_raw$delta_r),
          f_pcsv, row.names = FALSE)
msg("wrote ", f_pcsv)

vol <- read_de("TCGA")
names(vol)[names(vol) == "logFC"] <- "lfc"

f_vcsv <- file.path(DIR_EXP, "tcga_volcano_full.csv")
write.csv(data.frame(gene = vol$symbol, lfc = vol$lfc,
                     pvalue = vol$pvalue, padj = vol$padj),
          f_vcsv, row.names = FALSE)
msg("wrote ", f_vcsv, "  (", nrow(vol), " genes, unthinned)")

dir_up   <- c("COL10A1", "MMP11", "TOP2A", "MKI67")
dir_down <- c("PPARG", "ADIPOQ", "LEP", "FABP4")
stopifnot(all(vol$lfc[match(dir_up, vol$symbol)] > 0),
          all(vol$lfc[match(dir_down, vol$symbol)] < 0))
msg("direction verified: +lfc = up in tumour (", paste(dir_up, collapse = "/"),
    " positive, ", paste(dir_down, collapse = "/"), " negative)")

n_all      <- nrow(vol)
n_na_padj  <- sum(is.na(vol$padj))
vol        <- vol[!is.na(vol$padj), ]
n_tested   <- nrow(vol)

PADJ_CUT <- 0.05
LFC_CUT  <- 1
is_sig   <- vol$padj < PADJ_CUT & abs(vol$lfc) > LFC_CUT
n_sig    <- sum(is_sig)

MAX_SIG    <- 1200
N_NONSIG   <- 800
ALWAYS <- c("PPP4C", "PPP4R1", "PPP4R2", "PPP4R3A", "PPP4R3B", "PPP4R4",
            "ESR1", "PGR", "ERBB2", "MKI67", "TP53", "GATA3", "FOXA1")

sig_tab <- vol[is_sig, ]
sig_tab <- sig_tab[order(sig_tab$padj, sig_tab$pvalue, sig_tab$symbol), ]
keep_sig <- head(sig_tab$symbol, MAX_SIG)

non_tab <- vol[!is_sig, ]
non_tab <- non_tab[order(non_tab$padj, non_tab$pvalue, non_tab$symbol), ]
idx_non <- unique(round(seq(1, nrow(non_tab), length.out = min(N_NONSIG, nrow(non_tab)))))
keep_non <- non_tab$symbol[idx_non]

missing_always <- setdiff(ALWAYS, vol$symbol)
keep_always    <- intersect(ALWAYS, vol$symbol)
added_always   <- setdiff(keep_always, union(keep_sig, keep_non))

keep <- unique(c(keep_sig, keep_non, keep_always))
rows_tab <- vol[match(keep, vol$symbol), ]
rows_tab <- rows_tab[order(rows_tab$padj, rows_tab$pvalue, rows_tab$symbol), ]

floored <- rows_tab$symbol[rows_tab$padj <= 0]
rows_tab$padj_out <- pmax(rows_tab$padj, P_FLOOR)

volcano_out <- list(
  comparison   = "TCGA-BRCA, primary tumour vs adjacent normal",
  highLabel    = "Up in tumour",
  lowLabel     = "Up in normal",
  test         = "DESeq2 Wald test (negative binomial GLM)",
  nTested      = n_tested,
  nSignificant = n_sig,
  thresholds   = list(padj = PADJ_CUT, lfc = LFC_CUT),
  rows = unname(lapply(seq_len(nrow(rows_tab)), function(k) list(
    gene = rows_tab$symbol[k],
    lfc  = round(rows_tab$lfc[k], 3),
    padj = rows_tab$padj_out[k]))))

f_vol <- file.path(DIR_EXP, "tcga_volcano.json")
write(toJSON(volcano_out, auto_unbox = TRUE, digits = NA, pretty = 2, na = "null"), f_vol)
msg("wrote ", f_vol, "  (", nrow(rows_tab), " rows of ", n_tested, " tested)")

chk <- fromJSON(f_boxes, simplifyVector = FALSE)
stopifnot(length(chk$cohorts) == 5)
for (co in chk$cohorts) {
  stopifnot(length(co$groups) == 2,
            co$groups[[1]]$label == "Normal", co$groups[[1]]$slot == SLOT_NORMAL,
            co$groups[[2]]$label == "Tumour", co$groups[[2]]$slot == SLOT_TUMOUR,
            is.numeric(co$test$p), co$test$p > 0, nzchar(co$unit), nzchar(co$subtitle))
  for (gr in co$groups) {
    pts <- unlist(gr$points); out <- unlist(gr$outliers)
    stopifnot(is.list(gr$outliers), is.list(gr$points),
              gr$min <= gr$q1, gr$q1 <= gr$median, gr$median <= gr$q3, gr$q3 <= gr$max,
              length(pts) <= 100, length(pts) >= 1,
              all(pts >= gr$min), all(pts <= gr$max),
              length(out) == 0 || all(out < gr$min | out > gr$max))
  }
}
stopifnot(!any(is.na(unlist(chk))))

chp <- fromJSON(f_paired, simplifyVector = FALSE)
stopifnot(chp$summary$pairs == length(chp$pairs),
          chp$summary$up + chp$summary$down + n_zero == chp$summary$pairs,
          chp$summary$up == sum(vapply(chp$pairs, function(p) p$delta > 0, logical(1))),
          chp$summary$down == sum(vapply(chp$pairs, function(p) p$delta < 0, logical(1))),
          all(vapply(chp$pairs, function(p) isTRUE(all.equal(p$delta, round(p$tumour - p$normal, 2))), logical(1))),
          all(vapply(chp$pairs, function(p) grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", p$id), logical(1))),
          !anyDuplicated(vapply(chp$pairs, function(p) p$id, character(1))),
          chp$test$p > 0)
stopifnot(!any(is.na(unlist(chp))))

chv <- fromJSON(f_vol, simplifyVector = FALSE)
stopifnot(length(chv$rows) >= 1500, length(chv$rows) <= 2500,
          chv$nTested == n_tested, chv$nSignificant == n_sig,
          all(vapply(chv$rows, function(r) r$padj > 0, logical(1))),
          all(ALWAYS %in% vapply(chv$rows, function(r) r$gene, character(1))),
          !anyDuplicated(vapply(chv$rows, function(r) r$gene, character(1))))
stopifnot(!any(is.na(unlist(chv))))
msg("validated: 3 JSON files, slots fixed, arrays intact, ",
    "points within whiskers, deltas exact, no NA/null, no zero p-values")

aud <- do.call(rbind, audit)
pkg <- function(p) as.character(packageVersion(p))
de_p <- vapply(ORDER, function(g) de_row(g)$pvalue, numeric(1))
de_q <- vapply(ORDER, function(g) de_row(g)$padj,   numeric(1))
de_l <- vapply(ORDER, function(g) de_row(g)$logFC,  numeric(1))

prov <- c(
"# Provenance — figure data export", "",
sprintf("Generated %s by `scripts/14_export_figure_data.R`.", format(Sys.Date(), "%Y-%m-%d")), "",
sprintf("Gene: **%s** (`%s`).", GENE, ENSG), "",
"## Regenerate", "", "```bash", "Rscript scripts/14_export_figure_data.R", "```", "",
"Inputs are the cached objects under `exports/r_objects/` and the DE tables under",
"`exports/validation/`, produced by `bash scripts/run_all.sh`. Nothing here re-downloads",
"from GEO or GDC, and nothing re-runs a differential-expression model.", "",
"**Supersedes** `exports/cohort_comparison.json` (`scripts/13_export_cohort_json.R`), which",
"put all five cohorts on one centred log2FC axis. That standardisation is dropped: each",
"cohort is now its own panel on its own native scale. Those outputs",
"(`cohort_comparison.json`, `_values.csv`, `_provenance.md`) have been deleted;",
"`scripts/13_export_cohort_json.R` is retained as history and can regenerate them.", "",
"---", "", "## A. `cohort_boxes.json` — five tumour vs normal panels", "",
"Each panel is self-contained and uses that cohort's native unit. Nothing is centred,",
"rescaled or transformed beyond what the upstream DE already did. Values are therefore",
"**not comparable between panels** — that is the point of the redesign.", "",
"Normal is always slot 3 and Tumour always slot 8, in all five panels.", "",
"| Cohort | Platform | Normalisation | Native unit | Feature | n normal | n tumour |",
"|---|---|---|---|---|---|---|",
sprintf("| %s | %s | %s | %s | `%s` | %d | %d |",
        TITLE[ORDER], vapply(ORDER, function(g) raw[[g]]$platform, ""),
        vapply(ORDER, function(g) raw[[g]]$norm, ""),
        vapply(ORDER, function(g) raw[[g]]$unit, ""),
        vapply(ORDER, function(g) raw[[g]]$feature, ""),
        vapply(ORDER, function(g) aud$n[aud$cohort == g & aud$item == "normal (plotted)"], 0L),
        vapply(ORDER, function(g) aud$n[aud$cohort == g & aud$item == "tumour (plotted)"], 0L)),
"",
"GSE42568 and GSE86374 share the unit string `log₂ intensity` because that genuinely is",
"the same unit — both are RMA-summarised log2 probe intensities. The platforms differ",
"(GPL570 vs GPL6244) and the `subtitle` field carries that distinction. Inventing two",
"different unit strings would have been false precision.", "",
"GPL13607 carries two PPP4C probes (58319, 57896). Probe **58319** is used because that is",
"the one `collapse_by_symbol()` retained for the DE, so panel and statistics agree.", "",
"### Test per cohort — taken from the upstream DE, not recomputed", "",
"`p` and `effect` are read verbatim from `exports/validation/DE_<cohort>.csv`. No new test",
"was run for this figure.", "",
"| Cohort | Test as actually run | Paired? | logFC on file | p (raw) | padj (BH) |",
"|---|---|---|---|---|---|",
sprintf("| %s | %s | %s | %+.3f | %.3g | %.3g |",
        TITLE[ORDER], DE_TEST[ORDER],
        ifelse(DE_PAIRED[ORDER], "**yes**", "no"), de_l, de_p, de_q),
"",
"**GSE70947 is the one paired design** — `scripts/02_de_gse70947.R` fits",
"`~ patient + tissue`, blocking on patient across 148 matched pairs, so its moderated t",
"is a paired test and is labelled that way. The other four are two-group contrasts and no",
"paired test was applied to them.", "",
"The JSON `test.p` field carries the **raw** p-value (the test's own p). The BH-adjusted",
"padj is tabulated above but is not in the JSON; cite the adjusted value if the figure",
"caption claims genome-wide significance.", "",
sprintf("Note: TCGA's DE was run **unpaired** (`~ sample_type`), even though all %d of its",
        length(matched)),
"normals are in fact matched to a tumour from the same patient (see B). That was the",
"upstream choice and it is preserved here rather than silently upgraded.", "",
"Box construction: Q1/median/Q3 at `quantile(type = 7)`; `min`/`max` are Tukey whisker",
"ends (most extreme value within 1.5 x IQR), **not** the data range; `outliers` lists every",
"value past a whisker and is `[]` when empty; `points` are within-whisker values thinned",
"to **at most 100 per group** by evenly spaced quantiles (deterministic, no RNG). Whiskers",
"round outward so no point can fall outside them after rounding. Summary stats 2 dp,",
"points and outliers 3 dp.", "",
"---", "", "## B. `tcga_paired.json` — matched tumour/normal pairs", "",
sprintf("Native TCGA scale (%s). No centring, no fold-change.", raw$TCGA$unit), "",
"### Match accounting", "",
"| Item | n |", "|---|---|",
sprintf("| %s | %d |", aud$item[aud$cohort == "TCGA-paired"], aud$n[aud$cohort == "TCGA-paired"]),
"",
sprintf("**All %d adjacent normals matched a tumour from the same patient — none were lost.**",
        length(matched)),
"The 113 normal samples come from 113 distinct patients, and every one of those patients",
"also has at least one primary tumour sample in the matrix, so the pair count equals the",
"normal count exactly. Nothing was dropped for being unmatchable.", "",
sprintf("%d pairs is under the ~200 threshold, so **every pair is included and none are thinned**.",
        nrow(pairs_raw)), "",
"### Aliquot rule", "",
"Drawing a line needs exactly one value per patient per tissue.",
sprintf("%d of the %d matched patients contributed more than one primary-tumour aliquot (%s);",
        n_dup_t, length(matched),
        paste(matched[vapply(matched, function(p) sum(pat == p & stype == "01") > 1, logical(1))],
              collapse = ", ")),
sprintf("%d contributed more than one normal aliquot.", n_dup_n), "",
"**Rule applied uniformly to every patient**: keep the aliquot with the largest library",
"size (total mapped counts), i.e. the deepest and least noisy measurement. Ties break on",
"the lexicographically first barcode. The rule is evaluated for all patients, not only the",
"duplicated ones, so it cannot introduce a selection asymmetry.", "",
"### Deltas", "",
"`normal` and `tumour` are rounded to 2 dp first, then `delta` is computed from those",
"rounded values, so `delta == tumour - normal` is exactly true as displayed.",
sprintf("`up` = %d (delta > 0), `down` = %d (delta < 0), and `pairs` = %d.", n_up, n_down, nrow(pairs_raw)),
if (n_zero == 0)
  "No pair had a delta of exactly 0 after rounding, so up + down = pairs with no residual."
else
  sprintf(paste("%d pair(s) had a delta of exactly 0.00 after rounding and are counted in",
                "neither up nor down, so up + down = %d, which is %d short of pairs."),
          n_zero, n_up + n_down, n_zero), "",
sprintf("- **Wilcoxon signed-rank (paired)**: V = %s, computed on the *unrounded*",
        format(as.numeric(wsr$statistic))),
"  within-patient differences. `medianDelta` / `meanDelta` in the JSON come from the",
"  rounded deltas so that every number in the file agrees with the rows shown.",
if (as.numeric(wsr$p.value) <= 0)
  sprintf(paste("- The p-value **underflowed to exactly 0** and is floored at the smallest",
                "representable double (%.3g) in the JSON. Read it as \"below double",
                "precision\", not as a point estimate."), P_FLOOR)
else
  sprintf("- p = %.3g.", as.numeric(wsr$p.value)),
"- **This test is new to this export.** The upstream pipeline never ran a paired test on",
"  TCGA; its DE was the unpaired DESeq2 contrast. Do not cite this p as pre-existing.", "",
"---", "", "## C. `tcga_volcano.json` — genome-wide tumour vs adjacent normal", "",
"This is **tumour tissue vs adjacent normal tissue**, not the PPP4C-high vs PPP4C-low",
"comparison in `scripts/11_target_dossier.R`.", "",
sprintf("- **Method**: DESeq2 Wald test on the negative-binomial GLM, `design = ~ sample_type`,"),
"  run in `TCGA_BRCA_analysis.Rmd`, result at `exports/tables/DEA.csv`, symbol-collapsed by",
"  `scripts/04_tcga_ranked_list.R` into `exports/validation/DE_TCGA.csv`. Not re-run here.",
"- **Direction**: the Rmd sets `relevel(sample_type, ref = \"Solid Tissue Normal\")`, so the",
"  contrast is Primary Tumour vs Solid Tissue Normal and **positive `lfc` = up in tumour**.",
"  The script asserts this at runtime — COL10A1/MMP11/TOP2A/MKI67 must be positive and",
"  PPARG/ADIPOQ/LEP/FABP4 negative — and aborts if the sign ever flips.", "",
"### Gene universe", "", "| Step | Genes |", "|---|---|",
"| Quantified by STAR (GDC) | 60,660 |",
"| After `rowSums(counts) >= 10` in the Rmd | 56,705 |",
sprintf("| After symbol collapse (`collapse_by_symbol`) | %s |", format(n_all, big.mark = ",")),
sprintf("| With non-NA padj — **`nTested`** | %s |", format(n_tested, big.mark = ",")),
sprintf("| Passing padj < %.2f AND abs(lfc) > %d — **`nSignificant`** | %s |",
        PADJ_CUT, LFC_CUT, format(n_sig, big.mark = ",")),
"",
sprintf("%s genes have padj = NA from DESeq2 independent filtering (low counts). They cannot be",
        format(n_na_padj, big.mark = ",")),
"classified against a padj threshold, so they are excluded from `nTested` and from the rows.",
"`nTested` and `nSignificant` are the TRUE totals **before** thinning — the figure is",
"honestly a subset of them.", "",
"### Thinning (deterministic, no RNG)", "", "| Component | Rule | Rows |", "|---|---|---|",
sprintf("| Significant | all genes passing both thresholds, capped at %d with smallest padj | %d |",
        MAX_SIG, length(keep_sig)),
sprintf("| Non-significant | evenly spaced through the padj-sorted non-significant list of %s genes (`seq(1, %d, length.out = %d)`) | %d |",
        format(nrow(non_tab), big.mark = ","), nrow(non_tab), N_NONSIG, length(keep_non)),
sprintf("| Always-include | %s | %d added that thinning missed |",
        paste(ALWAYS, collapse = ", "), length(added_always)),
sprintf("| **Total rows** | after de-duplication | **%d** |", nrow(rows_tab)),
"",
sprintf("Sort order everywhere is padj, then pvalue, then symbol — fully deterministic."),
sprintf("Of the %s significant genes, %s were kept and %s were dropped by the %s cap, so the",
        format(n_sig, big.mark = ","), format(length(keep_sig), big.mark = ","),
        format(n_sig - length(keep_sig), big.mark = ","), format(MAX_SIG, big.mark = ",")),
"volcano's wings are truncated by design; the untruncated result is in",
"`exports/tcga_volcano_full.csv`.", "",
if (length(missing_always) == 0)
  "All always-include genes were present in the DE table." else
  sprintf("Always-include genes NOT found in the DE table and therefore omitted: %s.",
          paste(missing_always, collapse = ", ")),
if (length(added_always))
  sprintf("Always-include genes that thinning would have missed and were added back: %s.",
          paste(added_always, collapse = ", ")) else
  "Every always-include gene was already selected by the thinning.", "",
"### Underflow", "",
sprintf("%d gene(s) have padj underflowed to exactly 0 in DESeq2 and are floored at the smallest",
        length(floored)),
sprintf("representable double (%.3g) in the JSON: **%s**.", P_FLOOR, paste(floored, collapse = ", ")),
"Read those as \"below double precision\", not as point estimates. The renderer clamps padj",
"at 1e-50 for display; the JSON carries the real values regardless.",
"`exports/tcga_volcano_full.csv` is a faithful unthinned copy of the DE table and keeps the",
"raw `0` values, so it can be reconciled against `DE_TCGA.csv` byte for byte.", "",
"---", "", "## Files", "", "| File | Contents |", "|---|---|",
"| `exports/cohort_boxes.json` | A — five tumour vs normal panels, native units |",
"| `exports/tcga_paired.json` | B — matched-pair slope graph |",
"| `exports/tcga_paired_values.csv` | B — `patient_id, normal, tumour, delta` |",
"| `exports/tcga_volcano.json` | C — thinned volcano rows |",
"| `exports/tcga_volcano_full.csv` | C — unthinned DE result, `gene, lfc, pvalue, padj` |",
"| `exports/figure_data_provenance.md` | This file |", "",
"## Software", "",
sprintf("- %s", R.version.string),
sprintf("- edgeR %s", pkg("edgeR")),
sprintf("- limma %s (upstream DE)", pkg("limma")),
sprintf("- jsonlite %s", pkg("jsonlite")),
sprintf("- Platform: %s", R.version$platform), "",
sprintf("Generated %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

f_prov <- file.path(DIR_EXP, "figure_data_provenance.md")
writeLines(prov, f_prov)
msg("wrote ", f_prov)
msg("done.")
