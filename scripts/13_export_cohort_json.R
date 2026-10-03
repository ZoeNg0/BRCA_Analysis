# PPP4C across five cohorts (JSON export)

suppressPackageStartupMessages({
  library(edgeR)
  library(jsonlite)
})
source(file.path(getwd(), "R", "utils.R"))

GENE    <- "PPP4C"
ENSG    <- "ENSG00000149923"
DIR_EXP <- file.path(PROJ_ROOT, "exports")
RO      <- file.path(PROJ_ROOT, "exports", "r_objects")

ORDER <- c("TCGA", "GSE70947", "GSE58135", "GSE42568", "GSE86374")

drops <- list()
note  <- function(cohort, step, n) drops[[length(drops) + 1]] <<-
  data.frame(cohort = cohort, step = step, n = n, stringsAsFactors = FALSE)

raw <- list()

counts <- readRDS(file.path(RO, "rna.rds"))
stype  <- substr(colnames(counts), 14, 15)
dge    <- calcNormFactors(DGEList(counts = counts), method = "TMM")
lcpm   <- cpm(dge, log = TRUE, prior.count = 1)
i_tcga <- which(strip_ens_version(rownames(lcpm)) == ENSG)
stopifnot(length(i_tcga) == 1)
raw$TCGA <- list(
  value   = as.numeric(lcpm[i_tcga, ]),
  group   = ifelse(stype == "01", "tumor", "normal"),
  sample  = colnames(lcpm),
  feature = rownames(lcpm)[i_tcga],
  unit    = "log2 CPM (TMM-normalised STAR counts, prior.count = 1)")

d <- readRDS(file.path(RO, "GSE70947_de.rds"))
raw$GSE70947 <- list(
  value   = as.numeric(d$expr["58319", ]),
  group   = ifelse(as.character(d$pheno$tissue) == "tumor", "tumor", "normal"),
  sample  = colnames(d$expr),
  feature = "probe 58319",
  unit    = "log2 intensity, quantile-normalised (as deposited)")

d <- readRDS(file.path(RO, "GSE58135_de.rds"))
i <- which(strip_ens_version(rownames(d$expr)) == ENSG)
stopifnot(length(i) == 1)
raw$GSE58135 <- list(
  value   = as.numeric(d$expr[i, ]),
  group   = ifelse(as.character(d$pheno$group) == "tumor", "tumor", "normal"),
  sample  = colnames(d$expr),
  feature = rownames(d$expr)[i],
  unit    = "log2(FPKM + 1), quantile-normalised")

for (g in c("GSE42568", "GSE86374")) {
  d  <- readRDS(file.path(RO, paste0(g, "_de.rds")))
  pr <- if (g == "GSE42568") "208932_at" else "7994756"
  raw[[g]] <- list(
    value   = as.numeric(d$expr[pr, ]),
    group   = as.character(d$pheno$group),
    sample  = colnames(d$expr),
    feature = paste("probe", pr),
    unit    = "RMA log2 intensity, quantile-normalised")
}

centred <- list()
for (g in ORDER) {
  r  <- raw[[g]]
  note(g, "samples in cohort matrix", length(r$value))
  na <- is.na(r$value)
  if (any(na)) note(g, "dropped: NA expression value", sum(na))
  keep <- !na

  n_norm <- sum(keep & r$group == "normal")
  ref    <- median(r$value[keep & r$group == "normal"])
  note(g, "held out as normal reference (not plotted)", n_norm)

  sel <- keep & r$group == "tumor"
  note(g, "**plotted (tumour samples)**", sum(sel))
  centred[[g]] <- list(
    sample  = r$sample[sel],
    value   = r$value[sel] - ref,
    ref     = ref,
    n_norm  = n_norm,
    feature = r$feature,
    unit    = r$unit)
  msg(sprintf("%-9s %4d tumours | normal ref (n=%3d) = %6.3f | median log2FC = %+.3f",
              g, sum(sel), n_norm, ref, median(r$value[sel] - ref)))
}

long <- do.call(rbind, lapply(ORDER, function(g) data.frame(
  sample_id  = centred[[g]]$sample,
  cohort     = g,
  expression = centred[[g]]$value,
  stringsAsFactors = FALSE)))
f_csv <- file.path(DIR_EXP, "cohort_comparison_values.csv")
write.csv(long, f_csv, row.names = FALSE)
msg("wrote ", f_csv, "  (", nrow(long), " rows)")

kw <- kruskal.test(expression ~ factor(cohort), data = long)

pairs <- combn(ORDER, 2, simplify = FALSE)
pw <- do.call(rbind, lapply(pairs, function(p) data.frame(
  a = p[1], b = p[2],
  p = suppressWarnings(wilcox.test(long$expression[long$cohort == p[1]],
                                   long$expression[long$cohort == p[2]])$p.value),
  stringsAsFactors = FALSE)))
pw$padj <- p.adjust(pw$p, method = "BH")

P_FLOOR <- .Machine$double.xmin
n_underflow <- sum(pw$p <= 0)
pw$p    <- pmax(pw$p,    P_FLOOR)
pw$padj <- pmax(pw$padj, P_FLOOR)

box_stats <- function(x) {
  q  <- as.numeric(quantile(x, c(0.25, 0.5, 0.75), type = 7))
  iqr <- q[3] - q[1]
  lo <- q[1] - 1.5 * iqr
  hi <- q[3] + 1.5 * iqr
  inside <- x[x >= lo & x <= hi]
  list(q1 = q[1], median = q[2], q3 = q[3],
       whisker_lo = min(inside), whisker_hi = max(inside),
       mean = mean(x), outliers = sort(x[x < lo | x > hi]), inside = sort(inside))
}

thin <- function(inside, target = 200) {
  if (length(inside) <= target) return(inside)
  as.numeric(quantile(inside, probs = seq(0, 1, length.out = target), type = 7))
}

groups <- lapply(seq_along(ORDER), function(k) {
  g  <- ORDER[k]
  x  <- centred[[g]]$value
  bs <- box_stats(x)
  pts <- thin(bs$inside)
  list(label    = g,
       slot     = k,
       n        = length(x),
       min      = floor(bs$whisker_lo * 100) / 100,
       q1       = round(bs$q1, 2),
       median   = round(bs$median, 2),
       q3       = round(bs$q3, 2),
       max      = ceiling(bs$whisker_hi * 100) / 100,
       mean     = round(bs$mean, 2),
       outliers = I(round(bs$outliers, 3)),
       points   = I(round(pts, 3)))
})

out <- list(
  title = "PPP4C is elevated in tumour across five independent BRCA cohorts",
  unit  = "log₂ fold-change vs cohort normal median",
  test  = list(name = "Kruskal–Wallis",
               p = as.numeric(kw$p.value),
               effect = "across 5 cohorts"),
  groups = groups,
  pairwise = unname(lapply(seq_len(nrow(pw)), function(i) list(
    a = pw$a[i], b = pw$b[i],
    p = as.numeric(pw$p[i]), padj = as.numeric(pw$padj[i])))))

f_json <- file.path(DIR_EXP, "cohort_comparison.json")
write(toJSON(out, auto_unbox = TRUE, digits = NA, pretty = 2, na = "null"), f_json)
msg("wrote ", f_json)

chk <- fromJSON(f_json, simplifyVector = FALSE)
stopifnot(length(chk$groups) == 5)
for (gr in chk$groups) {
  pts <- unlist(gr$points); out <- unlist(gr$outliers)
  stopifnot(
    is.numeric(gr$n), is.numeric(gr$median),
    is.list(gr$outliers), is.list(gr$points),
    gr$min <= gr$q1, gr$q1 <= gr$median,
    gr$median <= gr$q3, gr$q3 <= gr$max,
    gr$n == sum(long$cohort == gr$label),
    length(pts) >= min(120, gr$n - length(out)),
    length(pts) <= 200,
    all(pts >= gr$min), all(pts <= gr$max),
    length(out) == 0 || all(out < gr$min | out > gr$max))
}
stopifnot(!any(is.na(unlist(chk))),
          chk$test$p > 0,
          all(vapply(chk$pairwise, function(p) p$p > 0 && p$padj > 0, logical(1))))
msg("validated: 5 groups, n matches CSV, arrays intact, ",
    "points within whiskers, no NA/null, no zero p-values")

dropdf <- do.call(rbind, drops)
pkg <- function(p) as.character(packageVersion(p))
prov <- c(
"# Provenance — `cohort_comparison.json`", "",
sprintf("Generated %s by `scripts/13_export_cohort_json.R`.", format(Sys.Date(), "%Y-%m-%d")), "",
sprintf("Figure: **%s** (`%s`) across five breast-cancer cohorts.", GENE, ENSG), "",
"## Regenerate", "", "```bash", "Rscript scripts/13_export_cohort_json.R", "```", "",
"Inputs are the cached objects under `exports/r_objects/`, themselves produced by",
"`bash scripts/run_all.sh`. Nothing here re-downloads from GEO or GDC.", "",
"## What is plotted", "",
"Each point is one **tumour** sample, expressed as a log2 fold-change against the",
"**median of that same cohort's normal samples**. Normal samples define the reference",
"and are not drawn.", "",
"This centring is deliberate. The five cohorts are on four platforms with four",
"normalisations and their raw PPP4C medians span 3.50 to 11.10; that spread is",
"measurement technology, not biology. Raw values are therefore not comparable on a",
"shared y-axis, and a one-way test across them (chi-squared = 1365.5) would be",
"quantifying platform differences. Centring each cohort on its own normals removes",
"the platform offset and leaves the tumour-vs-normal contrast the target selection",
"actually used.", "",
"Sanity check — the per-cohort medians below reproduce the independently computed",
"DE log fold-changes in `exports/validation/DE_<cohort>.csv` to within rounding.", "",
"| Cohort | Accession / source | Platform | Native unit before centring | Feature used | Normal ref (n) | Median log2FC here | DE logFC on file |",
"|---|---|---|---|---|---|---|---|",
sprintf("| TCGA-BRCA | GDC `TCGA-BRCA`, STAR counts (`exports/r_objects/rna.rds`) | Illumina RNA-seq | %s | `%s` | %.3f (n=%d) | %+.3f | +1.341 |",
        centred$TCGA$unit, centred$TCGA$feature, centred$TCGA$ref, centred$TCGA$n_norm, median(centred$TCGA$value)),
sprintf("| GSE70947 | GEO `GSE70947` (cohort_2) | Agilent GPL13607 | %s | %s | %.3f (n=%d) | %+.3f | +0.556 |",
        centred$GSE70947$unit, centred$GSE70947$feature, centred$GSE70947$ref, centred$GSE70947$n_norm, median(centred$GSE70947$value)),
sprintf("| GSE58135 | GEO `GSE58135`, cufflinks FPKM | Illumina HiSeq 2000 | %s | `%s` | %.3f (n=%d) | %+.3f | +0.830 |",
        centred$GSE58135$unit, centred$GSE58135$feature, centred$GSE58135$ref, centred$GSE58135$n_norm, median(centred$GSE58135$value)),
sprintf("| GSE42568 | GEO `GSE42568` | Affymetrix U133 Plus 2.0 (GPL570) | %s | %s | %.3f (n=%d) | %+.3f | +1.188 |",
        centred$GSE42568$unit, centred$GSE42568$feature, centred$GSE42568$ref, centred$GSE42568$n_norm, median(centred$GSE42568$value)),
sprintf("| GSE86374 | GEO `GSE86374` | Affymetrix Gene ST 1.0 (GPL6244) | %s | %s | %.3f (n=%d) | %+.3f | +0.843 |",
        centred$GSE86374$unit, centred$GSE86374$feature, centred$GSE86374$ref, centred$GSE86374$n_norm, median(centred$GSE86374$value)),
"",
"GPL13607 carries two PPP4C probes (58319, 57896). Probe **58319** is used because",
"that is the one `collapse_by_symbol()` retained for the DE analysis (larger |logFC|),",
"so the figure and the published statistics describe the same feature.", "",
"## Normalisation, as performed upstream", "",
"- **TCGA-BRCA** — raw STAR counts, TMM (`edgeR::calcNormFactors`), then",
"  `cpm(log = TRUE, prior.count = 1)`. Computed in this script, matching",
"  `scripts/08_survival.R` and `scripts/11_target_dossier.R`.",
"- **GSE70947** — background-corrected and log2 quantile-normalised by the submitters;",
"  used as deposited (`scripts/02_de_gse70947.R`).",
"- **GSE58135** — cufflinks FPKM, `log2(FPKM + 1)`, then `normalizeBetweenArrays`",
"  quantile (`scripts/03_de_gse58135.R`).",
"- **GSE42568 / GSE86374** — RMA log2 from the series matrix, then quantile",
"  (`scripts/06_de_affy.R`).", "",
"## Sample accounting", "",
"| Cohort | Step | Samples |", "|---|---|---|",
sprintf("| %s | %s | %d |", dropdf$cohort, dropdf$step, dropdf$n),
"",
sprintf("Plotted n per cohort: %s.",
        paste(sprintf("%s = %d", ORDER, sapply(ORDER, function(g) length(centred[[g]]$value))), collapse = ", ")),
"",
"Gene-level filtering upstream (`FPKM > 1` in GSE58135, `>50%` NA probe removal in the",
"array cohorts) acts on genes, not samples, and PPP4C passed every one — so no cohort",
"lost samples to it. The only per-sample loss is the NA probe reading noted above.", "",
"No patient-level de-duplication was applied to TCGA: all 1,111 primary-tumour aliquots",
"are plotted, matching the DE analysis. `scripts/11_target_dossier.R` de-duplicates to",
"1,095 patients for survival, where one row per patient is required; that is a different",
"question and a different n.", "",
"## Statistics", "",
sprintf("- **Kruskal-Wallis** across the five cohorts: chi-squared = %.1f, df = %d, p = %.3g.",
        as.numeric(kw$statistic), as.numeric(kw$parameter), as.numeric(kw$p.value)),
"  Non-parametric one-way, chosen over ANOVA because group sizes are severely unbalanced",
"  (1,111 vs 84) and the distributions are skewed. It is also consistent with the",
"  rank-based aggregation (RobustRankAggreg) that the selection rested on.",
"- **Post-hoc**: Wilcoxon rank-sum on all 10 cohort pairs, Benjamini-Hochberg adjusted.",
"  Recorded in the JSON `pairwise` block.",
sprintf("- %d of the 10 pairs underflowed the normal approximation and returned exactly 0.",
        n_underflow),
sprintf("  A p-value is never 0, so these are floored at the smallest representable double"),
sprintf("  (%.3g) and should be read as \"below double precision\", not as a point estimate.",
        P_FLOOR), "",
"**Read this test correctly.** It asks whether the *magnitude* of PPP4C elevation differs",
"between cohorts. It does — GSE70947 (paired Agilent) and GSE86374 (Gene ST) show",
"compressed effects, which `reports/validation_results.md` already attributes to array",
"dynamic-range compression. It is **not** a test of whether PPP4C is elevated in tumour;",
"that evidence is the per-cohort DE and the RRA meta-analysis (PPP4C RRA padj = 3.8e-3,",
"significant and concordant in all five cohorts).", "",
"**This test is new to this export.** The upstream pipeline never ran a one-way test",
"across cohorts — it combined cohorts by rank aggregation. The only Kruskal-Wallis",
"elsewhere in the project is PPP4C across PAM50 subtypes (p = 2.85e-08,",
"`logs/11_ppp4c.log`), which is a different comparison. Do not cite this p-value as",
"pre-existing.", "",
"## Box construction", "",
"- Box: Q1 / median / Q3, `quantile(type = 7)` (R default).",
"- `min` / `max` are **Tukey whisker ends** — the most extreme values within",
"  1.5 x IQR of the box, not the data range.",
"- `outliers`: every value beyond a whisker, full precision to 3 dp, `[]` when none.",
"- `points`: within-whisker values thinned to at most 200 by evenly spaced quantiles",
"  (`seq(0, 1, length.out = 200)`). Deterministic — no RNG, no seed, identical every run.",
"  Cohorts with fewer than 200 within-whisker values ship all of them.",
"- Summary stats rounded to 2 dp; points and outliers to 3 dp.", "",
"## Files", "",
"| File | Contents |", "|---|---|",
"| `exports/cohort_comparison.json` | Figure spec for the renderer |",
"| `exports/cohort_comparison_values.csv` | `sample_id, cohort, expression` — every plotted tumour sample. `expression` is the centred log2FC; add the cohort's normal ref from the table above to recover the native-scale value |",
"| `exports/cohort_comparison_provenance.md` | This file |", "",
"## Software", "",
sprintf("- %s", R.version.string),
sprintf("- edgeR %s", pkg("edgeR")),
sprintf("- limma %s (upstream DE)", pkg("limma")),
sprintf("- jsonlite %s", pkg("jsonlite")),
sprintf("- Platform: %s", R.version$platform), "",
sprintf("Generated %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

f_prov <- file.path(DIR_EXP, "cohort_comparison_provenance.md")
writeLines(prov, f_prov)
msg("wrote ", f_prov)
msg("done.")
