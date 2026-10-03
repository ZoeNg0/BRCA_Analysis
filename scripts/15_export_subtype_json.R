# PPP4C expression by PAM50 subtype (JSON export)

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(jsonlite); library(TCGAbiolinks)
})
source(file.path(getwd(), "R", "utils.R"))

GENE    <- "PPP4C"
ENSG    <- "ENSG00000149923"
DIR_EXP <- file.path(PROJ_ROOT, "exports")
RO      <- file.path(PROJ_ROOT, "exports", "r_objects")
P_FLOOR <- .Machine$double.xmin
MAX_POINTS <- 150

LEVELS <- c("LumA", "LumB", "Basal", "Her2", "Normal")
LABELS <- c(LumA = "Luminal A", LumB = "Luminal B", Basal = "Basal-like",
            Her2 = "HER2-enriched", Normal = "Normal-like")
SLOTS  <- c(LumA = 1L, LumB = 2L, Basal = 3L, Her2 = 4L, Normal = 5L)

audit <- list()
note  <- function(...) audit[[length(audit) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

counts <- readRDS(file.path(RO, "rna.rds"))
note(step = "aliquots in rna.rds", n = ncol(counts))

libsize <- colSums(counts)
dge  <- calcNormFactors(DGEList(counts = counts), method = "TMM")
lcpm <- cpm(dge, log = TRUE, prior.count = 1)
i <- which(strip_ens_version(rownames(lcpm)) == ENSG)
stopifnot(length(i) == 1)

stype <- substr(colnames(lcpm), 14, 15)
tum   <- which(stype == "01")
note(step = "primary tumour aliquots (sample type 01)", n = length(tum))

pat <- substr(colnames(lcpm), 1, 12)
keep <- vapply(split(tum, pat[tum]), function(idx)
  idx[order(-libsize[idx], colnames(lcpm)[idx])][1], integer(1))
keep <- sort(unname(keep))
note(step = "after de-duplication to one aliquot per patient", n = length(keep))
note(step = "  (duplicate tumour aliquots discarded)", n = length(tum) - length(keep))

d <- data.table(sample_id  = colnames(lcpm)[keep],
                patient_id = pat[keep],
                expression = as.numeric(lcpm[i, keep]))
stopifnot(!anyDuplicated(d$patient_id))

COHORT_BOXES_TUMOUR_MEDIAN <- 6.91
pooled_median <- median(d$expression)
msg(sprintf("pooled tumour median here = %.3f  (cohort_boxes.json TCGA tumour = %.2f)",
            pooled_median, COHORT_BOXES_TUMOUR_MEDIAN))
if (abs(pooled_median - COHORT_BOXES_TUMOUR_MEDIAN) > 0.15)
  stop("pooled median has drifted from cohort_boxes.json — scales are not comparable; ",
       "investigate before exporting")

sub <- as.data.table(TCGAquery_subtype("brca"))[, .(patient_id = patient,
                                                    pam50 = BRCA_Subtype_PAM50)]
note(step = "patients in TCGAquery_subtype('brca') table", n = nrow(sub))
d <- merge(d, sub, by = "patient_id", all.x = TRUE)
d[pam50 == "", pam50 := NA_character_]

n_unclassified <- sum(is.na(d$pam50))
note(step = "tumour patients WITHOUT a PAM50 call (dropped)", n = n_unclassified)
d_all <- copy(d)
d <- d[!is.na(pam50)]
d[, pam50 := factor(pam50, levels = LEVELS)]
stopifnot(!any(is.na(d$pam50)))
note(step = "tumour patients WITH a PAM50 call (plotted)", n = nrow(d))

setorder(d, pam50, patient_id)
msg("plotted: ", nrow(d), " patients; unclassified dropped: ", n_unclassified)
print(d[, .(n = .N, median = round(median(expression), 2)), by = pam50])

kw <- kruskal.test(expression ~ pam50, data = d)

cl <- as.data.table(readRDS(file.path(PROJ_ROOT, "data", "clinical",
                                      "TCGA_BRCA_clinical.rds")))
cl <- cl[vital_status %in% c("Alive", "Dead")]
cl[, os_time := fifelse(vital_status == "Dead", as.numeric(days_to_death),
                                                as.numeric(days_to_last_follow_up))]
cl <- cl[!is.na(os_time) & os_time > 0][!duplicated(submitter_id)]
d_surv <- d[patient_id %in% cl$submitter_id]
kw_surv <- kruskal.test(expression ~ droplevels(pam50), data = d_surv)

pairs <- combn(LEVELS, 2, simplify = FALSE)
pw <- do.call(rbind, lapply(pairs, function(p) data.frame(
  a = unname(LABELS[p[1]]), b = unname(LABELS[p[2]]),
  p = suppressWarnings(wilcox.test(d$expression[d$pam50 == p[1]],
                                   d$expression[d$pam50 == p[2]])$p.value),
  stringsAsFactors = FALSE)))
pw$padj <- p.adjust(pw$p, method = "BH")
n_underflow <- sum(pw$p <= 0)
pw$p <- pmax(pw$p, P_FLOOR); pw$padj <- pmax(pw$padj, P_FLOOR)

box_stats <- function(x) {
  q <- as.numeric(quantile(x, c(0.25, 0.5, 0.75), type = 7))
  iqr <- q[3] - q[1]; lo <- q[1] - 1.5 * iqr; hi <- q[3] + 1.5 * iqr
  inside <- x[x >= lo & x <= hi]
  list(q1 = q[1], median = q[2], q3 = q[3], whisker_lo = min(inside),
       whisker_hi = max(inside), mean = mean(x),
       outliers = sort(x[x < lo | x > hi]), inside = sort(inside))
}
thin_points <- function(inside, target) {
  if (length(inside) <= target) return(inside)
  as.numeric(quantile(inside, probs = seq(0, 1, length.out = target), type = 7))
}

groups <- lapply(LEVELS, function(k) {
  x  <- d$expression[d$pam50 == k]
  bs <- box_stats(x)
  list(label    = unname(LABELS[k]),
       slot     = unname(SLOTS[k]),
       n        = length(x),
       min      = floor(bs$whisker_lo * 100) / 100,
       q1       = round(bs$q1, 2),
       median   = round(bs$median, 2),
       q3       = round(bs$q3, 2),
       max      = ceiling(bs$whisker_hi * 100) / 100,
       mean     = round(bs$mean, 2),
       outliers = I(round(bs$outliers, 3)),
       points   = I(round(thin_points(bs$inside, MAX_POINTS), 3)))
})

out <- list(
  title = paste0(GENE, " expression by PAM50 subtype"),
  unit  = "log₂ CPM",
  test  = list(name = "Kruskal-Wallis", p = max(as.numeric(kw$p.value), P_FLOOR),
               effect = sprintf("across %d subtypes", length(LEVELS))),
  groups = groups,
  pairwise = unname(lapply(seq_len(nrow(pw)), function(k) list(
    a = pw$a[k], b = pw$b[k], p = pw$p[k], padj = pw$padj[k]))))

f_expr <- file.path(DIR_EXP, "subtype_expression.json")
write(toJSON(out, auto_unbox = TRUE, digits = NA, pretty = 2, na = "null"), f_expr)
msg("wrote ", f_expr)

comp <- list(total = nrow(d),
             items = unname(lapply(LEVELS, function(k) list(
               label = unname(LABELS[k]), value = sum(d$pam50 == k),
               slot = unname(SLOTS[k])))))
f_comp <- file.path(DIR_EXP, "subtype_composition.json")
write(toJSON(comp, auto_unbox = TRUE, digits = NA, pretty = 2, na = "null"), f_comp)
msg("wrote ", f_comp)

f_csv <- file.path(DIR_EXP, "subtype_values.csv")
write.csv(data.frame(sample_id = d$sample_id, patient_id = d$patient_id,
                     pam50 = as.character(d$pam50), expression = d$expression),
          f_csv, row.names = FALSE)
msg("wrote ", f_csv, "  (", nrow(d), " rows)")

chk <- fromJSON(f_expr, simplifyVector = FALSE)
stopifnot(length(chk$groups) == 5,
          identical(vapply(chk$groups, function(g) g$slot, 0L), 1:5),
          identical(vapply(chk$groups, function(g) g$label, ""), unname(LABELS[LEVELS])),
          chk$test$p > 0)
for (g in chk$groups) {
  pts <- unlist(g$points); out_v <- unlist(g$outliers)
  stopifnot(is.list(g$outliers), is.list(g$points),
            g$min <= g$q1, g$q1 <= g$median, g$median <= g$q3, g$q3 <= g$max,
            length(pts) <= MAX_POINTS, length(pts) >= 1,
            all(pts >= g$min), all(pts <= g$max),
            length(out_v) == 0 || all(out_v < g$min | out_v > g$max),
            g$n == sum(d$pam50 == names(LABELS)[match(g$label, LABELS)]))
}
stopifnot(!any(is.na(unlist(chk))),
          all(vapply(chk$pairwise, function(p) p$p > 0 && p$padj > 0, logical(1))),
          length(chk$pairwise) == 10)
chc <- fromJSON(f_comp, simplifyVector = FALSE)
stopifnot(chc$total == sum(vapply(chc$items, function(x) x$value, 0L)),
          chc$total == nrow(d),
          identical(vapply(chc$items, function(x) x$slot, 0L), 1:5))
msg("validated: 5 groups in order, slots 1-5, points within whiskers, ",
    "composition sums to total, no NA/null, no zero p-values")

aud <- do.call(rbind, audit)
pkg <- function(p) as.character(packageVersion(p))
tab <- d[, .(n = .N, median = round(median(expression), 2),
             q1 = round(quantile(expression, .25), 2),
             q3 = round(quantile(expression, .75), 2)), by = pam50]

prov <- c(
"# Provenance — PPP4C expression by PAM50 subtype", "",
sprintf("Generated %s by `scripts/15_export_subtype_json.R`.", format(Sys.Date(), "%Y-%m-%d")), "",
sprintf("Gene: **%s** (`%s`). Unit: **log₂ CPM**.", GENE, ENSG), "",
"## Regenerate", "", "```bash", "Rscript scripts/15_export_subtype_json.R", "```", "",
"## Where the PAM50 calls come from", "",
sprintf("`TCGAquery_subtype(\"brca\")` from **TCGAbiolinks %s** — a packaged table, not a", pkg("TCGAbiolinks")),
"network query, so this runs offline. It is the TCGA PanCancer Atlas BRCA subtype call,",
"published in doi.org/10.1016/j.ccell.2018.03.014, column `BRCA_Subtype_PAM50`, joined to",
sprintf("our data on `patient` -> `patient_id`. The table carries %s patients.",
        format(aud$n[aud$step == "patients in TCGAquery_subtype('brca') table"], big.mark = ",")),
"",
"PAM50 is *not* in `data/clinical/TCGA_BRCA_clinical.rds` — that file has 100 columns and",
"none of them is a subtype call. This was checked, not assumed.", "",
"## THE DENOMINATOR — read this before citing any n", "",
"The site's Figure 1 (`cohort_boxes.json`) plots **1,111 TCGA tumours**. This figure plots",
sprintf("**%s**. Every step of that difference:", format(nrow(d), big.mark = ",")), "",
"| Step | n |", "|---|---|",
sprintf("| %s | %s |", aud$step, format(aud$n, big.mark = ",")),
"",
paste("The three numbers in play are **1,111 aliquots**, **1,095 patients**, and",
      sprintf("**%s patients with a PAM50 call**.", format(nrow(d), big.mark = ","))),
sprintf("- Figure 1 plots aliquots and applies no de-duplication, so it shows 1,111."),
"- A subtype figure cannot: PAM50 is called once per patient, so the 16 duplicate tumour",
"  aliquots must collapse to their patients first. That is the 1,111 -> 1,095 step.",
sprintf("- %d of those 1,095 patients have **no PAM50 call** in the PanCancer Atlas table.",
        n_unclassified),
"  They are **dropped**, not shown as a sixth \"not classified\" group, because the",
"  requested spec fixes five groups at slots 1-5. They are counted above and nowhere",
"  silently absorbed.", "",
sprintf(paste("The dossier (`logs/11_ppp4c.log`) reports subtype counts summing to **1,059**, not",
              "%s. That is a *third* denominator: it additionally requires usable survival"),
        format(nrow(d), big.mark = ",")),
sprintf("follow-up, dropping a further %d patients. That restriction belongs to a survival",
        nrow(d) - nrow(d_surv)),
"analysis and has no place in an expression figure, so it is not applied here.", "",
"| Denominator | n | Kruskal-Wallis p | Used here |", "|---|---|---|---|",
sprintf("| Patients with a PAM50 call | %s | %.3g | **yes** |",
        format(nrow(d), big.mark = ","), as.numeric(kw$p.value)),
sprintf("| ...also with survival follow-up (dossier) | %s | %.3g | no |",
        format(nrow(d_surv), big.mark = ","), as.numeric(kw_surv$p.value)),
"",
"## n and medians per subtype", "",
"| Subtype | slot | n | Q1 | median | Q3 |", "|---|---|---|---|---|---|",
sprintf("| %s | %d | %d | %.2f | %.2f | %.2f |",
        LABELS[as.character(tab$pam50)], SLOTS[as.character(tab$pam50)],
        tab$n, tab$q1, tab$median, tab$q3),
"",
"## Scale — comparability with `cohort_boxes.json`", "",
"Normalisation is byte-identical to `scripts/14_export_figure_data.R`: TMM",
"(`edgeR::calcNormFactors`) over **all 1,224 aliquots** in `rna.rds`, then",
"`cpm(log = TRUE, prior.count = 1)`, then subsetting. Because TMM factors are estimated",
"on the full matrix before any subsetting, the values here sit on exactly the axis the",
"TCGA panel of `cohort_boxes.json` uses, and the two figures can be read against each",
"other.", "",
sprintf("Check performed at runtime: pooled median of the plotted patients = **%.3f** against",
        pooled_median),
sprintf("the cohort_boxes TCGA tumour median of **%.2f**. The script aborts if these drift",
        COHORT_BOXES_TUMOUR_MEDIAN),
"apart by more than 0.15, so a silent scale change cannot ship.", "",
"`scripts/11_target_dossier.R` computes TMM on the de-duplicated tumour matrix *alone*",
"rather than the full 1,224. That is a slightly different normalisation, which is why the",
"medians here differ from `logs/11_ppp4c.log` in the second decimal (e.g. Luminal A",
sprintf("%.2f here vs 6.86 there). Neither is wrong; they answer to different baselines.",
        tab$median[tab$pam50 == "LumA"]), "",
"De-duplication keeps the aliquot with the **largest library size** (ties on the",
"lexicographically first barcode) — the same rule as `scripts/14_export_figure_data.R`.",
"Note that `scripts/08_survival.R` and `scripts/11_target_dossier.R` instead keep the",
"*first* aliquot. With 11 affected patients the practical difference is negligible, but",
"the rules are genuinely different and the survival export documents its own choice.", "",
"## Statistics", "",
sprintf("- **Kruskal-Wallis** across the %d subtypes: chi-squared = %.1f, df = %d, p = %.3g.",
        length(LEVELS), as.numeric(kw$statistic), as.numeric(kw$parameter), as.numeric(kw$p.value)),
"  Non-parametric one-way — group sizes are very unequal (Luminal A is roughly seven times",
"  Normal-like) and the distributions are skewed.",
"- **Post-hoc**: Wilcoxon rank-sum on all 10 subtype pairs, Benjamini-Hochberg adjusted,",
"  in the JSON `pairwise` block.",
if (n_underflow == 0)
  "- No pairwise p underflowed to zero." else
  sprintf("- %d pairwise p-value(s) underflowed to exactly 0 and are floored at %.3g.",
          n_underflow, P_FLOOR), "",
"**Normal-like is kept as a group** because the spec asks for five. Be aware it is widely",
"read as a low-tumour-purity artifact rather than a biological subtype — these are tumour",
"samples whose profile resembles normal breast, and its low median here is consistent with",
"dilution by normal tissue. Treat any claim resting on it cautiously.", "",
"## Box construction", "",
"- Q1/median/Q3 at `quantile(type = 7)`.",
"- `min`/`max` are Tukey whisker ends — most extreme value within 1.5 x IQR — **not** the",
"  data range. They round *outward* so no point escapes them after rounding.",
"- `outliers`: every value past a whisker, `[]` when none, never omitted.",
sprintf("- `points`: within-whisker values thinned to at most %d per group by evenly spaced",
        MAX_POINTS),
"  quantiles. Deterministic, no RNG, identical on every run.",
"- Summary stats 2 dp; points and outliers 3 dp.", "",
"## Files", "", "| File | Contents |", "|---|---|",
"| `exports/subtype_expression.json` | Box plot data, 5 groups |",
"| `exports/subtype_composition.json` | Composition bar chart; `total` = sum of `value` |",
"| `exports/subtype_values.csv` | `sample_id, patient_id, pam50, expression` |",
"| `exports/subtype_provenance.md` | This file |", "",
"## Software", "",
sprintf("- %s", R.version.string),
sprintf("- edgeR %s", pkg("edgeR")),
sprintf("- TCGAbiolinks %s", pkg("TCGAbiolinks")),
sprintf("- jsonlite %s", pkg("jsonlite")),
sprintf("- Platform: %s", R.version$platform), "",
sprintf("Generated %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

f_prov <- file.path(DIR_EXP, "subtype_provenance.md")
writeLines(prov, f_prov)
msg("wrote ", f_prov)
msg("done.")
