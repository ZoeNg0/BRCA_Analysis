# PPP4C high vs low 5-year OS (JSON export)

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(survival); library(jsonlite)
  library(TCGAbiolinks); library(AnnotationDbi); library(org.Hs.eg.db)
})
source(file.path(getwd(), "R", "utils.R"))

GENE     <- "PPP4C"
ENSG     <- "ENSG00000149923"
DIR_EXP  <- file.path(PROJ_ROOT, "exports")
RO       <- file.path(PROJ_ROOT, "exports", "r_objects")
MAX_DAYS <- 5 * 365.25
MAX_TIME <- 60
P_FLOOR  <- .Machine$double.xmin
SLOT     <- c(Low = 1L, High = 8L)

audit <- list()
note  <- function(...) audit[[length(audit) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

counts <- readRDS(file.path(RO, "rna.rds"))
note(step = "aliquots in rna.rds", n = ncol(counts))
stype  <- substr(colnames(counts), 14, 15)
counts <- counts[, stype == "01", drop = FALSE]
note(step = "primary tumour aliquots (sample type 01)", n = ncol(counts))

patient <- substr(colnames(counts), 1, 12)
n_dup   <- sum(duplicated(patient))
counts  <- counts[, !duplicated(patient), drop = FALSE]
colnames(counts) <- patient[!duplicated(patient)]
note(step = "  duplicate aliquots dropped (first kept)", n = n_dup)
note(step = "unique patients with expression", n = ncol(counts))

dge  <- calcNormFactors(DGEList(counts = counts), method = "TMM")
lcpm <- cpm(dge, log = TRUE, prior.count = 1)
rownames(lcpm) <- strip_ens_version(rownames(lcpm))
i <- which(rownames(lcpm) == ENSG)
stopifnot(length(i) == 1)

cl <- as.data.table(readRDS(file.path(PROJ_ROOT, "data", "clinical",
                                      "TCGA_BRCA_clinical.rds")))
note(step = "clinical records on file", n = nrow(cl))
n0 <- nrow(cl)
cl <- cl[!is.na(vital_status) & vital_status %in% c("Alive", "Dead")]
note(step = "  dropped: vital_status not Alive/Dead", n = n0 - nrow(cl))

cl[, os_time  := fifelse(vital_status == "Dead", as.numeric(days_to_death),
                                                 as.numeric(days_to_last_follow_up))]
cl[, os_event := as.integer(vital_status == "Dead")]
n1 <- nrow(cl); cl <- cl[!is.na(os_time) & os_time > 0]
note(step = "  dropped: missing or non-positive follow-up time", n = n1 - nrow(cl))
n2 <- nrow(cl); cl <- cl[!duplicated(submitter_id)]
note(step = "  dropped: duplicate clinical rows", n = n2 - nrow(cl))
note(step = "clinical with usable follow-up", n = nrow(cl))

n_beyond5 <- sum(cl$os_time > MAX_DAYS & cl$os_event == 1L)
cl[os_time > MAX_DAYS, os_event := 0L]
cl[, os_time := pmin(os_time, MAX_DAYS)]
cl[, os_months := os_time / 30.4375]
note(step = "  deaths after 5 yr, administratively censored", n = n_beyond5)

shared <- intersect(cl$submitter_id, colnames(lcpm))
note(step = "  clinical rows with no expression match", n = nrow(cl) - length(shared))
cl <- cl[submitter_id %in% shared]
cl[, expr := as.numeric(lcpm[i, submitter_id])]
note(step = "ANALYSIS SET (expression + clinical)", n = nrow(cl))
note(step = "  deaths within 5 years", n = sum(cl$os_event))

cut <- median(cl$expr)
cl[, group := factor(ifelse(expr > cut, "High", "Low"), levels = c("Low", "High"))]
n_tie <- sum(cl$expr == cut)
note(step = "  patients exactly at the median (assigned Low)", n = n_tie)
msg(sprintf("median split at %.4f log2 CPM -> Low %d / High %d (%d exactly at median)",
            cut, sum(cl$group == "Low"), sum(cl$group == "High"), n_tie))

fit <- survfit(Surv(os_months, os_event) ~ group, data = cl)
lr  <- survdiff(Surv(os_months, os_event) ~ group, data = cl)
lr_p <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)

cox_uni <- coxph(Surv(os_months, os_event) ~ group, data = cl)
cu <- summary(cox_uni)
zph <- cox.zph(cox_uni)
zph_p <- zph$table["GLOBAL", "p"]

cl[, age_years := as.numeric(age_at_diagnosis) / 365.25]
cl[, stage := factor(fcase(
  grepl("Stage I($|[AB])",  ajcc_pathologic_stage), "I",
  grepl("Stage II($|[AB])", ajcc_pathologic_stage), "II",
  grepl("Stage III",        ajcc_pathologic_stage), "III",
  grepl("Stage IV",         ajcc_pathologic_stage), "IV",
  default = NA_character_), levels = c("I", "II", "III", "IV"))]
sub <- as.data.table(TCGAquery_subtype("brca"))[, .(submitter_id = patient,
                                                    pam50 = BRCA_Subtype_PAM50)]
cl <- merge(cl, sub, by = "submitter_id", all.x = TRUE)
cl[, pam50 := factor(pam50)]
cc <- cl[!is.na(age_years) & !is.na(stage) & !is.na(pam50)]
cox_adj <- coxph(Surv(os_months, os_event) ~ group + age_years + stage + pam50, data = cc)
ca <- summary(cox_adj)
note(step = "  complete cases for the adjusted model", n = nrow(cc))

strata_names <- sub("^group=", "", names(fit$strata))
idx <- rep(seq_along(fit$strata), fit$strata)
n_na_ci <- 0L

curves <- lapply(seq_along(strata_names), function(k) {
  g  <- strata_names[k]
  s  <- which(idx == k)
  dd <- cl[group == g]

  tt <- fit$time[s]; ss <- fit$surv[s]; nr <- fit$n.risk[s]
  lo <- fit$lower[s]; hi <- fit$upper[s]; nc <- fit$n.censor[s]

  na_here <- sum(is.na(lo) | is.na(hi))
  n_na_ci <<- n_na_ci + na_here
  lo[is.na(lo)] <- ss[is.na(lo)]; hi[is.na(hi)] <- ss[is.na(hi)]
  lo <- pmax(pmin(lo, 1), 0); hi <- pmax(pmin(hi, 1), 0)

  pts <- data.table(t = tt, s = ss, lo = lo, hi = hi,
                    atRisk = as.integer(nr), censored = nc > 0)
  pts <- rbind(data.table(t = 0, s = 1, lo = 1, hi = 1,
                          atRisk = nrow(dd), censored = FALSE), pts)
  last <- pts[.N]
  if (last$t < MAX_TIME)
    pts <- rbind(pts, data.table(t = MAX_TIME, s = last$s, lo = last$lo, hi = last$hi,
                                 atRisk = sum(dd$os_months >= MAX_TIME), censored = FALSE))
  pts <- pts[t <= MAX_TIME]

  list(name   = paste0(GENE, "-", tolower(g)),
       n      = nrow(dd),
       slot   = unname(SLOT[g]),
       events = sum(dd$os_event),
       points = unname(lapply(seq_len(nrow(pts)), function(j) list(
         t        = round(pts$t[j], 2),
         s        = round(pts$s[j], 4),
         lo       = round(pts$lo[j], 4),
         hi       = round(pts$hi[j], 4),
         atRisk   = pts$atRisk[j],
         censored = pts$censored[j]))))
})
names(curves) <- strata_names
curves <- curves[c("Low", "High")]

for (g in c("Low", "High"))
  msg(sprintf("%-5s n=%4d events=%3d  curve points=%d",
              g, curves[[g]]$n, curves[[g]]$events, length(curves[[g]]$points)))

out <- list(
  gene        = GENE,
  endpoint    = "Overall survival",
  split       = "median expression",
  maxTime     = MAX_TIME,
  timeUnit    = "months",
  pValue      = max(as.numeric(lr_p), P_FLOOR),
  hazardRatio = round(as.numeric(cu$conf.int[1, "exp(coef)"]), 2),
  hrCI        = I(round(as.numeric(cu$conf.int[1, c("lower .95", "upper .95")]), 2)),
  hazardRatioAdjusted = round(as.numeric(ca$conf.int[1, "exp(coef)"]), 2),
  hrCIAdjusted        = I(round(as.numeric(ca$conf.int[1, c("lower .95", "upper .95")]), 2)),
  adjustedFor         = "age + AJCC stage + PAM50 subtype",
  adjustedN           = nrow(cc),
  adjustedP           = signif(as.numeric(ca$coefficients[1, "Pr(>|z|)"]), 3),
  events      = list(high = curves$High$events, low = curves$Low$events),
  groups      = unname(curves))

f_json <- file.path(DIR_EXP, "survival_5yr.json")
write(toJSON(out, auto_unbox = TRUE, digits = NA, pretty = 2, na = "null"), f_json)
msg("wrote ", f_json)

f_csv <- file.path(DIR_EXP, "survival_5yr_values.csv")
write.csv(data.frame(patient_id = cl$submitter_id, expression = cl$expr,
                     group = as.character(cl$group), os_months = cl$os_months,
                     os_event = cl$os_event), f_csv, row.names = FALSE)
msg("wrote ", f_csv, "  (", nrow(cl), " rows)")

chk <- fromJSON(f_json, simplifyVector = FALSE)
stopifnot(length(chk$groups) == 2,
          chk$groups[[1]]$name == "PPP4C-low",  chk$groups[[1]]$slot == 1,
          chk$groups[[2]]$name == "PPP4C-high", chk$groups[[2]]$slot == 8,
          chk$maxTime == MAX_TIME, chk$pValue > 0)
for (g in chk$groups) {
  tv <- vapply(g$points, function(p) p$t, 0)
  sv <- vapply(g$points, function(p) p$s, 0)
  lv <- vapply(g$points, function(p) p$lo, 0)
  hv <- vapply(g$points, function(p) p$hi, 0)
  av <- vapply(g$points, function(p) p$atRisk, 0L)
  stopifnot(all(tv >= 0), all(tv <= MAX_TIME), !is.unsorted(tv),
            all(diff(sv) <= 1e-9),
            all(sv >= 0), all(sv <= 1),
            all(lv <= sv + 1e-9), all(sv <= hv + 1e-9),
            all(lv >= 0), all(hv <= 1),
            all(diff(av) <= 0),
            g$points[[1]]$t == 0, g$points[[1]]$s == 1,
            tv[length(tv)] == MAX_TIME)
}
stopifnot(chk$groups[[1]]$n + chk$groups[[2]]$n == nrow(cl),
          chk$events$low + chk$events$high == sum(cl$os_event),
          !any(is.na(unlist(chk))))

pub <- fread(file.path(PROJ_ROOT, "exports", "validation", "survival_5yr_summary.csv"))
pub <- pub[gene == GENE]
stopifnot(nrow(pub) == 1,
          pub$n == nrow(cl), pub$events == sum(cl$os_event),
          abs(pub$HR - as.numeric(cu$conf.int[1, "exp(coef)"])) < 1e-6,
          abs(pub$logrank_p - lr_p) < 1e-9)
msg("validated: curves monotone, CI brackets s, at-risk non-increasing, ",
    "t in [0,60], n/events/HR/log-rank match survival_5yr_summary.csv")

aud <- do.call(rbind, audit)
pkg <- function(p) as.character(packageVersion(p))
tp  <- summary(fit, times = c(12, 24, 36, 48, 60))
tpd <- data.table(stratum = sub("^group=", "", as.character(tp$strata)),
                  t = tp$time, s = tp$surv, nrisk = tp$n.risk)

prov <- c(
"# Provenance — PPP4C 5-year survival curves", "",
sprintf("Generated %s by `scripts/16_export_survival_json.R`.", format(Sys.Date(), "%Y-%m-%d")), "",
"## Regenerate", "", "```bash", "Rscript scripts/16_export_survival_json.R", "```", "",
"## What this is", "",
"The same analysis `scripts/08_survival.R` already plots to",
"`plots/survival/KM_PPP4C.png`, exported as data so the site can draw it interactively.",
"The script **aborts** unless the recomputed n, events, hazard ratio and log-rank p match",
"`exports/validation/survival_5yr_summary.csv`, so the page and the PNG cannot drift.", "",
"## Which n this figure uses — 1,111 vs 1,095 vs 1,071", "",
"The distinction you flagged does apply here, and this figure uses the **smallest** of the",
"three. Full accounting:", "",
"| Step | n |", "|---|---|",
sprintf("| %s | %s |", aud$step, format(aud$n, big.mark = ",")),
"",
"- Figure 1 (`cohort_boxes.json`) plots **1,111 aliquots** with no de-duplication.",
sprintf("- Survival needs one row per patient, so 16 duplicate aliquots collapse -> **1,095**."),
sprintf("- Of those, %s have usable follow-up and expression -> **%s patients, %d deaths**",
        format(nrow(cl), big.mark = ","), format(nrow(cl), big.mark = ","), sum(cl$os_event)),
"  within 5 years. That is the number on this chart.", "",
"**The page should say 1,111 for Figure 1 and",
sprintf("%s for the survival curve.** They are not the same cohort and there is no", format(nrow(cl), big.mark = ",")),
"single n that honestly covers both.", "",
sprintf("De-duplication here keeps the **first** aliquot per patient, replicating"),
"`scripts/08_survival.R`. Note `scripts/14`/`15` keep the aliquot with the largest library",
"size instead. Eleven patients are affected and the two rules pick the same aliquot for",
"most of them, but they are genuinely different rules and the choice is made explicitly in",
"each script rather than assumed.", "",
"Likewise TMM is computed here on the **de-duplicated tumour matrix alone** (1,095",
"columns), as `08_survival.R` does — not on all 1,224 aliquots as `scripts/14`/`15` do.",
"That is why a patient's log2 CPM here can differ in the second decimal from",
"`subtype_values.csv`. Within this figure everything is on one internally consistent scale.", "",
"## Median split", "",
sprintf("- Cut-point: **%.4f log₂ CPM**, the median of the %s analysis-set patients.",
        cut, format(nrow(cl), big.mark = ",")),
sprintf("- Low **n = %d** (%d deaths), High **n = %d** (%d deaths).",
        curves$Low$n, curves$Low$events, curves$High$n, curves$High$events),
sprintf("- The rule is `expr > median`, so ties at the median fall in **Low**. %s",
        if (n_tie == 0) "No patient sat exactly on the median." else
          sprintf("%d patient(s) sat exactly on the median.", n_tie)),
"- Pre-specified, **not** an optimised cut-point. Choosing the cut that minimises p inflates",
"  significance badly and would need permutation correction; the median needs none.", "",
"## Censoring", "",
sprintf("- Administrative censoring at 5 years: `os_time > %.2f days` is truncated and the",
        MAX_DAYS),
"  event flag cleared.",
sprintf("- %d death(s) occurred after 5 years and are censored at 60 months.", n_beyond5),
"- `os_months = days / 30.4375`, so the horizon is exactly 60.00 months and every `t`",
"  in the JSON satisfies 0 <= t <= 60.", "",
"## Curve data", "",
"- `points` carries **every distinct time survfit returns** — deaths and censorings alike.",
"  Death times are the staircase; censoring-only times are where the ticks go, flagged by",
"  `censored: true`. Several censorings at one time collapse to a single point.",
sprintf("- Point counts: Low **%d**, High **%d**. Both are under the ~400 you set as the",
        length(curves$Low$points), length(curves$High$points)),
"  threshold for reconsidering, so nothing is thinned — these are all the steps.",
"- `t = 0, s = 1.0` is prepended; a final point at `t = 60` carries the last estimate out to",
"  the axis end so the curve does not stop short.",
sprintf("- `lo`/`hi` are survfit's real Greenwood 95%% bounds (`$lower`/`$upper`), not a"),
"  synthetic band.",
if (n_na_ci == 0)
  "- **No bound came back NA**, so nothing was substituted." else
  sprintf(paste("- %d point(s) had an NA bound — survfit returns NA while the estimate is still",
                "exactly 1 and the variance is 0. Those were set to `s` (a zero-width interval,",
                "which is what survfit means there). This is the only substitution made and it",
                "affects no point after the first death."), n_na_ci),
"- `atRisk` is survfit's `n.risk`, for the risk table.", "",
"## Survival at fixed times", "",
"| Stratum | 12 mo | 24 mo | 36 mo | 48 mo | 60 mo |", "|---|---|---|---|---|---|",
sprintf("| %s | %s |", c("Low", "High"),
        vapply(c("Low", "High"), function(g)
          paste(sprintf("%.3f", tpd[stratum == g][order(t)]$s), collapse = " | "), "")),
"",
"## Statistics", "",
sprintf("- **Log-rank** within the 5-year window: chi-squared = %.2f, p = **%.4g**.",
        as.numeric(lr$chisq), as.numeric(lr_p)),
sprintf("- **Cox, univariate** (median split): HR **%.2f** (95%% CI %.2f-%.2f), p = %.3g.",
        cu$conf.int[1, "exp(coef)"], cu$conf.int[1, "lower .95"],
        cu$conf.int[1, "upper .95"], cu$coefficients[1, "Pr(>|z|)"]),
sprintf("- **Cox, adjusted** for age + stage + PAM50 (n = %d): HR **%.2f** (95%% CI %.2f-%.2f), p = %.3g.",
        nrow(cc), ca$conf.int[1, "exp(coef)"], ca$conf.int[1, "lower .95"],
        ca$conf.int[1, "upper .95"], ca$coefficients[1, "Pr(>|z|)"]), "",
"### Which HR to show", "",
sprintf(paste("**Show the univariate HR (%.2f) as the headline, because it is the one the curve",
              "itself depicts** — a median split with no covariates. But the adjusted HR"),
        cu$conf.int[1, "exp(coef)"]),
sprintf(paste("(%.2f, p = %.3g) must appear alongside it, because it is the number that decides",
              "whether this is a real independent signal, and it sits on the wrong side of"),
        ca$conf.int[1, "exp(coef)"], ca$coefficients[1, "Pr(>|z|)"]),
"0.05. Showing only the univariate value would overstate the result. Both are in the JSON",
"(`hazardRatio`/`hrCI` and `hazardRatioAdjusted`/`hrCIAdjusted`).", "",
"### Proportional hazards", "",
sprintf("`cox.zph` on the univariate 5-year model: global p = **%.3g**.", as.numeric(zph_p)),
if (as.numeric(zph_p) < 0.05)
  paste("That is a **violation** — the hazard ratio is not constant over the 5 years, so the",
        "single HR is a time-averaged summary rather than a constant multiplier. Since the",
        "page displays a hazard ratio, say so. The curves themselves remain valid; it is the",
        "single-number summary that is weaker than it looks.") else
  paste("No evidence against proportional hazards over this window, so a single HR is a fair",
        "summary."),
sprintf("The dossier reports global p = 0.016 for the fully *adjusted* model (`logs/11_ppp4c.log`),"),
"a separate check on a different model; both point the same way.", "",
"## Files", "", "| File | Contents |", "|---|---|",
"| `exports/survival_5yr.json` | Curve data for the interactive chart |",
"| `exports/survival_5yr_values.csv` | `patient_id, expression, group, os_months, os_event` |",
"| `exports/survival_provenance.md` | This file |", "",
"## Software", "",
sprintf("- %s", R.version.string),
sprintf("- survival %s", pkg("survival")),
sprintf("- edgeR %s", pkg("edgeR")),
sprintf("- TCGAbiolinks %s", pkg("TCGAbiolinks")),
sprintf("- jsonlite %s", pkg("jsonlite")),
sprintf("- Platform: %s", R.version$platform), "",
sprintf("Generated %s.", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

writeLines(prov, file.path(DIR_EXP, "survival_provenance.md"))
msg("wrote ", file.path(DIR_EXP, "survival_provenance.md"))

msg("\n=== curve summary ===")
print(as.data.frame(dcast(tpd, stratum ~ t, value.var = "s")))
msg("log-rank p = ", signif(lr_p, 4))
msg("done.")
