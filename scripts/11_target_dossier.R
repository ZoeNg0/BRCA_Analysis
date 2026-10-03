# Per-gene summary, e.g. GENE=PPP4C Rscript scripts/11_target_dossier.R

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(survival); library(survminer)
  library(TCGAbiolinks); library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2)
})
source(file.path(getwd(), "R", "utils.R"))

GENE     <- Sys.getenv("GENE", "PPP4C")
ann1 <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = GENE, keytype = "SYMBOL", columns = c("ENSEMBL", "MAP")))
ENSG     <- ann1$ENSEMBL[1]
LOCUS    <- ann1$MAP[1]
CHR      <- sub("[pq].*", "", LOCUS)
MAX_DAYS <- 5 * 365.25
DIR_D    <- file.path(DIR_PLOTS, tolower(GENE))
dir.create(DIR_D, recursive = TRUE, showWarnings = FALSE)
PAL <- c(Low = "#2a78d6", High = "#eb6834")

counts <- readRDS(file.path(PROJ_ROOT, "exports", "r_objects", "rna.rds"))
counts <- counts[, substr(colnames(counts), 14, 15) == "01", drop = FALSE]
counts <- counts[, !duplicated(substr(colnames(counts), 1, 12)), drop = FALSE]
colnames(counts) <- substr(colnames(counts), 1, 12)
dge  <- calcNormFactors(DGEList(counts = counts), method = "TMM")
lcpm <- cpm(dge, log = TRUE, prior.count = 1)
rownames(lcpm) <- strip_ens_version(rownames(lcpm))
lcpm <- lcpm[rowMeans(lcpm) > 0, , drop = FALSE]
msg("tumours: ", ncol(lcpm), " patients")

cl <- as.data.table(readRDS(file.path(PROJ_ROOT, "data", "clinical", "TCGA_BRCA_clinical.rds")))
cl <- cl[vital_status %in% c("Alive", "Dead")]
cl[, os_time  := fifelse(vital_status == "Dead", as.numeric(days_to_death),
                                                 as.numeric(days_to_last_follow_up))]
cl[, os_event := as.integer(vital_status == "Dead")]
cl <- cl[!is.na(os_time) & os_time > 0][!duplicated(submitter_id)]
cl[os_time > MAX_DAYS, os_event := 0L]
cl[, os_time := pmin(os_time, MAX_DAYS)]
cl[, os_months := os_time / 30.4375]
cl[, age_years := as.numeric(age_at_diagnosis) / 365.25]

cl[, stage := fcase(
  grepl("Stage I($|[AB])",   ajcc_pathologic_stage), "I",
  grepl("Stage II($|[AB])",  ajcc_pathologic_stage), "II",
  grepl("Stage III",         ajcc_pathologic_stage), "III",
  grepl("Stage IV",          ajcc_pathologic_stage), "IV",
  default = NA_character_)]
cl[, stage := factor(stage, levels = c("I", "II", "III", "IV"))]

sub <- as.data.table(TCGAquery_subtype("brca"))[, .(submitter_id = patient,
                                                    pam50 = BRCA_Subtype_PAM50)]
cl <- merge(cl, sub, by = "submitter_id", all.x = TRUE)
cl[, pam50 := factor(pam50, levels = c("LumA", "LumB", "Her2", "Basal", "Normal"))]

shared <- intersect(cl$submitter_id, colnames(lcpm))
cl <- cl[submitter_id %in% shared]
cl[, expr := as.numeric(lcpm[ENSG, submitter_id])]
cl[, group := factor(ifelse(expr > median(expr), "High", "Low"), levels = c("Low", "High"))]
cl[, expr_z := as.numeric(scale(expr))]
msg("analysis set: ", nrow(cl), " patients, ", sum(cl$os_event), " deaths in 5 yr")

fit_report <- function(f, label, data) {
  m  <- coxph(f, data = data)
  cs <- summary(m)
  i  <- grep("^group|^expr_z", rownames(cs$conf.int))[1]
  data.table(model = label, n = cs$n, events = cs$nevent,
             term = rownames(cs$conf.int)[i],
             HR = cs$conf.int[i, "exp(coef)"],
             lo = cs$conf.int[i, "lower .95"],
             hi = cs$conf.int[i, "upper .95"],
             p  = cs$coefficients[i, "Pr(>|z|)"])
}

cc <- cl[!is.na(stage) & !is.na(pam50) & !is.na(age_years)]
msg("complete cases for adjustment: ", nrow(cc), " (", sum(cc$os_event), " deaths)")

models <- rbindlist(list(
  fit_report(Surv(os_months, os_event) ~ group, "univariate (median split)", cl),
  fit_report(Surv(os_months, os_event) ~ group, "univariate, complete cases", cc),
  fit_report(Surv(os_months, os_event) ~ group + age_years, "+ age", cc),
  fit_report(Surv(os_months, os_event) ~ group + age_years + stage, "+ age + stage", cc),
  fit_report(Surv(os_months, os_event) ~ group + age_years + stage + pam50,
             "+ age + stage + PAM50", cc),
  fit_report(Surv(os_months, os_event) ~ expr_z + age_years + stage + pam50,
             "continuous (per SD), fully adjusted", cc)))
fwrite(models, file.path(DIR_OUT, paste0(GENE, "_cox_models.csv")))

msg("\n=== ", GENE, " Cox models (5-year OS) ===")
print(as.data.frame(models[, .(model, n, events, HR = round(HR, 2),
                               CI = sprintf("%.2f-%.2f", lo, hi),
                               p = signif(p, 3))]))

zph <- cox.zph(coxph(Surv(os_months, os_event) ~ group + age_years + stage + pam50, data = cc))
msg("\nPH assumption (global p): ", signif(zph$table["GLOBAL", "p"], 3))

msg("\n=== expression by PAM50 ===")
print(cl[!is.na(pam50), .(n = .N, median_lcpm = round(median(expr), 2)), by = pam50][order(-median_lcpm)])
kw <- kruskal.test(expr ~ pam50, data = cl[!is.na(pam50)])
msg("Kruskal-Wallis across subtypes: p = ", signif(kw$p.value, 3))

p_sub <- ggplot(cl[!is.na(pam50)], aes(x = pam50, y = expr)) +
  geom_boxplot(width = 0.55, outlier.size = 0.5, fill = "#cde2fb",
               colour = "#1c5cab", linewidth = 0.4) +
  labs(title = paste0(GENE, " expression by PAM50 subtype"),
       subtitle = sprintf("TCGA-BRCA primary tumours (n = %d) · Kruskal-Wallis p = %s",
                          sum(!is.na(cl$pam50)), format.pval(kw$p.value, digits = 2, eps = 1e-16)),
       x = NULL, y = expression(log[2]~"CPM")) +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 12.5, colour = "#0b0b0b"),
        plot.subtitle = element_text(size = 9, colour = "#52514e", margin = margin(b = 6)),
        axis.text = element_text(colour = "#52514e"))
ggsave(file.path(DIR_D, paste0(GENE, "_by_pam50.png")), p_sub, width = 7, height = 4.6,
       dpi = 300, bg = "#fcfcfb")

km <- function(d, ttl, sfx) {
  f <- survfit(Surv(os_months, os_event) ~ group, data = d)
  lr <- survdiff(Surv(os_months, os_event) ~ group, data = d)
  pv <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  g <- ggsurvplot(f, data = d, conf.int = TRUE, conf.int.alpha = 0.13,
                  palette = unname(PAL[c("Low","High")]), legend.title = "",
                  legend.labs = c("Low","High"), legend = "top",
                  xlab = "Months from diagnosis", ylab = "Overall survival",
                  xlim = c(0, 60), break.time.by = 12, ylim = c(0, 1),
                  risk.table = TRUE, risk.table.height = 0.25, risk.table.fontsize = 3,
                  censor.size = 2.5, size = 0.7,
                  ggtheme = theme_minimal(base_size = 11) +
                    theme(panel.grid.minor = element_blank(),
                          plot.background = element_rect(fill = "#fcfcfb", colour = NA)),
                  tables.theme = theme_minimal(base_size = 9) +
                    theme(panel.grid = element_blank(),
                          plot.background = element_rect(fill = "#fcfcfb", colour = NA)))
  g$plot <- g$plot + labs(title = ttl,
      subtitle = sprintf("n = %d, %d deaths · log-rank p = %s", nrow(d), sum(d$os_event),
                         format.pval(pv, digits = 2, eps = 1e-4))) +
    theme(plot.title = element_text(face = "bold", size = 12, colour = "#0b0b0b"),
          plot.subtitle = element_text(size = 9, colour = "#52514e"))
  png(file.path(DIR_D, paste0("KM_", GENE, "_", sfx, ".png")), width = 1800, height = 1700, res = 250)
  print(g, newpage = FALSE); dev.off()
  data.table(stratum = sfx, n = nrow(d), events = sum(d$os_event), logrank_p = pv)
}

strata <- list(km(cl, paste0(GENE, " — all tumours"), "all"))
for (s in c("LumA", "LumB", "Her2", "Basal")) {
  d <- cl[pam50 == s]
  d[, group := factor(ifelse(expr > median(expr), "High", "Low"), levels = c("Low","High"))]
  if (nrow(d) > 40 && sum(d$os_event) >= 5)
    strata[[length(strata) + 1]] <- km(d, paste0(GENE, " — ", s), s)
}
st <- rbindlist(strata)
fwrite(st, file.path(DIR_OUT, paste0(GENE, "_subtype_survival.csv")))
msg("\n=== KM by stratum ==="); print(as.data.frame(st))

up <- fread(file.path(DIR_OUT, "RRA_upregulated.csv"))
msg("\n=== cross-cohort effect (from the five-cohort aggregation) ===")
print(as.data.frame(up[symbol == GENE, .SD, .SDcols = patterns("symbol|logFC_|padj")]))

rk <- t(apply(lcpm, 1, rank)); rk <- t(rk)
r  <- sort(cor(rk[, ENSG], rk)[1, ], decreasing = TRUE)
r  <- r[names(r) != ENSG]

ann <- as.data.table(suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = names(r), keytype = "ENSEMBL",
  columns = c("SYMBOL", "MAP"))))[!duplicated(ENSEMBL)]
ann[, chr := sub("[pq].*", "", MAP)]

nb300 <- names(head(r, 300))
same  <- ann[ENSEMBL %in% nb300 & chr == CHR, ENSEMBL]
msg("\n", GENE, " is at ", LOCUS)
msg("same-chromosome genes in the 300-gene neighbourhood: ", length(same),
    " (", round(100 * length(same) / 300, 1), "%)")
msg("same-chromosome genes in the top 30: ",
    sum(ann[ENSEMBL %in% names(head(r, 30)), chr] == CHR, na.rm = TRUE), "/30")

top <- data.table(ensembl = names(head(r, 30)), rho = round(head(r, 30), 3))
top[, symbol := ann$SYMBOL[match(ensembl, ann$ENSEMBL)]]
top[, locus  := ann$MAP[match(ensembl, ann$ENSEMBL)]]
top[, same_chr := ann$chr[match(ensembl, ann$ENSEMBL)] == CHR]
fwrite(top, file.path(DIR_OUT, paste0(GENE, "_top_coexpressed.csv")))
msg("top 15 co-expressed: ", paste(head(top$symbol, 15), collapse = ", "))

r_off <- r[!(names(r) %in% ann[chr == CHR, ENSEMBL])]
nb_off <- names(head(r_off, 300))
top_off <- data.table(ensembl = names(head(r_off, 30)),
                      rho = round(head(r_off, 30), 3))
top_off[, symbol := ann$SYMBOL[match(ensembl, ann$ENSEMBL)]]
top_off[, locus  := ann$MAP[match(ensembl, ann$ENSEMBL)]]
fwrite(top_off, file.path(DIR_OUT, paste0(GENE, "_top_coexpressed_offchrom.csv")))
msg("top 15 OFF-chromosome co-expressed: ",
    paste(head(top_off$symbol, 15), collapse = ", "))
saveRDS(list(nb_raw = nb300, nb_off = nb_off, locus = LOCUS, chr = CHR),
        file.path(DIR_OUT, paste0(GENE, "_neighbourhoods.rds")))
msg("done.")
