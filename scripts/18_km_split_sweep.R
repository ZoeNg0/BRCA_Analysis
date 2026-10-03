# PPP4C 5-year KM across a range of split cutoffs

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(survival); library(survminer); library(ggplot2)
  library(TCGAbiolinks)
})
source(file.path(getwd(), "R", "utils.R"))

GENE     <- "PPP4C"
ENSG     <- "ENSG00000149923"
MAX_DAYS <- 5 * 365.25
PAL      <- c(Low = "#2a78d6", High = "#eb6834")
DIR_S    <- file.path(DIR_PLOTS, tolower(GENE), "sweep")
dir.create(DIR_S, recursive = TRUE, showWarnings = FALSE)

SPLITS <- list(
  list(lo = 0.50, hi = 0.50, key = "50_50", label = "50 / 50  (median)"),
  list(lo = 0.40, hi = 0.60, key = "60_40", label = "60 / 40"),
  list(lo = 0.30, hi = 0.70, key = "70_30", label = "70 / 30"),
  list(lo = 0.20, hi = 0.80, key = "80_20", label = "80 / 20"),
  list(lo = 0.10, hi = 0.90, key = "90_10", label = "90 / 10"))

counts  <- readRDS(file.path(PROJ_ROOT, "exports", "r_objects", "rna.rds"))
counts  <- counts[, substr(colnames(counts), 14, 15) == "01", drop = FALSE]
patient <- substr(colnames(counts), 1, 12)
counts  <- counts[, !duplicated(patient), drop = FALSE]
colnames(counts) <- patient[!duplicated(patient)]
lcpm <- cpm(calcNormFactors(DGEList(counts = counts), method = "TMM"),
            log = TRUE, prior.count = 1)
rownames(lcpm) <- strip_ens_version(rownames(lcpm))
i <- which(rownames(lcpm) == ENSG); stopifnot(length(i) == 1)

cl <- as.data.table(readRDS(file.path(PROJ_ROOT, "data", "clinical",
                                      "TCGA_BRCA_clinical.rds")))
cl <- cl[vital_status %in% c("Alive", "Dead")]
cl[, os_time  := fifelse(vital_status == "Dead", as.numeric(days_to_death),
                                                 as.numeric(days_to_last_follow_up))]
cl[, os_event := as.integer(vital_status == "Dead")]
cl <- cl[!is.na(os_time) & os_time > 0][!duplicated(submitter_id)]
cl[os_time > MAX_DAYS, os_event := 0L]
cl[, os_time := pmin(os_time, MAX_DAYS)]
cl[, os_months := os_time / 30.4375]
cl <- cl[submitter_id %in% colnames(lcpm)]
cl[, expr := as.numeric(lcpm[i, submitter_id])]
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
msg("analysis set: ", nrow(cl), " patients, ", sum(cl$os_event), " deaths in 5 yr")

curves <- list(); summ <- list()

for (s in SPLITS) {
  q_lo <- as.numeric(quantile(cl$expr, s$lo, type = 7))
  q_hi <- as.numeric(quantile(cl$expr, s$hi, type = 7))
  d <- copy(cl)
  d[, group := fcase(expr <= q_lo, "Low", expr >= q_hi, "High", default = NA_character_)]
  n_mid <- sum(is.na(d$group))
  d <- d[!is.na(group)]
  d[, group := factor(group, levels = c("Low", "High"))]

  fit <- survfit(Surv(os_months, os_event) ~ group, data = d)
  lr  <- survdiff(Surv(os_months, os_event) ~ group, data = d)
  p   <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  cs  <- summary(coxph(Surv(os_months, os_event) ~ group, data = d))

  cc <- d[!is.na(age_years) & !is.na(stage) & !is.na(pam50)]
  ca <- tryCatch(summary(coxph(Surv(os_months, os_event) ~ group + age_years + stage + pam50,
                               data = cc)), error = function(e) NULL)

  summ[[s$key]] <- data.table(
    split = s$label, n = nrow(d), dropped = n_mid,
    n_low = sum(d$group == "Low"), n_high = sum(d$group == "High"),
    ev_low = sum(d$os_event[d$group == "Low"]), ev_high = sum(d$os_event[d$group == "High"]),
    HR = cs$conf.int[1, "exp(coef)"], lo = cs$conf.int[1, "lower .95"],
    hi = cs$conf.int[1, "upper .95"], logrank_p = p,
    HR_adj = if (is.null(ca)) NA_real_ else ca$conf.int[1, "exp(coef)"],
    lo_adj = if (is.null(ca)) NA_real_ else ca$conf.int[1, "lower .95"],
    hi_adj = if (is.null(ca)) NA_real_ else ca$conf.int[1, "upper .95"],
    p_adj  = if (is.null(ca)) NA_real_ else ca$coefficients[1, "Pr(>|z|)"],
    n_adj  = nrow(cc))

  st <- rep(sub("^group=", "", names(fit$strata)), fit$strata)
  cd <- data.table(split = s$label, group = st, t = fit$time, surv = fit$surv,
                   lower = fit$lower, upper = fit$upper, ncens = fit$n.censor)
  cd[is.na(lower), lower := surv]; cd[is.na(upper), upper := surv]
  zero <- data.table(split = s$label, group = c("Low", "High"), t = 0, surv = 1,
                     lower = 1, upper = 1, ncens = 0)
  end  <- cd[, .SD[which.max(t)], by = group][, .(split = s$label, group, t = 60,
                                                  surv, lower, upper, ncens = 0)]
  curves[[s$key]] <- rbind(zero, cd, end)

  ttl <- sprintf("%s — 5-year overall survival, TCGA-BRCA", GENE)
  sb  <- if (s$lo == s$hi)
    sprintf("Median split (nothing excluded)\nn = %d · HR %.2f (95%% CI %.2f-%.2f) · log-rank p = %s",
            nrow(d), cs$conf.int[1,1], cs$conf.int[1,3], cs$conf.int[1,4],
            format.pval(p, digits = 2, eps = 1e-4))
  else
    sprintf("Bottom %d%% vs top %d%% of log2 CPM (middle %d%% excluded, n = %d)\nn = %d · HR %.2f (95%% CI %.2f-%.2f) · log-rank p = %s",
            round(100*s$lo), round(100*(1-s$hi)), round(100*(s$hi-s$lo)), n_mid,
            nrow(d), cs$conf.int[1,1], cs$conf.int[1,3], cs$conf.int[1,4],
            format.pval(p, digits = 2, eps = 1e-4))

  gp <- ggsurvplot(fit, data = d, conf.int = TRUE, conf.int.alpha = 0.13,
    palette = unname(PAL[c("Low","High")]), legend.title = "",
    legend.labs = c("Low","High"), legend = "top",
    xlab = "Months from diagnosis", ylab = "Overall survival",
    xlim = c(0, 60), break.time.by = 12, ylim = c(0, 1),
    risk.table = TRUE, risk.table.height = 0.24, risk.table.title = "Patients at risk",
    risk.table.fontsize = 3.1, censor.size = 2.6, size = 0.7,
    ggtheme = theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.25),
            plot.background = element_rect(fill = "#fcfcfb", colour = NA),
            axis.text = element_text(colour = "#52514e"),
            axis.title = element_text(colour = "#52514e", size = 9.5),
            legend.text = element_text(colour = "#52514e")),
    tables.theme = theme_minimal(base_size = 9) +
      theme(panel.grid = element_blank(),
            plot.background = element_rect(fill = "#fcfcfb", colour = NA),
            axis.text = element_text(colour = "#52514e"), axis.title = element_blank()))
  gp$plot <- gp$plot + labs(title = ttl, subtitle = sb) +
    theme(plot.title = element_text(face = "bold", size = 12.5, colour = "#0b0b0b"),
          plot.subtitle = element_text(size = 8.6, colour = "#52514e", margin = margin(b = 6)))
  f <- file.path(DIR_S, sprintf("KM_%s_%s.png", GENE, s$key))
  png(f, width = 1900, height = 1750, res = 260); print(gp, newpage = FALSE); dev.off()
}

res <- rbindlist(summ)
res[, split := factor(split, levels = sapply(SPLITS, `[[`, "label"))]
fwrite(res, file.path(DIR_OUT, "PPP4C_split_sweep.csv"))

msg("\n=== split sweep ===")
print(as.data.frame(res[, .(split, n, dropped,
                            arms = sprintf("%d/%d", n_low, n_high),
                            deaths = sprintf("%d/%d", ev_low, ev_high),
                            HR = round(HR, 2), CI = sprintf("%.2f-%.2f", lo, hi),
                            p = signif(logrank_p, 3),
                            HR_adj = round(HR_adj, 2),
                            p_adj = signif(p_adj, 3))]))

cv <- rbindlist(curves)
cv[, split := factor(split, levels = sapply(SPLITS, `[[`, "label"))]
cv[, group := factor(group, levels = c("Low", "High"))]
lab <- res[, .(split, txt = sprintf("HR %.2f (%.2f-%.2f)\np = %s\nn = %d",
                                    HR, lo, hi, format.pval(logrank_p, digits = 2, eps = 1e-4), n))]

p_panel <- ggplot(cv, aes(t, surv, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.13, colour = NA) +
  geom_step(linewidth = 0.65) +
  geom_point(data = cv[ncens > 0], shape = 124, size = 1.1, show.legend = FALSE) +
  geom_text(data = lab, aes(x = 2, y = 0.625, label = txt), inherit.aes = FALSE,
            hjust = 0, vjust = 0, size = 2.9, colour = "#52514e", lineheight = 1.05) +
  facet_wrap(~ split, nrow = 1) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_fill_manual(values = PAL, name = NULL) +
  scale_x_continuous(breaks = seq(0, 60, 12), limits = c(0, 60)) +
  coord_cartesian(ylim = c(0.6, 1)) +
  labs(title = paste0(GENE, " — 5-year overall survival across expression-split thresholds"),
       subtitle = paste("Low vs High arms only; the discarded middle fraction grows left to right.",
                        "y-axis truncated at 0.6 to show the separation."),
       x = "Months from diagnosis", y = "Overall survival") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.25),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        strip.text = element_text(face = "bold", size = 10, colour = "#0b0b0b"),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 13, colour = "#0b0b0b"),
        plot.subtitle = element_text(size = 8.8, colour = "#52514e", margin = margin(b = 8)),
        axis.text = element_text(colour = "#52514e"),
        axis.title = element_text(colour = "#52514e", size = 9.5))
ggsave(file.path(DIR_S, "KM_sweep_panel.png"), p_panel, width = 15, height = 4.6,
       dpi = 240, bg = "#fcfcfb")

fo <- rbind(
  res[, .(split, model = "univariate", HR, lo, hi, p = logrank_p)],
  res[, .(split, model = "adjusted (age + stage + PAM50)", HR = HR_adj,
          lo = lo_adj, hi = hi_adj, p = p_adj)])
fo[, model := factor(model, levels = c("univariate", "adjusted (age + stage + PAM50)"))]

p_forest <- ggplot(fo, aes(x = HR, y = split, colour = model)) +
  geom_vline(xintercept = 1, colour = "#9a9992", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.7,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 2.3, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c("univariate" = "#2a78d6",
                                 "adjusted (age + stage + PAM50)" = "#9a9992"), name = NULL) +
  scale_x_continuous(trans = "log2", breaks = c(0.75, 1, 1.5, 2, 3, 4, 6)) +
  labs(title = paste0(GENE, " — hazard ratio by expression-split threshold"),
       subtitle = paste("Wider CIs to the right: more extreme splits discard more patients.",
                        "Vertical line = no effect."),
       x = "Hazard ratio (High vs Low, log scale)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 12.5, colour = "#0b0b0b"),
        plot.subtitle = element_text(size = 8.8, colour = "#52514e", margin = margin(b = 8)),
        axis.text = element_text(colour = "#52514e"),
        axis.title = element_text(colour = "#52514e", size = 9.5))
ggsave(file.path(DIR_S, "KM_sweep_forest.png"), p_forest, width = 8.6, height = 4.4,
       dpi = 260, bg = "#fcfcfb")

msg("\nwrote ", DIR_S, "/  (5 KM plots + panel + forest)")
msg("wrote ", file.path(DIR_OUT, "PPP4C_split_sweep.csv"))
