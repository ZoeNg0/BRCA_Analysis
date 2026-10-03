# Five-year overall-survival analysis for the replicated candidate targets.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
source(file.path(getwd(), "R", "utils.R"))

GENES <- c("CBX3","FANCD2","GGCT","ESRP1","PPP4C","DLG3","CCT6A","SPINT2",
           "TMED3","SDC1","SLC7A5","RAD1","TUBA1C","SMS","DCTPP1","PTK6",
           "TPD52","RCC1","MAPK13","TOMM34","AP1S1","UNC5B","TRAF2","DYNLT1")

MAX_DAYS <- 5 * 365.25
DIR_SURV <- file.path(DIR_PLOTS, "survival")
dir.create(DIR_SURV, recursive = TRUE, showWarnings = FALSE)

PAL <- c(Low = "#2a78d6", High = "#eb6834")

counts <- readRDS(file.path(PROJ_ROOT, "exports", "r_objects", "rna.rds"))
msg("counts: ", nrow(counts), " genes x ", ncol(counts), " samples")

stype <- substr(colnames(counts), 14, 15)
counts <- counts[, stype == "01", drop = FALSE]
msg("primary tumour samples: ", ncol(counts))

patient <- substr(colnames(counts), 1, 12)
counts  <- counts[, !duplicated(patient), drop = FALSE]
patient <- patient[!duplicated(patient)]
colnames(counts) <- patient
msg("unique patients with expression: ", ncol(counts))

dge  <- DGEList(counts = counts)
dge  <- calcNormFactors(dge, method = "TMM")
lcpm <- cpm(dge, log = TRUE, prior.count = 1)
rownames(lcpm) <- strip_ens_version(rownames(lcpm))

map <- as.data.table(AnnotationDbi::select(
  org.Hs.eg.db, keys = GENES, keytype = "SYMBOL", columns = "ENSEMBL"))
setnames(map, c("SYMBOL", "ENSEMBL"), c("symbol", "ensembl"))
map <- map[!is.na(ensembl) & ensembl %in% rownames(lcpm)]
map[, mean_expr := rowMeans(lcpm[ensembl, , drop = FALSE])]
setorder(map, symbol, -mean_expr)
map <- map[!duplicated(symbol)]
msg("genes resolved to expressed ENSEMBL ids: ", nrow(map), "/", length(GENES))
if (length(setdiff(GENES, map$symbol)))
  warning("unresolved: ", paste(setdiff(GENES, map$symbol), collapse = ", "))

cl <- readRDS(file.path(PROJ_ROOT, "data", "clinical", "TCGA_BRCA_clinical.rds"))
setDT(cl)
cl <- cl[!is.na(vital_status) & vital_status %in% c("Alive", "Dead")]

cl[, os_time  := fifelse(vital_status == "Dead",
                         as.numeric(days_to_death),
                         as.numeric(days_to_last_follow_up))]
cl[, os_event := as.integer(vital_status == "Dead")]
cl <- cl[!is.na(os_time) & os_time > 0]
cl <- cl[!duplicated(submitter_id)]
msg("clinical with usable follow-up: ", nrow(cl),
    "  (deaths: ", sum(cl$os_event), ")")

cl[os_time > MAX_DAYS, os_event := 0L]
cl[, os_time := pmin(os_time, MAX_DAYS)]
cl[, os_months := os_time / 30.4375]
msg("events within 5 years: ", sum(cl$os_event))

shared <- intersect(cl$submitter_id, colnames(lcpm))
cl     <- cl[submitter_id %in% shared]
lcpm   <- lcpm[, cl$submitter_id, drop = FALSE]
msg("matched expression + clinical: ", nrow(cl), " patients, ",
    sum(cl$os_event), " events in 5 years")

summ <- list()

for (i in seq_len(nrow(map))) {
  g   <- map$symbol[i]
  eid <- map$ensembl[i]
  e   <- lcpm[eid, ]

  d <- data.frame(time   = cl$os_months,
                  event  = cl$os_event,
                  expr   = as.numeric(e))
  d$group <- factor(ifelse(d$expr > median(d$expr), "High", "Low"),
                    levels = c("Low", "High"))

  fit <- survfit(Surv(time, event) ~ group, data = d)
  lr  <- survdiff(Surv(time, event) ~ group, data = d)
  p   <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  cox <- coxph(Surv(time, event) ~ group, data = d)
  cs  <- summary(cox)

  summ[[g]] <- data.table(
    gene     = g,
    n        = nrow(d),
    events   = sum(d$event),
    HR       = cs$conf.int[1, "exp(coef)"],
    HR_lo    = cs$conf.int[1, "lower .95"],
    HR_hi    = cs$conf.int[1, "upper .95"],
    logrank_p = p)

  ttl <- sprintf("%s — 5-year overall survival, TCGA-BRCA", g)
  sub <- sprintf("Median split of log2 CPM · HR %.2f (95%% CI %.2f–%.2f) · log-rank p = %s",
                 cs$conf.int[1, "exp(coef)"], cs$conf.int[1, "lower .95"],
                 cs$conf.int[1, "upper .95"],
                 format.pval(p, digits = 2, eps = 1e-4))

  gp <- ggsurvplot(
    fit, data = d,
    conf.int = TRUE, conf.int.alpha = 0.13,
    palette = unname(PAL[c("Low", "High")]),
    legend.title = "", legend.labs = c("Low", "High"),
    legend = "top",
    xlab = "Months from diagnosis", ylab = "Overall survival",
    xlim = c(0, 60), break.time.by = 12,
    ylim = c(0, 1),
    risk.table = TRUE, risk.table.height = 0.24,
    risk.table.title = "Patients at risk",
    risk.table.fontsize = 3.1,
    censor.size = 2.6, size = 0.7,
    ggtheme = theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.25),
            plot.background  = element_rect(fill = "#fcfcfb", colour = NA),
            axis.text  = element_text(colour = "#52514e"),
            axis.title = element_text(colour = "#52514e", size = 9.5),
            legend.text = element_text(colour = "#52514e")),
    tables.theme = theme_minimal(base_size = 9) +
      theme(panel.grid = element_blank(),
            plot.background = element_rect(fill = "#fcfcfb", colour = NA),
            axis.text  = element_text(colour = "#52514e"),
            axis.title = element_blank()))

  gp$plot <- gp$plot +
    labs(title = ttl, subtitle = sub) +
    theme(plot.title    = element_text(face = "bold", size = 12.5, colour = "#0b0b0b"),
          plot.subtitle = element_text(size = 9, colour = "#52514e",
                                       margin = margin(b = 6)))

  f <- file.path(DIR_SURV, paste0("KM_", g, ".png"))
  png(f, width = 1900, height = 1750, res = 260)
  print(gp, newpage = FALSE)
  dev.off()
  msg(sprintf("  %-9s HR %.2f (%.2f-%.2f)  p=%.3g  -> %s",
              g, summ[[g]]$HR, summ[[g]]$HR_lo, summ[[g]]$HR_hi, p, basename(f)))
}

res <- rbindlist(summ)
res[, logrank_fdr := p.adjust(logrank_p, method = "BH")]
setorder(res, logrank_p)
fwrite(res, file.path(DIR_OUT, "survival_5yr_summary.csv"))
msg("wrote ", file.path(DIR_OUT, "survival_5yr_summary.csv"))

msg("\n=== 5-year OS summary (ranked by log-rank p) ===")
print(as.data.frame(res[, .(gene, n, events,
                            HR = round(HR, 2),
                            CI = sprintf("%.2f-%.2f", HR_lo, HR_hi),
                            p = signif(logrank_p, 3),
                            FDR = signif(logrank_fdr, 3))]))

fp <- copy(res)
setorder(fp, HR)
fp[, gene := factor(gene, levels = gene)]
fp[, sig  := fifelse(logrank_fdr < 0.05, "FDR < 0.05", "not significant")]

p_forest <- ggplot(fp, aes(x = HR, y = gene, colour = sig)) +
  geom_vline(xintercept = 1, colour = "#9a9992", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = HR_lo, xmax = HR_hi), height = 0, linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c("FDR < 0.05" = "#2a78d6",
                                 "not significant" = "#9a9992"), name = NULL) +
  scale_x_continuous(trans = "log2", breaks = c(0.75, 1, 1.5, 2, 3),
                     limits = c(0.66, 3.8), expand = expansion(mult = 0.01)) +
  labs(title = "Five-year overall survival by expression, TCGA-BRCA",
       subtitle = paste0("Univariate Cox hazard ratio (high vs low, median split) with 95% CI · ",
                         nrow(cl), " patients, ", sum(cl$os_event), " deaths within 5 years"),
       x = "Hazard ratio (high expression vs low)", y = NULL,
       caption = "HR > 1 means higher expression is associated with worse survival.") +
  theme_minimal(base_size = 11) +
  theme(plot.background  = element_rect(fill = "#fcfcfb", colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#e6e5e1", linewidth = 0.25),
        legend.position  = "top", legend.justification = "left",
        plot.title    = element_text(face = "bold", size = 13, colour = "#0b0b0b"),
        plot.subtitle = element_text(size = 9, colour = "#52514e", margin = margin(b = 8)),
        plot.caption  = element_text(size = 8, colour = "#52514e", hjust = 0),
        axis.text.y   = element_text(face = "italic", colour = "#0b0b0b"),
        axis.text.x   = element_text(colour = "#52514e"),
        axis.title.x  = element_text(colour = "#52514e", size = 9))

f_forest <- file.path(DIR_PLOTS, "survival_5yr_forest.png")
ggsave(f_forest, p_forest, width = 7.5, height = 6.4, dpi = 300, bg = "#fcfcfb")
msg("wrote ", f_forest)

pdf_f <- file.path(DIR_SURV, "KM_all_genes.pdf")
pdf(pdf_f, width = 7.3, height = 6.7, onefile = TRUE)
print(p_forest)
for (g in res$gene) {
  img <- file.path(DIR_SURV, paste0("KM_", g, ".png"))
  if (file.exists(img)) {
    r <- png::readPNG(img)
    grid::grid.newpage(); grid::grid.raster(r)
  }
}
invisible(dev.off())
msg("wrote ", pdf_f)
