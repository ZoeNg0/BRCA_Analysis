# PPP4C 5-year KM with extreme-group split

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(survival); library(survminer); library(ggplot2)
})
source(file.path(getwd(), "R", "utils.R"))

GENE     <- "PPP4C"
ENSG     <- "ENSG00000149923"
MAX_DAYS <- 5 * 365.25
LOW_Q    <- as.numeric(Sys.getenv("LOW_Q",  "0.40"))
HIGH_Q   <- as.numeric(Sys.getenv("HIGH_Q", "0.60"))
stopifnot(LOW_Q > 0, HIGH_Q < 1, LOW_Q <= HIGH_Q)
PAL      <- c(Low = "#2a78d6", High = "#eb6834")
DIR_D    <- file.path(DIR_PLOTS, tolower(GENE))
dir.create(DIR_D, recursive = TRUE, showWarnings = FALSE)

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
msg("analysis set: ", nrow(cl), " patients, ", sum(cl$os_event), " deaths in 5 yr")

q_lo <- as.numeric(quantile(cl$expr, LOW_Q,  type = 7))
q_hi <- as.numeric(quantile(cl$expr, HIGH_Q, type = 7))
cl[, group := fcase(expr <= q_lo, "Low", expr >= q_hi, "High", default = NA_character_)]
d <- cl[!is.na(group)]
d[, group := factor(group, levels = c("Low", "High"))]
n_mid <- nrow(cl) - nrow(d)
msg(sprintf("cutoffs: p%02d = %.4f, p%02d = %.4f log2 CPM",
            round(100 * LOW_Q), q_lo, round(100 * HIGH_Q), q_hi))
msg(sprintf("Low n=%d (%d deaths) | High n=%d (%d deaths) | middle %d%% dropped: n=%d",
            sum(d$group == "Low"),  sum(d$os_event[d$group == "Low"]),
            sum(d$group == "High"), sum(d$os_event[d$group == "High"]),
            round(100 * (HIGH_Q - LOW_Q)), n_mid))

fit <- survfit(Surv(os_months, os_event) ~ group, data = d)
lr  <- survdiff(Surv(os_months, os_event) ~ group, data = d)
p   <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
cs  <- summary(coxph(Surv(os_months, os_event) ~ group, data = d))
hr  <- cs$conf.int[1, "exp(coef)"]; lo <- cs$conf.int[1, "lower .95"]
hi  <- cs$conf.int[1, "upper .95"]

m <- copy(cl)
m[, group := factor(ifelse(expr > median(expr), "High", "Low"), levels = c("Low", "High"))]
lr_m <- survdiff(Surv(os_months, os_event) ~ group, data = m)
p_m  <- pchisq(lr_m$chisq, df = 1, lower.tail = FALSE)
cs_m <- summary(coxph(Surv(os_months, os_event) ~ group, data = m))

cmp <- data.frame(
  split   = c("median (50/50)", sprintf("extreme (p%02d vs p%02d)",
                                        round(100*LOW_Q), round(100*HIGH_Q))),
  n       = c(nrow(m), nrow(d)),
  events  = c(sum(m$os_event), sum(d$os_event)),
  HR      = c(cs_m$conf.int[1, "exp(coef)"], hr),
  CI      = c(sprintf("%.2f-%.2f", cs_m$conf.int[1, "lower .95"], cs_m$conf.int[1, "upper .95"]),
              sprintf("%.2f-%.2f", lo, hi)),
  logrank = c(p_m, p))
msg("\n=== median split vs extreme split ===")
print(transform(cmp, HR = round(HR, 2), logrank = signif(logrank, 3)))

ttl <- sprintf("%s — 5-year overall survival, TCGA-BRCA", GENE)
sub <- sprintf(paste0("Extreme split: bottom %d%% vs top %d%% of log2 CPM (middle %d%% ",
                      "excluded, n = %d)\nn = %d · HR %.2f (95%% CI %.2f-%.2f) · log-rank p = %s"),
               round(100 * LOW_Q), round(100 * (1 - HIGH_Q)),
               round(100 * (HIGH_Q - LOW_Q)), n_mid,
               nrow(d), hr, lo, hi, format.pval(p, digits = 2, eps = 1e-4))

gp <- ggsurvplot(
  fit, data = d, conf.int = TRUE, conf.int.alpha = 0.13,
  palette = unname(PAL[c("Low", "High")]),
  legend.title = "", legend.labs = c("Low", "High"), legend = "top",
  xlab = "Months from diagnosis", ylab = "Overall survival",
  xlim = c(0, 60), break.time.by = 12, ylim = c(0, 1),
  risk.table = TRUE, risk.table.height = 0.24,
  risk.table.title = "Patients at risk", risk.table.fontsize = 3.1,
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
          axis.text  = element_text(colour = "#52514e"), axis.title = element_blank()))

gp$plot <- gp$plot + labs(title = ttl, subtitle = sub) +
  theme(plot.title    = element_text(face = "bold", size = 12.5, colour = "#0b0b0b"),
        plot.subtitle = element_text(size = 8.6, colour = "#52514e", margin = margin(b = 6)))

f <- file.path(DIR_D, sprintf("KM_%s_p%02d_p%02d.png", GENE,
                              round(100 * LOW_Q), round(100 * HIGH_Q)))
png(f, width = 1900, height = 1750, res = 260)
print(gp, newpage = FALSE)
dev.off()
msg("\nwrote ", f)
