# Bar plot of the top 10 replicated ("agreed") upregulated targets.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source(file.path(getwd(), "R", "utils.R"))

COHORTS <- c("TCGA", "GSE70947", "GSE58135", "GSE42568", "GSE86374")
N_TOP   <- 10

PAL_LIGHT <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
PAL_DARK  <- c("#3987e5", "#d95926", "#199e70", "#c98500", "#d55181")
names(PAL_LIGHT) <- names(PAL_DARK) <- COHORTS

up <- fread(file.path(DIR_OUT, "RRA_upregulated.csv"))

agreed <- up[padj_RRA < 0.05 & n_cohorts_sig == length(COHORTS)][seq_len(N_TOP)]
msg("top ", N_TOP, " agreed targets: ", paste(agreed$symbol, collapse = ", "))

tab <- agreed[, c("symbol", "padj_RRA", paste0("logFC_", COHORTS)), with = FALSE]
setnames(tab, paste0("logFC_", COHORTS), COHORTS)
fwrite(tab, file.path(DIR_OUT, "top10_agreed_targets.csv"))

long <- melt(agreed,
             id.vars = "symbol",
             measure.vars = paste0("logFC_", COHORTS),
             variable.name = "cohort", value.name = "log2FC")
long[, cohort := factor(sub("^logFC_", "", cohort), levels = COHORTS)]
long[, symbol := factor(symbol, levels = rev(agreed$symbol))]

make_plot <- function(pal, surface, ink, ink2, grid_col) {
  ggplot(long, aes(x = log2FC, y = symbol, fill = cohort)) +
    geom_col(width = 0.86,
             position = position_dodge2(padding = 0.34, reverse = TRUE)) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04)),
                       breaks = seq(0, 7, 1)) +
    labs(
      title    = "Top 10 replicated upregulated targets in breast tumour vs normal",
      subtitle = paste0("Ordered by robust rank aggregation. Every gene shown is ",
                        "significant (padj < 0.05) and concordant in all five cohorts."),
      x = expression(log[2]~"fold change (tumour vs normal)"),
      y = NULL,
      caption = "TCGA-BRCA (RNA-seq) · GSE70947 (Agilent, paired) · GSE58135 (RNA-seq) · GSE42568 (Affy U133) · GSE86374 (Affy Gene ST)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.background    = element_rect(fill = surface, colour = NA),
      panel.background   = element_rect(fill = surface, colour = NA),
      legend.position    = "top",
      legend.justification = "left",
      legend.key.size    = unit(9, "pt"),
      legend.text        = element_text(colour = ink2, size = 9),
      legend.margin      = margin(b = 2),
      plot.title         = element_text(colour = ink,  face = "bold", size = 13,
                                        margin = margin(b = 3)),
      plot.subtitle      = element_text(colour = ink2, size = 9.5,
                                        margin = margin(b = 10)),
      plot.caption       = element_text(colour = ink2, size = 7.5, hjust = 0,
                                        margin = margin(t = 10)),
      axis.text.y        = element_text(colour = ink, face = "italic", size = 10),
      axis.text.x        = element_text(colour = ink2, size = 9),
      axis.title.x       = element_text(colour = ink2, size = 9,
                                        margin = margin(t = 6)),
      panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin        = margin(14, 18, 10, 12)
    )
}

p_light <- make_plot(PAL_LIGHT, "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1")
p_dark  <- make_plot(PAL_DARK,  "#1a1a19", "#ffffff", "#c3c2b7", "#383835")

f1 <- file.path(DIR_PLOTS, "top10_agreed_targets.png")
f2 <- file.path(DIR_PLOTS, "top10_agreed_targets_dark.png")
ggsave(f1, p_light, width = 9, height = 6.2, dpi = 300, bg = "#fcfcfb")
ggsave(f2, p_dark,  width = 9, height = 6.2, dpi = 300, bg = "#1a1a19")
msg("wrote ", f1)
msg("wrote ", f2)
