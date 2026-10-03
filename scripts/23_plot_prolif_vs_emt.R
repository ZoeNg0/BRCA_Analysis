# Plot proliferation vs EMT

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(getwd(), "R", "utils.R"))

DIR_EN <- file.path(DIR_OUT, "hub_enrichment")
d <- fread(file.path(DIR_EN, "directed_emt_vs_proliferation.csv"))
d <- d[!is.na(overlap) & theme != "other"]

pretty <- function(x) {
  x <- sub("^(HALLMARK|GOBP|REACTOME|WP|PID|BIOCARTA)_", "", x)
  x <- tools::toTitleCase(tolower(gsub("_", " ", x)))
  x <- sub("Epithelial Mesenchymal Transition", "EMT (Hallmark)", x)
  x <- sub("Epithelial to Mesenchymal Transition", "Epithelial-to-Mesenchymal Transition", x)
  x <- gsub("\\bDna\\b", "DNA", x); x <- gsub("\\bG2m\\b", "G2M", x)
  x <- gsub("\\bMyc\\b", "MYC", x); x <- gsub("\\bE2f\\b", "E2F", x)
  x <- gsub("Tgf Beta", "TGF-beta", x); x <- gsub("\\bIntegrin1\\b", "Integrin-1", x)
  x <- gsub("\\bCellcycle\\b", "Cell Cycle", x)
  gsub("Cell Cell Adhesion", "Cell-Cell Adhesion", x)
}
d[, label := pretty(term)]
d[, src := sub("_.*$", "", term)]
d[, label := paste0(label, "  [", src, "]")]
d[, fold_plot := fifelse(overlap > 0, fold, NA_real_)]
d[, sig := fifelse(padj < 0.05, "FDR < 0.05", "not significant")]
setorder(d, theme, fold)
d[, label := factor(label, levels = label)]

stars <- function(p) fifelse(is.na(p), "", fifelse(p < 1e-20, "FDR<1e-20",
          fifelse(p < 0.05, sprintf("FDR=%.2g", p), sprintf("ns (FDR=%.2g)", p))))
d[, ann := stars(padj)]
d[overlap == 0, ann := "0 hub genes"]
d[, ann := paste0(overlap, "/200  ", ann)]
d[, ann_x := fifelse(overlap > 0, fold, 1)]

pal <- c("proliferation" = "#2a78d6", "EMT / motility" = "#eb6834")

p <- ggplot(d, aes(x = fold_plot, y = label, fill = theme, alpha = sig)) +
  geom_col(width = 0.68, na.rm = TRUE) +
  geom_vline(xintercept = 1, colour = "#52514e", linewidth = 0.4, linetype = "22") +
  geom_text(aes(x = ann_x, label = ann), hjust = -0.08, size = 2.5,
            colour = "#52514e", alpha = 1, show.legend = FALSE) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 25),
                labels = c("0.5x", "1x", "2x", "5x", "10x", "25x"),
                limits = c(0.4, 60),
                expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "not significant" = 0.35), name = NULL) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title = "The TCGA-BRCA hub programme is proliferative, not motile",
    x = "Fold enrichment over expectation (log scale)", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background   = element_rect(fill = "#fcfcfb", colour = NA),
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#e6e5e1", linewidth = 0.25),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 9, colour = "#0b0b0b"),
    plot.title        = element_text(face = "bold", size = 13, colour = "#0b0b0b",
                                     margin = margin(b = 10)),
    axis.text.y       = element_text(size = 7.6, colour = "#0b0b0b"),
    axis.text.x       = element_text(colour = "#52514e"),
    legend.position   = "top",
    legend.justification = "left",
    legend.margin     = margin(b = -4))

f <- file.path(DIR_PLOTS, "hub_proliferation_vs_emt.png")
ggsave(f, p, width = 9.6, height = 6.9, dpi = 300, bg = "#fcfcfb")
msg("wrote ", f)
