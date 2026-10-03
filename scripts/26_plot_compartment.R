# Plot module/compartment correlation

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(getwd(), "R", "utils.R"))

R <- fread(file.path(DIR_OUT, "hub_enrichment", "module_compartment_correlation.csv"))
LAB <- c("1" = "Module 1 — ECM / EMT", "2" = "Module 2 — Proliferation\n(where the hubs are)",
         "3" = "Module 3 — Endothelial", "4" = "Module 4 — Immune")
CLAB <- c(fibroblast = "Fibroblast / stroma", epithelial = "Epithelium / tumour cell",
          endothelial = "Endothelium", immune = "Immune")
PAL <- c("Fibroblast / stroma" = "#eb6834", "Epithelium / tumour cell" = "#2a78d6",
         "Endothelium" = "#1baf7a", "Immune" = "#eda100")

R <- R[module %in% names(LAB) | module %in% as.integer(names(LAB))]
R[, mlab := factor(LAB[as.character(module)], levels = LAB)]
R[, clab := factor(CLAB[compartment], levels = CLAB)]
R[, marker := factor(marker, levels = unique(marker[order(clab, -rho)]))]

p <- ggplot(R, aes(x = rho, y = marker, colour = clab)) +
  geom_vline(xintercept = 0, colour = "#52514e", linewidth = 0.4, linetype = "22") +
  geom_segment(aes(x = 0, xend = rho, yend = marker), linewidth = 0.7, alpha = 0.55) +
  geom_point(size = 2.3) +
  facet_grid(clab ~ mlab, scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(limits = c(-0.62, 1), breaks = c(-0.5, 0, 0.5, 1),
                     expand = expansion(mult = 0.03)) +
  scale_colour_manual(values = PAL, name = NULL, guide = "none") +
  labs(
    title = "The EMT signal in TCGA-BRCA is stromal fibroblast content, not tumour-cell EMT",
    subtitle = paste0(
      "Spearman correlation between each module's score (PC1) and marker genes of known cellular origin, across 1,095 TCGA-BRCA\n",
      "tumours. Module 1 tracks fibroblasts (DCN 0.93, LUM 0.89, PDGFRB 0.89, FAP 0.87) while being flat to negative against\n",
      "epithelial markers (EPCAM -0.38, CDH1 -0.15, KRT8 0.00). Carcinoma-cell EMT would track epithelial content, not oppose it."),
    x = "Spearman rho with module score", y = NULL,
    caption = paste0(
      "TCGA-BRCA primary tumours only. Module 2, which holds 133 of the 200 hub genes, runs the other way: it moves against stromal and endothelial\n",
      "content and with EPCAM (0.53) — a programme of the carcinoma cells themselves, which is what a monoculture assay can actually test.")) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background    = element_rect(fill = "#fcfcfb", colour = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#e6e5e1", linewidth = 0.25),
    panel.spacing.x    = unit(0.7, "lines"),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 90, size = 7.6, colour = "#0b0b0b", face = "bold"),
    strip.text.x       = element_text(size = 8.6, colour = "#0b0b0b", face = "bold", lineheight = 1.1),
    plot.title         = element_text(face = "bold", size = 12.5, colour = "#0b0b0b"),
    plot.subtitle      = element_text(size = 8, colour = "#52514e", margin = margin(b = 9)),
    plot.caption       = element_text(size = 7, colour = "#52514e", hjust = 0),
    axis.text.y        = element_text(size = 7.8, colour = "#0b0b0b"),
    axis.text.x        = element_text(colour = "#52514e", size = 8))

f <- file.path(DIR_PLOTS, "module_compartment_correlation.png")
ggsave(f, p, width = 11, height = 6.2, dpi = 300, bg = "#fcfcfb")
msg("wrote ", f)
