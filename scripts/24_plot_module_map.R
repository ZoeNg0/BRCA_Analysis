# Module vs hallmark map

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source(file.path(getwd(), "R", "utils.R"))

DIR_EN <- file.path(DIR_OUT, "hub_enrichment")
R <- fread(file.path(DIR_EN, "module_enrichment_all_databases.csv"))
H <- fread(file.path(DIR_OUT, "network", "module_intra_degree.csv"))
net <- readRDS(file.path(DIR_OUT, "network", "tcga_network.rds"))
hub <- as.data.table(net$hub)

LAB <- c("1" = "Module 1\nECM / stromal", "2" = "Module 2\nProliferation",
         "3" = "Module 3\nEndothelial", "4" = "Module 4\nImmune",
         "5" = "Module 5\nBasal epithelial", "9" = "Module 9\nMitochondrial")

hm <- R[database == "MSigDB_Hallmark"]
hm[, module := as.integer(sub("module_", "", gene_set))]
hm[, term := tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", Description))))]
fix <- c("E2f"="E2F", "G2m"="G2M", "Il2"="IL2", "Stat5"="STAT5", "Kras"="KRAS",
         "Tnfa"="TNFA", "Nfkb"="NF-kB", "Myc"="MYC", "Dna"="DNA", "Mtorc1"="mTORC1")
for (k in names(fix)) hm[, term := gsub(paste0("\\b", k, "\\b"), fix[[k]], term)]
keep <- hm[, .(best = min(p.adjust)), by = term][order(best)][1:16, term]
hp <- hm[term %in% keep]

grid <- CJ(module = as.integer(names(LAB)), term = keep, unique = TRUE)
hp <- merge(grid, hp[, .(module, term, p.adjust)], by = c("module", "term"), all.x = TRUE)
hp[, nlp := -log10(p.adjust)]

n_hub <- hub[is_hub == TRUE, .N, by = module][!is.na(module)]
sizes <- H[, .N, by = module]
sizes <- merge(sizes, n_hub, by = "module", all.x = TRUE, suffixes = c("_genes", "_hub"))
sizes[is.na(N_hub), N_hub := 0]
sizes[, lab := sprintf("%s\n%d genes · %d of the 200 hubs",
                       sub("\n", " — ", LAB[as.character(module)]), N_genes, N_hub)]

hub_n0 <- setNames(sizes$N_hub, sizes$module)
LAB2 <- setNames(sprintf("%s\n%d hubs", LAB, hub_n0[names(LAB)]), names(LAB))
hp[, mlab := factor(LAB2[as.character(module)], levels = LAB2)]
hp[, term := factor(term, levels = rev(keep))]
hub_n <- setNames(sizes$N_hub, sizes$module)
hp[, hub_count := hub_n[as.character(module)]]

p <- ggplot(hp, aes(x = mlab, y = term, fill = nlp)) +
  geom_tile(colour = "#fcfcfb", linewidth = 1.6) +
  geom_text(data = hp[!is.na(nlp) & nlp > 12], aes(label = sprintf("%.0f", nlp)),
            colour = "#fcfcfb", size = 2.7, fontface = "bold") +
  scale_fill_gradient(low = "#cde2fb", high = "#184f95", na.value = "#f2f1ed",
                      name = expression(-log[10]~"FDR")) +
  scale_x_discrete(position = "top") +
  labs(
    title = "Proliferation and EMT are separate modules of the same TCGA-BRCA network",
    subtitle = paste0(
      "Hallmark enrichment of each co-expression module. All 200 degree hubs fall in modules 1-4, and 133 of them\n",
      "in module 2 alone. Module 1 carries the EMT signal but contributes no proliferation biology — and vice versa.\n",
      "Grey = not enriched (FDR > 0.05)."),
    x = NULL, y = NULL,
    caption = "TCGA-BRCA only: DEGs (padj<0.05, |log2FC|>1) co-expressed at |Spearman rho| >= 0.7 across 1,095 tumours; Louvain modules.") +
  theme_minimal(base_size = 10) +
  theme(
    plot.background   = element_rect(fill = "#fcfcfb", colour = NA),
    panel.grid        = element_blank(),
    plot.title        = element_text(face = "bold", size = 13, colour = "#0b0b0b"),
    plot.subtitle     = element_text(size = 8.2, colour = "#52514e", margin = margin(b = 10)),
    plot.caption      = element_text(size = 7, colour = "#52514e", hjust = 0),
    axis.text.x.top   = element_text(size = 8.4, colour = "#0b0b0b", face = "bold", lineheight = 1.1),
    axis.text.y       = element_text(size = 8.6, colour = "#0b0b0b"),
    legend.position   = "right",
    legend.key.height = unit(1.1, "cm"))

f <- file.path(DIR_PLOTS, "tcga_module_hallmark_map.png")
ggsave(f, p, width = 10.4, height = 6.6, dpi = 300, bg = "#fcfcfb")
msg("wrote ", f)

print(as.data.frame(sizes[order(-N_hub), .(module, N_genes, N_hub)]))
