# TCGA-BRCA download, mutation summary and tumour vs normal DESeq2

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(maftools)
  library(EnhancedVolcano)
  library(dplyr)
})

ro_dir    <- "exports/r_objects"
tab_dir   <- "exports/tables"
plots_dir <- "plots"
for (d in c(ro_dir, tab_dir, plots_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

snv_query <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Simple Nucleotide Variation",
                      access = "open",
                      data.type = "Masked Somatic Mutation",
                      data.format = "MAF")
GDCdownload(snv_query, files.per.chunk = 50)
maf <- GDCprepare(snv_query)
saveRDS(maf, file.path(ro_dir, "maf.rds"))

rna_query <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      access = "open",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor", "Solid Tissue Normal"))
GDCdownload(rna_query)
rna <- GDCprepare(rna_query)
rna_matrix <- assay(rna, "unstranded")
saveRDS(rna_matrix, file.path(ro_dir, "rna.rds"))

maf <- read.maf(maf = maf)

png(file.path(plots_dir, "mutationOncoplot.png"), width = 2400, height = 1800, res = 300)
oncoplot(maf = maf, top = 20)
dev.off()

png(file.path(plots_dir, "titvPlot.png"), width = 2400, height = 1800, res = 300)
titv(maf = maf, plot = TRUE, useSyn = TRUE)
dev.off()

dds <- DESeqDataSetFromMatrix(countData = rna_matrix,
                              colData = colData(rna),
                              design = ~ sample_type)
dds$sample_type <- relevel(factor(dds$sample_type), ref = "Solid Tissue Normal")
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
saveRDS(dds, file.path(ro_dir, "dds.rds"))

res_df <- as.data.frame(results(dds))
write.csv(res_df, file.path(tab_dir, "DEA.csv"))

res_df$gene_id <- rownames(res_df)
res_df <- left_join(res_df, as.data.frame(rowData(rna)), by = "gene_id")

png(file.path(plots_dir, "VolcanoPlot.png"), width = 2400, height = 1800, res = 300)
EnhancedVolcano(res_df, lab = res_df$gene_name,
                x = "log2FoldChange", y = "padj",
                pCutoff = 0.05, FCcutoff = 1)
dev.off()
