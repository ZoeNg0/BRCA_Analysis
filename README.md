# BRCA Analysis

R scripts for a breast cancer (TCGA-BRCA) analysis: differential expression in TCGA, replication across four GEO cohorts, survival analysis, and co-expression hub gene / pathway enrichment.

## Data

| Cohort | Source | Platform |
|---|---|---|
| TCGA-BRCA | GDC via TCGAbiolinks | RNA-seq (STAR counts), MAF |
| GSE70947 | GEO | Microarray, paired |
| GSE58135 | GEO | RNA-seq |
| GSE42568 | GEO | Affymetrix GPL570 |
| GSE86374 | GEO | Affymetrix GPL6244 |

Raw data is not included in the repo; the scripts download it.

## Scripts

- `tcga_brca_dea.R` - TCGA download, mutation plots, DESeq2 tumour vs normal
- `01`-`06` - GEO download and DE per cohort, cross-cohort rank aggregation
- `07`-`12` - top targets, 5-year survival, pathway enrichment, per-gene summaries
- `13`-`18` - PPP4C figures: cohorts, PAM50 subtypes, Kaplan-Meier
- `19`-`27` - TCGA co-expression hubs, enrichment, proliferation vs EMT

`R/utils.R` holds shared helpers.

## Running

```bash
Rscript scripts/00_install_deps.R
bash scripts/run_all.sh
```

Scripts from `09` onwards are run individually. Outputs go to `exports/` and `plots/`.
