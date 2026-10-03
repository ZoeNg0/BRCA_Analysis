#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs

run () {
  echo "=== $1 ${GSE:+($GSE)} ==="
  Rscript "scripts/$1" 2>&1 | tee "logs/${1%.R}${GSE:+_$GSE}.log"
}

run 00_install_deps.R
run tcga_brca_dea.R
run 01_download_geo.R
run 02_de_gse70947.R
run 03_de_gse58135.R
GSE=GSE42568 run 06_de_affy.R
GSE=GSE86374 run 06_de_affy.R
run 04_tcga_ranked_list.R
run 05_rank_aggregation.R
run 07_plot_top_targets.R
run 08_survival.R
