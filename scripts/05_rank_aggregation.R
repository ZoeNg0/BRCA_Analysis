# Cross-cohort rank aggregation.

suppressPackageStartupMessages({
  library(RobustRankAggreg)
  library(data.table)
})
source(file.path(getwd(), "R", "utils.R"))

COHORTS <- strsplit(Sys.getenv("COHORTS",
            "TCGA,GSE70947,GSE58135,GSE42568,GSE86374"), ",")[[1]]
msg("aggregating over: ", paste(COHORTS, collapse = ", "))

de <- lapply(COHORTS, read_de)
names(de) <- COHORTS

for (n in names(de)) {
  d <- de[[n]]
  msg(sprintf("  %-10s %6d genes | up(padj<.05,lfc>1)=%5d down=%5d", n, nrow(d),
              sum(d$padj < 0.05 & d$logFC >  1, na.rm = TRUE),
              sum(d$padj < 0.05 & d$logFC < -1, na.rm = TRUE)))
}

common <- Reduce(intersect, lapply(de, `[[`, "symbol"))
msg("genes measured in all ", length(de), " cohorts: ", length(common))
if (length(common) < 1000)
  warning("suspiciously small common gene set — check symbol mapping")

de <- lapply(de, function(d) d[match(common, d$symbol), ])

score <- sapply(de, function(d) {
  p <- pmax(d$pvalue, .Machine$double.xmin)
  sign(d$logFC) * -log10(p)
})
rownames(score) <- common

aggregate_dir <- function(score, direction = c("up", "down")) {
  direction <- match.arg(direction)
  s <- if (direction == "up") score else -score
  lists <- lapply(seq_len(ncol(s)), function(j) rownames(s)[order(s[, j], decreasing = TRUE)])
  names(lists) <- colnames(s)
  ag <- aggregateRanks(glist = lists, N = nrow(s), method = "RRA", full = TRUE)
  ag <- as.data.table(ag)
  setnames(ag, c("Name", "Score"), c("symbol", "RRA_score"))
  ag[, padj_RRA := p.adjust(RRA_score, method = "BH")]
  for (n in names(de)) {
    d <- de[[n]]
    i <- match(ag$symbol, d$symbol)
    set(ag, j = paste0("logFC_", n), value = d$logFC[i])
    set(ag, j = paste0("padj_",  n), value = d$padj[i])
  }
  lfc <- as.matrix(ag[, paste0("logFC_", names(de)), with = FALSE])
  pad <- as.matrix(ag[, paste0("padj_",  names(de)), with = FALSE])
  want <- if (direction == "up") 1 else -1
  ag[, n_cohorts_sig := rowSums(sign(lfc) == want & pad < 0.05, na.rm = TRUE)]
  ag[, n_cohorts_dir := rowSums(sign(lfc) == want, na.rm = TRUE)]
  setorder(ag, RRA_score)
  ag[]
}

up   <- aggregate_dir(score, "up")
down <- aggregate_dir(score, "down")

f_up   <- file.path(DIR_OUT, "RRA_upregulated.csv")
f_down <- file.path(DIR_OUT, "RRA_downregulated.csv")
fwrite(up,   f_up)
fwrite(down, f_down)
msg("wrote ", f_up)
msg("wrote ", f_down)

msg("RRA-significant (padj_RRA<0.05): up=", sum(up$padj_RRA < 0.05),
    " down=", sum(down$padj_RRA < 0.05))

msg("\n=== top 25 replicated UPREGULATED genes ===")
print(as.data.frame(head(up[, c("symbol", "padj_RRA",
                  paste0("logFC_", names(de)), "n_cohorts_sig"), with = FALSE], 25)))

msg("\n=== top 15 replicated DOWNREGULATED genes ===")
print(as.data.frame(head(down[, c("symbol", "padj_RRA",
                    paste0("logFC_", names(de)), "n_cohorts_sig"), with = FALSE], 15)))

short <- up[padj_RRA < 0.05 & n_cohorts_sig == length(de)]
fwrite(short, file.path(DIR_OUT, "shortlist_upregulated.csv"))
msg("\nshortlist (RRA padj<0.05 AND significant+concordant in all ",
    length(de), " cohorts): ", nrow(short), " genes")
