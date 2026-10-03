# Shared helpers

proj_dir  <- function() {
  d <- getwd()
  while (!file.exists(file.path(d, "R/utils.R")) && dirname(d) != d) d <- dirname(d)
  if (file.exists(file.path(d, "R/utils.R"))) d else getwd()
}

PROJ_ROOT <- proj_dir()
DIR_DATA  <- file.path(PROJ_ROOT, "data", "geo")
DIR_OUT   <- file.path(PROJ_ROOT, "exports", "validation")
DIR_PLOTS <- file.path(PROJ_ROOT, "plots")
for (d in c(DIR_DATA, DIR_OUT, DIR_PLOTS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

strip_ens_version <- function(x) sub("\\..*$", "", x)

collapse_by_symbol <- function(df, symbol_col = "symbol", stat_col = "logFC") {
  df <- df[!is.na(df[[symbol_col]]) & nzchar(df[[symbol_col]]), , drop = FALSE]
  ord <- order(abs(df[[stat_col]]), decreasing = TRUE)
  df  <- df[ord, , drop = FALSE]
  df[!duplicated(df[[symbol_col]]), , drop = FALSE]
}

write_de <- function(df, cohort) {
  stopifnot(all(c("symbol", "logFC", "pvalue", "padj") %in% names(df)))
  df <- df[order(df$pvalue), c("symbol", "logFC", "pvalue", "padj")]
  f  <- file.path(DIR_OUT, paste0("DE_", cohort, ".csv"))
  write.csv(df, f, row.names = FALSE)
  msg("wrote ", f, "  (", nrow(df), " genes)")
  invisible(df)
}

read_de <- function(cohort) {
  f <- file.path(DIR_OUT, paste0("DE_", cohort, ".csv"))
  if (!file.exists(f)) stop("missing DE result: ", f, " — run the earlier scripts first")
  read.csv(f, stringsAsFactors = FALSE)
}
