# Download GEO validation cohorts to data/geo/

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
})
source(file.path(getwd(), "R", "utils.R"))

options(timeout = 3600)

gse70947_rds <- file.path(DIR_DATA, "GSE70947_eset.rds")

if (file.exists(gse70947_rds)) {
  msg("GSE70947: cached, skipping download")
} else {
  msg("GSE70947: downloading series matrix (~150 MB uncompressed) ...")
  es <- getGEO("GSE70947", GSEMatrix = TRUE, AnnotGPL = FALSE,
               destdir = DIR_DATA)
  if (is.list(es)) es <- es[[1]]
  msg("GSE70947: ", nrow(es), " features x ", ncol(es), " samples")
  saveRDS(es, gse70947_rds)
  msg("GSE70947: cached -> ", gse70947_rds)
}

soft_gz <- file.path(DIR_DATA, "GSE58135_family.soft.gz")
if (!file.exists(soft_gz)) {
  msg("GSE58135: fetching SOFT record ...")
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58135/soft/GSE58135_family.soft.gz",
    soft_gz, mode = "wb", quiet = TRUE)
}

soft <- readLines(gzfile(soft_gz))

recs <- list(); cur <- NULL
for (ln in soft) {
  if (startsWith(ln, "^SAMPLE")) {
    if (!is.null(cur)) recs[[length(recs) + 1]] <- cur
    cur <- list(gsm = trimws(sub("^\\^SAMPLE\\s*=\\s*", "", ln)),
                title = NA_character_, source = NA_character_, fpkm = NA_character_)
  } else if (!is.null(cur)) {
    if (startsWith(ln, "!Sample_title"))
      cur$title <- trimws(sub("^[^=]+=\\s*", "", ln))
    else if (startsWith(ln, "!Sample_source_name_ch1"))
      cur$source <- trimws(sub("^[^=]+=\\s*", "", ln))
    else if (startsWith(ln, "!Sample_supplementary_file") && grepl("fpkm_tracking", ln))
      cur$fpkm <- trimws(sub("^[^=]+=\\s*", "", ln))
  }
}
if (!is.null(cur)) recs[[length(recs) + 1]] <- cur

meta <- rbindlist(lapply(recs, as.data.table))
msg("GSE58135: ", nrow(meta), " samples in SOFT record")

meta[, group := fifelse(grepl("Reduction Mammoplasty", source), "healthy_normal",
                fifelse(grepl("^Uninvolved", source),          "adjacent_normal",
                fifelse(grepl("Cell Line", source),            "cell_line",
                fifelse(grepl("Primary Tumor", source),        "tumor",
                        NA_character_))))]
if (anyNA(meta$group)) stop("unclassified sample source: ",
                            paste(unique(meta$source[is.na(meta$group)]), collapse = " | "))
meta[, subtype := fifelse(grepl("Triple Negative|TNBC", source), "TNBC",
                  fifelse(grepl("ER\\+", source), "ERpos", NA_character_))]

print(meta[, .N, by = .(group, subtype)])

keep <- meta[group != "cell_line"]
msg("GSE58135: keeping ", nrow(keep), " tissue samples (dropped ",
    nrow(meta) - nrow(keep), " cell lines)")

fwrite(keep, file.path(DIR_DATA, "GSE58135_sample_metadata.csv"))

supp_dir <- file.path(DIR_DATA, "GSE58135_fpkm")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(keep))) {
  url  <- keep$fpkm[i]
  if (is.na(url)) { warning("no FPKM url for ", keep$gsm[i]); next }
  dest <- file.path(supp_dir, basename(url))
  if (file.exists(dest) && file.size(dest) > 1000) next
  url_https <- sub("^ftp://ftp\\.ncbi\\.nlm\\.nih\\.gov", "https://ftp.ncbi.nlm.nih.gov", url)
  tryCatch({
    download.file(url_https, dest, mode = "wb", quiet = TRUE)
  }, error = function(e) warning("failed ", keep$gsm[i], ": ", conditionMessage(e)))
  if (i %% 20 == 0) msg("GSE58135: ", i, "/", nrow(keep), " files")
}

got <- list.files(supp_dir, pattern = "fpkm_tracking")
msg("GSE58135: ", length(got), "/", nrow(keep), " FPKM files present")
msg("done.")
