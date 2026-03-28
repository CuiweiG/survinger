.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")

cat("=== VIGNETTES ===\n")
vigs <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
for (v in vigs) {
  content <- readLines(v, warn = FALSE)
  fname <- basename(v)
  has_index <- any(grepl("VignetteIndexEntry", content))
  has_engine <- any(grepl("VignetteEngine", content))
  has_enc <- any(grepl("VignetteEncoding", content))
  
  fig_refs <- regmatches(content, gregexpr("figures/[^\"')\\s]+\\.png", content))
  fig_refs <- unique(unlist(fig_refs))
  missing <- 0L
  for (fr in fig_refs) {
    if (!file.exists(file.path("vignettes", fr))) missing <- missing + 1L
  }
  
  cat(sprintf("%s: index=%s engine=%s enc=%s figs=%d missing=%d\n",
              fname, has_index, has_engine, has_enc, length(fig_refs), missing))
}

cat("\n=== S3 METHODS ===\n")
ns <- readLines("NAMESPACE")
s3 <- grep("^S3method", ns, value = TRUE)
cat("Registered:", length(s3), "\n")
for (m in s3) cat(" ", m, "\n")

cat("\n=== IMPORTS ACTUALLY USED ===\n")
r_code <- paste(unlist(lapply(list.files("R", full.names = TRUE, pattern = "\\.R$"),
                              readLines, warn = FALSE)), collapse = "\n")
declared_imports <- c("checkmate", "cli", "dplyr", "generics", "ggplot2",
                      "methods", "purrr", "rlang", "stats", "tibble")
for (pkg in declared_imports) {
  pattern <- paste0(pkg, "::|@importFrom\\s+", pkg)
  used <- grepl(pattern, r_code)
  cat(sprintf("  %s: %s\n", pkg, ifelse(used, "USED", "NOT USED")))
}

cat("\n=== FINAL R CMD CHECK ===\n")
devtools::document(quiet = TRUE)
chk <- devtools::check(error_on = "never", args = c("--as-cran", "--no-manual"), quiet = TRUE)
cat("E:", length(chk[["errors"]]), "\n")
cat("W:", length(chk[["warnings"]]), "\n")
cat("N:", length(chk[["notes"]]), "\n")
if (length(chk[["errors"]]) > 0) cat("ERRORS:\n", chk[["errors"]], "\n")
if (length(chk[["warnings"]]) > 0) cat("WARNINGS:\n", chk[["warnings"]], "\n")
if (length(chk[["notes"]]) > 0) cat("NOTES:\n", chk[["notes"]], "\n")

if (length(chk[["errors"]]) == 0 && length(chk[["warnings"]]) == 0 && length(chk[["notes"]]) == 0) {
  cat("\n=== PERFECT: 0E 0W 0N ===\n")
} else if (length(chk[["errors"]]) == 0 && length(chk[["warnings"]]) == 0) {
  cat("\n=== ACCEPTABLE: 0E 0W ===\n")
}
