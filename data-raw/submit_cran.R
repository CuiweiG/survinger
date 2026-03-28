.libPaths(Sys.getenv("R_LIBS_USER"))

pkg_path <- "C:/Users/openclaw/survinger_0.1.0.tar.gz"
cat("Package:", pkg_path, "\n")
cat("Exists:", file.exists(pkg_path), "\n")
cat("Size:", round(file.info(pkg_path)$size / 1024^2, 2), "MB\n\n")

comment_text <- paste(
  "New submission. survinger provides design-adjusted inference for",
  "pathogen genomic surveillance under unequal sequencing rates and",
  "reporting delays. Implements Horvitz-Thompson/Hajek/post-stratified",
  "estimators with Wilson score confidence intervals, constrained Neyman",
  "allocation optimization, and right-truncation-corrected nowcasting.",
  "Validated against survey::svymean (exact match) and on real ECDC data",
  "(99,093 sequences, 5 EU countries).",
  "R CMD check --as-cran: 0 errors, 0 warnings, 0 notes.",
  "Fills the methodological gap between phylosamp (sample size calculation)",
  "and epinowcast (Bayesian nowcasting) with a lightweight, pure-R solution."
)

cat("Uploading to CRAN...\n")
resp <- httr::POST(
  "https://xmpalantir.wu.ac.at/cransubmit/index2.php",
  body = list(
    name = "CuiweiG",
    email = "48gaocuiwei@gmail.com",
    uploaded_file = httr::upload_file(pkg_path, type = "application/x-gzip"),
    comment = comment_text,
    upload = "Upload the package"
  ),
  encode = "multipart"
)

cat("HTTP status:", httr::status_code(resp), "\n")
body <- httr::content(resp, "text", encoding = "UTF-8")

# Save full response
writeLines(body, "C:/Users/openclaw/cran_step1.html")
cat("Response saved to cran_step1.html\n")

# Check what step we're at
if (grepl("Step 2", body)) {
  cat("\n=== UPLOAD SUCCESSFUL — Moved to Step 2 ===\n")
  cat("CRAN will send a confirmation email to 48gaocuiwei@gmail.com\n")

  # Try to extract pkg_id and auto-confirm
  id_matches <- regmatches(body, gregexpr('value="[0-9a-zA-Z_-]+"', body))
  cat("Hidden values found:", paste(unlist(id_matches), collapse = ", "), "\n")

  # Extract the pkg_id
  pkg_id_line <- regmatches(body, regexpr('name="pkg_id"[^>]+value="[^"]+"', body))
  if (length(pkg_id_line) > 0) {
    pkg_id <- sub('.*value="([^"]+)"', "\\1", pkg_id_line)
    cat("pkg_id:", pkg_id, "\n")

    # Step 2: Submit confirmation
    cat("\nSubmitting Step 2 confirmation...\n")
    resp2 <- httr::POST(
      "https://xmpalantir.wu.ac.at/cransubmit/index2.php",
      body = list(
        pkg_id = pkg_id,
        name = "CuiweiG",
        email = "48gaocuiwei@gmail.com",
        policy_check = "1",
        submit = "Submit package"
      ),
      encode = "multipart"
    )
    body2 <- httr::content(resp2, "text", encoding = "UTF-8")
    writeLines(body2, "C:/Users/openclaw/cran_step2.html")

    if (grepl("Step 3", body2) || grepl("confirmation", body2, ignore.case = TRUE)) {
      cat("\n=== STEP 2 COMPLETE ===\n")
      cat("Check your email at 48gaocuiwei@gmail.com for confirmation link.\n")
      cat("Click the link to finalize CRAN submission.\n")
    } else {
      cat("Step 2 response saved to cran_step2.html\n")
    }
  }
} else if (grepl("error|Error", body, ignore.case = TRUE)) {
  errors <- regmatches(body, gregexpr("<font[^>]*color[^>]*red[^>]*>[^<]+</font>", body))
  cat("ERRORS:", paste(unlist(errors), collapse = "\n"), "\n")
} else {
  cat("Unexpected response. Check cran_step1.html\n")
}
