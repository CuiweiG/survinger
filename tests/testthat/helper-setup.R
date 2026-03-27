# Shared test fixtures - auto-loaded by testthat

make_test_data <- function(n_regions = 3, n_weeks = 4, seed = 42) {
  set.seed(seed)
  regions <- paste0("Region_", LETTERS[seq_len(n_regions)])
  lineages <- c("BA.5", "XBB.1.5", "BA.2.86", "Other")
  sources <- c("clinical", "wastewater", "sentinel")
  records <- vector("list", n_regions * n_weeks)
  idx <- 0L
  for (w in seq_len(n_weeks)) {
    week_start <- as.Date("2024-01-01") + (w - 1L) * 7L
    for (r in seq_len(n_regions)) {
      idx <- idx + 1L
      n <- sample(20:50, 1L)
      records[[idx]] <- tibble::tibble(
        sequence_id = paste0("s_", w, "_", r, "_", seq_len(n)),
        region = regions[r],
        source_type = sample(sources, n, replace = TRUE, prob = c(.7, .2, .1)),
        lineage = sample(lineages, n, replace = TRUE, prob = c(.4, .3, .2, .1)),
        collection_date = week_start + sample(0:6, n, replace = TRUE),
        report_date = week_start + sample(0:6, n, replace = TRUE) +
          stats::rnbinom(n, mu = 10, size = 3)
      )
    }
  }
  sequences <- dplyr::bind_rows(records)
  seq_rates <- c(0.10, 0.02, 0.005)[seq_len(n_regions)]
  n_pos <- c(500, 300, 200)[seq_len(n_regions)] * n_weeks
  population <- tibble::tibble(
    region = regions, n_positive = n_pos,
    n_sequenced = as.integer(n_pos * seq_rates),
    seq_rate = seq_rates,
    pop_total = c(50000, 30000, 20000)[seq_len(n_regions)]
  )
  list(sequences = sequences, population = population)
}
