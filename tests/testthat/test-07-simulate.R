test_that("surv_simulate produces expected structure", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 1)
  expect_type(sim, "list")
  expect_named(sim, c("sequences", "population", "truth", "parameters"))
  expect_s3_class(sim$sequences, "tbl_df")
  expect_true(all(c("sequence_id", "region", "source_type", "lineage",
                     "collection_date", "report_date", "epiweek") %in%
                    names(sim$sequences)))
  expect_true(nrow(sim$sequences) > 0)
  expect_s3_class(sim$population, "tbl_df")
  expect_equal(nrow(sim$population), 3)
  expect_true(all(sim$population$seq_rate > 0))
  expect_true(all(sim$population$seq_rate <= 1))
})

test_that("surv_simulate is reproducible with seed", {
  sim1 <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 42)
  sim2 <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 42)
  expect_identical(sim1$sequences, sim2$sequences)
})

test_that("surv_simulate uses custom sequencing rates", {
  rates <- c(0.5, 0.01, 0.001)
  sim <- surv_simulate(n_regions = 3, n_weeks = 10,
                       sequencing_rates = rates, seed = 1)
  by_region <- table(sim$sequences$region)
  expect_true(by_region["Region_A"] > by_region["Region_C"] * 3)
})

test_that("simulated delays are non-negative", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
  delays <- as.integer(sim$sequences$report_date - sim$sequences$collection_date)
  expect_true(all(delays >= 0))
})

test_that("true prevalence sums to ~1 per region-week", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 1)
  sums <- sim$truth |>
    dplyr::group_by(.data$epiweek, .data$region) |>
    dplyr::summarise(total = sum(.data$true_prevalence), .groups = "drop")
  expect_true(all(abs(sums$total - 1) < 0.01))
})
