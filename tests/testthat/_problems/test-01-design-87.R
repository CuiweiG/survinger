# Extracted from test-01-design.R:87

# test -------------------------------------------------------------------------
sim <- surv_simulate(n_regions = 2, n_weeks = 2, seed = 20)
expect_error(
    surv_design(
      data = sim$sequences, strata = ~ nonexistent,
      sequencing_rate = sim$population[c("region", "seq_rate")],
      population = sim$population
    ), "not found"
  )
