# Reproduction script for sarscov2_surveillance dataset
devtools::load_all()
sarscov2_surveillance <- surv_simulate(
  n_regions = 5, n_weeks = 26,
  sequencing_rates = c(0.15, 0.08, 0.03, 0.01, 0.005),
  delay_params = list(mu = 10, size = 3),
  sources = c("clinical", "wastewater", "sentinel"),
  source_weights = c(0.7, 0.2, 0.1),
  seed = 20240101
)
usethis::use_data(sarscov2_surveillance, overwrite = TRUE)
