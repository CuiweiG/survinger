.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all(quiet = TRUE)

cols <- .surv_colors()
fig_dir <- "man/figures"

ecdc <- readRDS("data-raw/ecdc_surveillance.rds")
d <- surv_design(
  data = ecdc$sequences, strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population, source_type = "source_type"
)

H <- d$n_strata
strats <- c("equal", "proportional", "min_mse")
ad <- purrr::map_dfr(strats, function(s) {
  if (s == "equal") {
    n <- rep(1000L %/% H, H); n[1] <- n[1] + 1000L - sum(n)
    tibble::tibble(region = d$strata_info$region, n = n, strategy = s)
  } else if (s == "proportional") {
    pi <- d$population$pop_total / sum(d$population$pop_total)
    n <- as.integer(round(pi * 1000)); n[1] <- n[1] + 1000L - sum(n)
    tibble::tibble(region = d$population$region, n = n, strategy = s)
  } else {
    r <- surv_optimize_allocation(d, s, total_capacity = 1000)
    tibble::tibble(region = r$allocation$region, n = r$allocation$n_allocated, strategy = s)
  }
})
ad$strategy <- factor(ad$strategy, levels = c("equal","proportional","min_mse"),
                      labels = c("Equal","Population-proportional","MSE-optimal (Neyman)"))

cat("Allocation by strategy:\n")
print(tidyr::pivot_wider(ad, names_from = strategy, values_from = n))

p5 <- ggplot2::ggplot(ad, ggplot2::aes(x = .data$region, y = .data$n, fill = .data$strategy)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(0.8), width = 0.7) +
  ggplot2::scale_fill_manual(values = c(cols[["neutral"]], cols[["light_blue"]], cols[["primary"]])) +
  ggplot2::labs(
    title = "Optimal allocation redistributes sequences across countries",
    subtitle = "1,000 weekly sequences across 5 EU countries | ECDC data",
    x = NULL, y = "Sequences allocated per week", fill = "Strategy",
    caption = paste0("MSE-optimal: Neyman allocation for lineage prevalence | ",
                     "Data: ECDC (n = ", formatC(nrow(ecdc$sequences), big.mark=","), ")")
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

grDevices::png(file.path(fig_dir, "fig5_allocation.png"), width = 2800, height = 1600, res = 300)
print(p5); grDevices::dev.off()
file.copy(file.path(fig_dir, "fig5_allocation.png"),
          "vignettes/figures/fig5_allocation.png", overwrite = TRUE)
cat("=== fig5: ECDC 5-country allocation ===\n")
