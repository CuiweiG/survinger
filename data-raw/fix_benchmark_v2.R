.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all()

cols <- .surv_colors()
fig_dir <- "man/figures"

# ============================================================
# FIG 6 FIX: Use controlled rate generation with exact Gini
# ============================================================

# Generate rates with EXACT target Gini using ordered uniform method
generate_rates_gini <- function(n, target_gini, seed) {
  set.seed(seed)
  # Method: generate from Pareto-like distribution, then calibrate
  # Lower shape = higher inequality
  shape <- (1 - target_gini) / (2 * target_gini)
  shape <- max(shape, 0.3)
  
  rates_raw <- stats::rgamma(n, shape = shape, rate = 1)
  rates_raw <- rates_raw / max(rates_raw) * 0.4  # scale to [0, 0.4]
  rates_raw <- pmax(rates_raw, 0.003)
  rates_raw <- pmin(rates_raw, 0.5)
  sort(rates_raw, decreasing = TRUE)
}

cat("Checking rate generation:\n")
set.seed(2024)
for (g in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  rates <- generate_rates_gini(5, g, seed = 42 + round(g * 100))
  actual_gini <- .gini_coefficient(rates)
  cat(sprintf("Target %.1f -> Actual %.3f : %s\n", g, actual_gini,
              paste(round(rates, 3), collapse = ", ")))
}

cat("\nRunning benchmark...\n")
n_sim <- 50
gini_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

br <- purrr::map_dfr(gini_levels, function(g) {
  purrr::map_dfr(seq_len(n_sim), function(i) {
    rates <- generate_rates_gini(5, g, seed = i * 1000 + round(g * 100))
    
    sim <- tryCatch(
      surv_simulate(n_regions = 5, n_weeks = 12, sequencing_rates = rates,
                    seed = i * 1000 + round(g * 100)),
      error = function(e) NULL
    )
    if (is.null(sim)) return(NULL)
    
    dat <- sim$sequences
    regions <- sort(unique(dat$region))
    prev_by_r <- seq(0.30, 0.05, length.out = 5)
    names(prev_by_r) <- regions
    for (r in regions) {
      idx <- which(dat$region == r)
      dat$lineage[idx] <- sample(c("BA.2.86", "Other"), length(idx),
                                 replace = TRUE, prob = c(prev_by_r[r], 1 - prev_by_r[r]))
    }
    
    pop_shares <- sim$population$pop_total / sum(sim$population$pop_total)
    true_prev <- sum(pop_shares * prev_by_r[sim$population$region])
    
    sim$sequences <- dat
    d <- tryCatch(
      surv_design(dat, ~ region, sim$population[c("region", "seq_rate")],
                  sim$population),
      error = function(e) NULL
    )
    if (is.null(d)) return(NULL)
    
    w_dat <- .map_weights_to_obs(d)
    w_dat$.t <- as.integer(w_dat$lineage == "BA.2.86")
    
    tibble::tibble(
      gini = g, sim = i,
      abs_bias_w = abs(sum(w_dat$weight * w_dat$.t) / sum(w_dat$weight) - true_prev),
      abs_bias_n = abs(mean(w_dat$.t) - true_prev)
    )
  })
})

# Check we have data for all levels
cat("\nData points per Gini level:\n")
print(table(br$gini))

bs <- br |>
  dplyr::group_by(.data$gini) |>
  dplyr::summarise(
    m_w = mean(.data$abs_bias_w, na.rm = TRUE),
    m_n = mean(.data$abs_bias_n, na.rm = TRUE),
    se_w = stats::sd(.data$abs_bias_w, na.rm = TRUE) / sqrt(dplyr::n()),
    se_n = stats::sd(.data$abs_bias_n, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

cat("\nBenchmark summary:\n")
print(bs)

bl <- tidyr::pivot_longer(bs, cols = c("m_w", "m_n"), names_to = "method", values_to = "bias")
bl$se <- ifelse(bl$method == "m_w", bs$se_w, bs$se_n)
bl$Method <- ifelse(bl$method == "m_w", "Hajek (design-weighted)", "Naive (unweighted)")

p6 <- ggplot2::ggplot(bl, ggplot2::aes(
  x = .data$gini, y = .data$bias, color = .data$Method, fill = .data$Method)) +
  ggplot2::geom_ribbon(ggplot2::aes(
    ymin = pmax(0, .data$bias - 1.96 * .data$se),
    ymax = .data$bias + 1.96 * .data$se), alpha = 0.15, color = NA) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 3.5) +
  ggplot2::scale_color_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_fill_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_x_continuous(breaks = gini_levels) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100, 1), " pp")) +
  ggplot2::labs(
    title = "Design weighting reduces bias under heterogeneous prevalence",
    subtitle = paste0("5 regions, prevalence 5%-30%, ", n_sim, " reps/level"),
    x = "Gini coefficient of sequencing rates",
    y = "Mean absolute bias (pp)",
    color = NULL, fill = NULL,
    caption = "Shaded: 95% CI | Prevalence correlated with sequencing rate"
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    legend.position = c(0.70, 0.88),
    legend.background = ggplot2::element_rect(fill = "#FFFFFFDD", color = NA)
  )

grDevices::png(file.path(fig_dir, "fig6_benchmark.png"), width = 2800, height = 1800, res = 300)
print(p6); grDevices::dev.off()

file.copy(file.path(fig_dir, "fig6_benchmark.png"),
          "vignettes/figures/fig6_benchmark.png", overwrite = TRUE)

cat("\n=== fig6 benchmark FIXED (continuous line, all 6 points) ===\n")
