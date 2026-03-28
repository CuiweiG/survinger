.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all()

cols <- .surv_colors()
fig_dir <- "man/figures"

# ============================================================
# FIG 6 FIX: Benchmark with REALISTIC heterogeneous prevalence
#
# The key insight: naive is biased when high-sequencing regions
# have DIFFERENT prevalence than low-sequencing regions.
# We simulate this by correlating prevalence with sequencing rate.
# ============================================================

cat("Running corrected simulation benchmark...\n")
set.seed(2024)
n_sim <- 50
gini_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

bench_results <- purrr::map_dfr(gini_levels, function(g) {
  purrr::map_dfr(seq_len(n_sim), function(i) {
    # Generate rates with target Gini
    rates <- stats::rbeta(5, shape1 = 1, shape2 = max(0.5, (1 - g) / g))
    rates <- pmax(rates, 0.005)
    rates <- pmin(rates, 0.5)
    rates <- sort(rates, decreasing = TRUE)  # Region_A = highest rate

    # KEY FIX: Create heterogeneous prevalence that correlates with seq rate
    # High-sequencing regions see MORE of the target variant (realistic:
    # sentinel sites in urban areas detect emerging variants earlier)
    sim <- surv_simulate(n_regions = 5, n_weeks = 12,
                         sequencing_rates = rates,
                         seed = i * 1000 + round(g * 100))

    # Modify truth: make prevalence systematically different by region
    # High-seq regions get higher BA.2.86 prevalence
    dat <- sim$sequences
    n_total <- nrow(dat)
    regions <- sort(unique(dat$region))

    # For each region, resample lineages with region-specific probabilities
    # Higher-seq regions get more BA.2.86
    prev_by_region <- seq(0.3, 0.05, length.out = 5)  # Region_A: 30%, Region_E: 5%
    names(prev_by_region) <- regions

    for (r in regions) {
      idx <- which(dat$region == r)
      p_target <- prev_by_region[r]
      n_r <- length(idx)
      new_lineages <- sample(c("BA.2.86", "Other"), n_r, replace = TRUE,
                             prob = c(p_target, 1 - p_target))
      dat$lineage[idx] <- new_lineages
    }

    # True population prevalence (weighted by population, not sequencing)
    pop_shares <- sim$population$pop_total / sum(sim$population$pop_total)
    true_prev <- sum(pop_shares * prev_by_region[sim$population$region])

    # Build design
    sim$sequences <- dat
    d <- tryCatch(
      surv_design(dat, ~ region,
                  sim$population[c("region", "seq_rate")], sim$population),
      error = function(e) NULL
    )
    if (is.null(d)) return(NULL)

    # Weighted estimate (pooled across all time)
    w_dat <- .map_weights_to_obs(d)
    w_dat$.target <- as.integer(w_dat$lineage == "BA.2.86")
    est_weighted <- sum(w_dat$weight * w_dat$.target) / sum(w_dat$weight)
    est_naive <- mean(w_dat$.target)

    tibble::tibble(
      gini = g, sim = i,
      bias_weighted = est_weighted - true_prev,
      bias_naive = est_naive - true_prev,
      abs_bias_weighted = abs(est_weighted - true_prev),
      abs_bias_naive = abs(est_naive - true_prev)
    )
  })
})

# Summarize
bench_summary <- bench_results |>
  dplyr::group_by(.data$gini) |>
  dplyr::summarise(
    mean_absbias_w = mean(.data$abs_bias_weighted, na.rm = TRUE),
    mean_absbias_n = mean(.data$abs_bias_naive, na.rm = TRUE),
    se_absbias_w = stats::sd(.data$abs_bias_weighted, na.rm = TRUE) / sqrt(dplyr::n()),
    se_absbias_n = stats::sd(.data$abs_bias_naive, na.rm = TRUE) / sqrt(dplyr::n()),
    mean_bias_w = mean(.data$bias_weighted, na.rm = TRUE),
    mean_bias_n = mean(.data$bias_naive, na.rm = TRUE),
    .groups = "drop"
  )

cat("Benchmark summary:\n")
print(bench_summary)

# Plot: Mean Absolute Bias vs Gini
bias_long <- tidyr::pivot_longer(
  bench_summary,
  cols = c("mean_absbias_w", "mean_absbias_n"),
  names_to = "method", values_to = "abs_bias"
)
bias_long$se <- ifelse(bias_long$method == "mean_absbias_w",
                       bench_summary$se_absbias_w, bench_summary$se_absbias_n)
bias_long$Method <- ifelse(bias_long$method == "mean_absbias_w",
                           "Hajek (design-weighted)", "Naive (unweighted)")

p6 <- ggplot2::ggplot(bias_long, ggplot2::aes(
  x = .data$gini, y = .data$abs_bias,
  color = .data$Method, fill = .data$Method
)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(0, .data$abs_bias - 1.96 * .data$se),
                                    ymax = .data$abs_bias + 1.96 * .data$se),
                       alpha = 0.15, color = NA) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 3.5) +
  ggplot2::scale_color_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_fill_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_x_continuous(breaks = gini_levels) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100, 1), " pp")) +
  ggplot2::labs(
    title = "Design weighting reduces bias under heterogeneous prevalence",
    subtitle = paste0("Simulation: 5 regions, heterogeneous true prevalence (5%\u201330%), ",
                      n_sim, " replicates per Gini level"),
    x = "Gini coefficient of sequencing rates (higher = more unequal)",
    y = "Mean absolute bias (percentage points)",
    color = NULL, fill = NULL,
    caption = "Prevalence correlated with sequencing rate (realistic scenario) | Shaded: 95% CI"
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    legend.position = c(0.70, 0.88),
    legend.background = ggplot2::element_rect(fill = "#FFFFFFDD", color = NA)
  )

grDevices::png(file.path(fig_dir, "fig6_benchmark.png"),
               width = 2800, height = 1800, res = 300)
print(p6); grDevices::dev.off()
cat("OK fig6_benchmark (corrected)\n")

# ============================================================
# FIG 4 FIX: Nowcast with actual right-truncation
#
# Truncate data to simulate real-world incomplete reporting
# ============================================================

ecdc <- readRDS("data-raw/ecdc_surveillance.rds")
design <- surv_design(
  data = ecdc$sequences, strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population, source_type = "source_type"
)

lineages <- sort(table(ecdc$sequences$lineage), decreasing = TRUE)
target_lin <- names(lineages)[!grepl("Other", names(lineages))][1]

delay <- surv_estimate_delay(design)

# Create a truncated version: pretend the reference date is 3 weeks earlier
# This simulates looking at the data when recent weeks are incomplete
ref_date <- max(ecdc$sequences$report_date, na.rm = TRUE)
truncation_date <- ref_date - 21  # 3 weeks back

# Filter: only include sequences reported by truncation_date
trunc_data <- ecdc$sequences[ecdc$sequences$report_date <= truncation_date, ]
cat("Original:", nrow(ecdc$sequences), "-> Truncated:", nrow(trunc_data), "\n")

design_trunc <- surv_design(
  data = trunc_data, strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population, source_type = "source_type"
)

nc_trunc <- surv_nowcast_lineage(design_trunc, delay, target_lin,
                                 ref_date = truncation_date)

nc_df <- nc_trunc$estimates
nc_df$time_short <- gsub("2023-", "", nc_df$time)

# Also get the "full" counts for comparison
full_counts <- ecdc$sequences[ecdc$sequences$lineage == target_lin, ] |>
  dplyr::mutate(ew_short = gsub("2023-", "", .data$epiweek)) |>
  dplyr::group_by(.data$ew_short) |>
  dplyr::summarise(n_full = dplyr::n(), .groups = "drop")

nc_df <- dplyr::left_join(nc_df, full_counts, by = c("time_short" = "ew_short"))

show_labels <- sort(unique(nc_df$time_short))
show_labels <- show_labels[seq(1, length(show_labels), by = 4)]

p4 <- ggplot2::ggplot(nc_df, ggplot2::aes(x = .data$time_short)) +
  # Full data (what we'd see later)
  ggplot2::geom_col(ggplot2::aes(y = .data$n_full),
                    fill = "#E8E8E8", width = 0.7) +
  # Truncated observed
  ggplot2::geom_col(ggplot2::aes(y = .data$n_observed),
                    fill = cols[["neutral"]], alpha = 0.7, width = 0.5) +
  # Nowcast estimate
  ggplot2::geom_line(ggplot2::aes(y = .data$n_estimated, group = 1),
                     color = cols[["quaternary"]], linewidth = 1.1) +
  ggplot2::geom_point(ggplot2::aes(y = .data$n_estimated,
                                   shape = ifelse(.data$is_nowcast, "Nowcasted", "Complete")),
                      color = cols[["quaternary"]], size = 2.5) +
  ggplot2::scale_shape_manual(values = c("Complete" = 16, "Nowcasted" = 17)) +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::labs(
    title = paste0("Nowcasting recovers right-truncated counts for ", target_lin),
    subtitle = paste0("Simulated reporting cutoff: 3 weeks before data end | ",
                      "NegBin(\u03bc=", round(delay$parameters$mu[1], 1),
                      ", k=", round(delay$parameters$size[1], 1), ")"),
    x = "Epiweek (2023)", y = "Sequence count",
    shape = NULL,
    caption = paste0("Light grey: eventual full count | Dark grey: reported by cutoff | ",
                     "Orange: delay-corrected nowcast | \u25b2 = incomplete | Data: ECDC")
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    legend.position = c(0.85, 0.85)
  )

grDevices::png(file.path(fig_dir, "fig4_nowcast.png"),
               width = 3000, height = 1800, res = 300)
print(p4); grDevices::dev.off()
cat("OK fig4_nowcast (corrected with truncation)\n")

# Copy to vignettes
file.copy(file.path(fig_dir, "fig6_benchmark.png"),
          "vignettes/figures/fig6_benchmark.png", overwrite = TRUE)
file.copy(file.path(fig_dir, "fig4_nowcast.png"),
          "vignettes/figures/fig4_nowcast.png", overwrite = TRUE)

cat("\n=== CORRECTED FIGURES DONE ===\n")
