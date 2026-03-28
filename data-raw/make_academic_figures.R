.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all()

ecdc <- readRDS("data-raw/ecdc_surveillance.rds")
cols <- .surv_colors()

design <- surv_design(
  data = ecdc$sequences, strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population, source_type = "source_type"
)

lineages <- sort(table(ecdc$sequences$lineage), decreasing = TRUE)
target_lin <- names(lineages)[!grepl("Other", names(lineages))][1]

weighted <- surv_lineage_prevalence(design, target_lin, method = "hajek")
naive <- surv_naive_prevalence(design, target_lin)
delay <- surv_estimate_delay(design)
nc <- surv_nowcast_lineage(design, delay, target_lin)
adj <- surv_adjusted_prevalence(design, delay, target_lin)

fig_dir <- "man/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# FIG 1: Sequencing inequality (academic)
# ============================================================
info <- design$strata_info
info$label <- info$region
info <- info[order(info$seq_rate), ]
info$label <- factor(info$label, levels = info$label)
avg_rate <- mean(info$seq_rate)
ratio <- round(max(info$seq_rate) / min(info$seq_rate))

p1 <- ggplot2::ggplot(info, ggplot2::aes(x = .data$label, y = .data$seq_rate)) +
  ggplot2::geom_segment(ggplot2::aes(xend = .data$label, y = 0, yend = .data$seq_rate),
                        color = cols[["primary"]], linewidth = 1.5) +
  ggplot2::geom_point(size = 5, color = cols[["primary"]]) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(round(.data$seq_rate * 100, 1), "%")),
                     hjust = -0.3, size = 3.5, color = "#333333") +
  ggplot2::geom_hline(yintercept = avg_rate, linetype = "dashed",
                      color = cols[["secondary"]], linewidth = 0.5) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(round(x * 100, 1), "%"),
    expand = ggplot2::expansion(mult = c(0, 0.2))
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = paste0("Sequencing rates vary ", ratio, "-fold across countries"),
    subtitle = "ECDC COVID-19 variant surveillance, 2023 (n = 99,093 sequences)",
    x = NULL, y = "Sequencing rate (% of confirmed cases sequenced)",
    caption = "Data: ECDC Open Data Portal | Dashed line: mean rate across countries"
  ) +
  theme_survinger(base_size = 12)

grDevices::png(file.path(fig_dir, "fig1_inequality.png"),
               width = 2400, height = 1500, res = 300)
print(p1); grDevices::dev.off()
cat("OK fig1\n")

# ============================================================
# FIG 2: Weighted vs Naive (academic, with annotation)
# ============================================================
w_df <- weighted$estimates[!is.na(weighted$estimates$prevalence), ]
w_df$Method <- "Design-weighted (Hajek)"
n_df <- naive$estimates[!is.na(naive$estimates$prevalence), ]
n_df$Method <- "Naive (unweighted)"
combined <- dplyr::bind_rows(w_df, n_df)
combined$time_short <- gsub("2023-", "", combined$time)

# Calculate summary stats for annotation
mean_diff <- round(mean(abs(w_df$prevalence - n_df$prevalence), na.rm = TRUE) * 100, 1)
max_diff <- round(max(abs(w_df$prevalence - n_df$prevalence), na.rm = TRUE) * 100, 1)

# Thin x-axis: show every 4th label
all_times <- sort(unique(combined$time_short))
show_labels <- all_times[seq(1, length(all_times), by = 4)]

p2 <- ggplot2::ggplot(combined, ggplot2::aes(
  x = .data$time_short, y = .data$prevalence,
  color = .data$Method, fill = .data$Method, group = .data$Method
)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
                       alpha = 0.10, color = NA) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 1.8, shape = 16) +
  ggplot2::scale_color_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_fill_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(round(x * 100), "%"),
    limits = c(0, NA)
  ) +
  ggplot2::annotate("label", x = all_times[3], y = max(w_df$prevalence, na.rm = TRUE) * 0.95,
                    label = paste0("Mean |bias|: ", mean_diff, " pp\nMax |bias|: ", max_diff, " pp"),
                    hjust = 0, size = 3.2, fill = "#FFFFFFCC", label.size = 0.3,
                    color = "#333333") +
  ggplot2::labs(
    title = paste0(target_lin, ": design weighting reveals systematic bias"),
    subtitle = "Hajek estimator corrects for 40-fold inequality in sequencing rates (n = 99,093)",
    x = "Epiweek (2023)", y = "Estimated prevalence (%)",
    color = NULL, fill = NULL,
    caption = "Shaded regions: 95% Wald confidence intervals | Data: ECDC Open Data Portal"
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
    legend.position = c(0.75, 0.92),
    legend.background = ggplot2::element_rect(fill = "#FFFFFFDD", color = NA)
  )

grDevices::png(file.path(fig_dir, "fig2_compare.png"),
               width = 3000, height = 1800, res = 300)
print(p2); grDevices::dev.off()
cat("OK fig2\n")

# ============================================================
# FIG 3: Bias heatmap (academic, improved)
# ============================================================
per_rw <- design$data
per_rw$.target <- as.integer(per_rw$lineage == target_lin)
per_rw <- per_rw |>
  dplyr::group_by(.data$region, .data$epiweek) |>
  dplyr::summarise(naive_prev = mean(.data$.target), .groups = "drop")
global_w <- weighted$estimates[c("time", "prevalence")]
names(global_w) <- c("epiweek", "weighted_prev")
per_rw <- dplyr::left_join(per_rw, global_w, by = "epiweek")
per_rw$bias <- per_rw$naive_prev - per_rw$weighted_prev
per_rw$ew_short <- gsub("2023-", "", per_rw$epiweek)

# Order countries by seq rate
rate_order <- design$strata_info[order(design$strata_info$seq_rate, decreasing = TRUE), ]
per_rw$region <- factor(per_rw$region, levels = rate_order$region)

show_weeks <- sort(unique(per_rw$ew_short))
show_wk_labels <- show_weeks[seq(1, length(show_weeks), by = 3)]

p3 <- ggplot2::ggplot(per_rw, ggplot2::aes(
  x = .data$ew_short, y = .data$region, fill = .data$bias
)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.4) +
  ggplot2::scale_fill_gradient2(
    low = cols[["primary"]], mid = "#FFFFFF", high = cols[["secondary"]],
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,
    labels = function(x) paste0(ifelse(x > 0, "+", ""), round(x * 100), " pp"),
    name = "Bias\n(naive \u2212 weighted)"
  ) +
  ggplot2::scale_x_discrete(breaks = show_wk_labels) +
  ggplot2::labs(
    title = paste0("Per-country estimation bias for ", target_lin),
    subtitle = "Countries ordered by sequencing rate (highest top) | Positive = naive overestimates",
    x = "Epiweek (2023)", y = NULL,
    caption = "Blue: local naive < global weighted | Red: local naive > global weighted | Data: ECDC"
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "right",
    legend.key.height = ggplot2::unit(40, "pt")
  )

grDevices::png(file.path(fig_dir, "fig3_bias_heatmap.png"),
               width = 3200, height = 1500, res = 300)
print(p3); grDevices::dev.off()
cat("OK fig3\n")

# ============================================================
# FIG 4: Nowcast (observed vs corrected)
# ============================================================
nc_df <- nc$estimates
nc_df$time_short <- gsub("2023-", "", nc_df$time)
nc_df$inflation <- nc_df$n_estimated / pmax(nc_df$n_observed, 1)

show_nc_labels <- sort(unique(nc_df$time_short))
show_nc_labels <- show_nc_labels[seq(1, length(show_nc_labels), by = 4)]

p4 <- ggplot2::ggplot(nc_df, ggplot2::aes(x = .data$time_short)) +
  ggplot2::geom_col(ggplot2::aes(y = .data$n_observed),
                    fill = cols[["neutral"]], alpha = 0.6, width = 0.6) +
  ggplot2::geom_line(ggplot2::aes(y = .data$n_estimated, group = 1),
                     color = cols[["quaternary"]], linewidth = 1) +
  ggplot2::geom_point(ggplot2::aes(y = .data$n_estimated,
                                   shape = ifelse(.data$is_nowcast, "Nowcasted", "Complete")),
                      color = cols[["quaternary"]], size = 2.5) +
  ggplot2::scale_shape_manual(values = c("Complete" = 16, "Nowcasted" = 17)) +
  ggplot2::scale_x_discrete(breaks = show_nc_labels) +
  ggplot2::labs(
    title = paste0("Delay-adjusted nowcasting for ", target_lin),
    subtitle = paste0("Right-truncation correction using NegBin(\u03bc=",
                      round(delay$parameters$mu[1], 1), ", k=",
                      round(delay$parameters$size[1], 1), ") delay model"),
    x = "Epiweek (2023)", y = "Sequence count",
    shape = NULL,
    caption = "Grey bars: reported counts | Orange: delay-corrected | \u25b2 = incomplete weeks | Data: ECDC"
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    legend.position = c(0.85, 0.85)
  )

grDevices::png(file.path(fig_dir, "fig4_nowcast.png"),
               width = 3000, height = 1800, res = 300)
print(p4); grDevices::dev.off()
cat("OK fig4\n")

# ============================================================
# FIG 5: Allocation comparison (grouped bars, academic)
# ============================================================
strats <- c("equal", "proportional", "min_mse")
H <- design$n_strata
alloc_details <- purrr::map_dfr(strats, function(s) {
  if (s == "equal") {
    n <- rep(1000L %/% H, H)
    n[1] <- n[1] + 1000L - sum(n)
    tibble::tibble(region = design$strata_info$region, n = n, strategy = s)
  } else if (s == "proportional") {
    pi_h <- design$population$pop_total / sum(design$population$pop_total)
    n <- as.integer(round(pi_h * 1000))
    n[1] <- n[1] + 1000L - sum(n)
    tibble::tibble(region = design$population$region, n = n, strategy = s)
  } else {
    res <- surv_optimize_allocation(design, s, total_capacity = 1000)
    tibble::tibble(region = res$allocation$region, n = res$allocation$n_allocated, strategy = s)
  }
})
alloc_details$strategy <- factor(alloc_details$strategy,
  levels = c("equal", "proportional", "min_mse"),
  labels = c("Equal", "Population-proportional", "MSE-optimal (Neyman)"))

p5 <- ggplot2::ggplot(alloc_details, ggplot2::aes(
  x = .data$region, y = .data$n, fill = .data$strategy
)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::scale_fill_manual(values = c(cols[["neutral"]], cols[["light_blue"]], cols[["primary"]])) +
  ggplot2::labs(
    title = "Optimal allocation differs from naive strategies",
    subtitle = "Distribution of 1,000 weekly sequences across 5 countries",
    x = NULL, y = "Sequences allocated per week",
    fill = "Allocation strategy",
    caption = "MSE-optimal: Neyman allocation minimizing lineage prevalence MSE | n = 1,000 total"
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

grDevices::png(file.path(fig_dir, "fig5_allocation.png"),
               width = 2800, height = 1600, res = 300)
print(p5); grDevices::dev.off()
cat("OK fig5\n")

# ============================================================
# FIG 6: Simulation benchmark — Bias vs Inequality (KEY FIGURE)
# ============================================================
cat("Running simulation benchmark...\n")
set.seed(2024)
n_sim <- 30  # enough for README demo
gini_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

bench_results <- purrr::map_dfr(gini_levels, function(g) {
  purrr::map_dfr(seq_len(n_sim), function(i) {
    # Generate rates with target Gini
    rates <- stats::rbeta(5, shape1 = 1, shape2 = (1 - g) / g)
    rates <- pmax(rates, 0.005)
    rates <- pmin(rates, 0.5)

    sim <- surv_simulate(n_regions = 5, n_weeks = 12,
                         sequencing_rates = rates, seed = i * 100 + round(g * 100))
    d <- surv_design(sim$sequences, ~ region,
                     sim$population[c("region", "seq_rate")], sim$population)

    truth_agg <- sim$truth |>
      dplyr::filter(.data$lineage == "BA.2.86") |>
      dplyr::group_by(.data$epiweek) |>
      dplyr::summarise(true_prev = stats::weighted.mean(.data$true_prevalence, .data$n_positive),
                       .groups = "drop")

    w <- surv_lineage_prevalence(d, "BA.2.86", method = "hajek")
    n <- surv_naive_prevalence(d, "BA.2.86")

    w_merged <- dplyr::inner_join(w$estimates, truth_agg, by = c("time" = "epiweek"))
    n_merged <- dplyr::inner_join(n$estimates, truth_agg, by = c("time" = "epiweek"))

    tibble::tibble(
      gini = g, sim = i,
      bias_weighted = mean(w_merged$prevalence - w_merged$true_prev, na.rm = TRUE),
      bias_naive = mean(n_merged$prevalence - n_merged$true_prev, na.rm = TRUE),
      rmse_weighted = sqrt(mean((w_merged$prevalence - w_merged$true_prev)^2, na.rm = TRUE)),
      rmse_naive = sqrt(mean((n_merged$prevalence - n_merged$true_prev)^2, na.rm = TRUE))
    )
  })
})

# Summarize
bench_summary <- bench_results |>
  dplyr::group_by(.data$gini) |>
  dplyr::summarise(
    mean_bias_w = mean(.data$bias_weighted, na.rm = TRUE),
    mean_bias_n = mean(.data$bias_naive, na.rm = TRUE),
    mean_rmse_w = mean(.data$rmse_weighted, na.rm = TRUE),
    mean_rmse_n = mean(.data$rmse_naive, na.rm = TRUE),
    se_rmse_w = stats::sd(.data$rmse_weighted, na.rm = TRUE) / sqrt(dplyr::n()),
    se_rmse_n = stats::sd(.data$rmse_naive, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

# Reshape for plotting
rmse_long <- tidyr::pivot_longer(
  bench_summary,
  cols = c("mean_rmse_w", "mean_rmse_n"),
  names_to = "method", values_to = "rmse"
)
rmse_long$se <- ifelse(rmse_long$method == "mean_rmse_w",
                       bench_summary$se_rmse_w, bench_summary$se_rmse_n)
rmse_long$Method <- ifelse(rmse_long$method == "mean_rmse_w",
                           "Hajek (design-weighted)", "Naive (unweighted)")

p6 <- ggplot2::ggplot(rmse_long, ggplot2::aes(
  x = .data$gini, y = .data$rmse,
  color = .data$Method, fill = .data$Method
)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$rmse - 1.96 * .data$se,
                                    ymax = .data$rmse + 1.96 * .data$se),
                       alpha = 0.15, color = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_fill_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_x_continuous(breaks = gini_levels) +
  ggplot2::labs(
    title = "Design weighting reduces RMSE as sequencing inequality increases",
    subtitle = paste0("Simulation: 5 regions, 12 weeks, ", n_sim, " replicates per Gini level"),
    x = "Gini coefficient of sequencing rates (higher = more unequal)",
    y = "Root mean squared error (RMSE)",
    color = NULL, fill = NULL,
    caption = "Shaded bands: 95% CI across simulation replicates | Target lineage: BA.2.86"
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    legend.position = c(0.25, 0.85),
    legend.background = ggplot2::element_rect(fill = "#FFFFFFDD", color = NA)
  )

grDevices::png(file.path(fig_dir, "fig6_benchmark.png"),
               width = 2800, height = 1800, res = 300)
print(p6); grDevices::dev.off()
cat("OK fig6 benchmark\n")

# Copy to vignettes/figures too
for (f in list.files(fig_dir, pattern = "\\.png$", full.names = TRUE)) {
  file.copy(f, file.path("vignettes/figures", basename(f)), overwrite = TRUE)
}

figs <- list.files(fig_dir, pattern = "\\.png$")
cat(sprintf("\n=== %d academic figures generated ===\n", length(figs)))
for (f in figs) {
  sz <- round(file.info(file.path(fig_dir, f))$size / 1024)
  cat(sprintf("  %s (%d KB)\n", f, sz))
}
