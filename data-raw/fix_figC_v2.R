.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all(quiet = TRUE)

cols <- .surv_colors()
fig_dir <- "man/figures"

# Simulate a scenario with a rare emerging variant (2% prevalence)
# Small sample sizes per week -> many p-hat = 0 weeks
set.seed(2024)
sim <- surv_simulate(
  n_regions = 5, n_weeks = 15,
  sequencing_rates = c(0.10, 0.05, 0.02, 0.005, 0.002),
  total_positive_per_week = 500,
  seed = 2024
)

d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)

# BA.2.86 is the rare emerging one in surv_simulate default dynamics
w <- surv_lineage_prevalence(d, "BA.2.86", method = "hajek")
est <- w$estimates[!is.na(w$estimates$prevalence), ]

# Compute Wald for comparison
z <- stats::qnorm(0.975)
est$wald_lower <- pmax(0, est$prevalence - z * est$se)
est$wald_upper <- pmin(1, est$prevalence + z * est$se)
est$wald_width <- est$wald_upper - est$wald_lower
est$wilson_width <- est$ci_upper - est$ci_lower

n_wald_zero <- sum(est$wald_width == 0)
n_wilson_zero <- sum(est$wilson_width == 0)
cat("Wald zero-width:", n_wald_zero, "/", nrow(est), "\n")
cat("Wilson zero-width:", n_wilson_zero, "/", nrow(est), "\n")
cat("Weeks with p-hat = 0:", sum(est$prevalence == 0), "\n")

est$wald_fails <- est$wald_width == 0 & est$wilson_width > 0

est_long <- dplyr::bind_rows(
  est |> dplyr::transmute(
    time = .data$time, prevalence = .data$prevalence,
    ci_lo = .data$wald_lower, ci_hi = .data$wald_upper,
    Method = "Wald (standard)", fails = .data$wald_fails
  ),
  est |> dplyr::transmute(
    time = .data$time, prevalence = .data$prevalence,
    ci_lo = .data$ci_lower, ci_hi = .data$ci_upper,
    Method = "Wilson (survinger)", fails = FALSE
  )
)
est_long$Method <- factor(est_long$Method,
  levels = c("Wald (standard)", "Wilson (survinger)"))

show_labels <- sort(unique(est$time))
show_labels <- show_labels[seq(1, length(show_labels), by = 2)]

pC <- ggplot2::ggplot(est_long, ggplot2::aes(x = .data$time, y = .data$prevalence * 100)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$ci_lo * 100, ymax = .data$ci_hi * 100),
    width = 0.3, linewidth = 0.7, color = cols[["primary"]]
  ) +
  ggplot2::geom_point(ggplot2::aes(color = .data$fails), size = 2.5) +
  ggplot2::scale_color_manual(values = c("FALSE" = cols[["primary"]], "TRUE" = "#CC0000"),
                              guide = "none") +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggplot2::facet_wrap(~ Method, ncol = 1, scales = "free_y") +
  ggplot2::labs(
    title = paste0("Wilson CI prevents zero-width intervals when p\u0302 = 0"),
    subtitle = paste0("Emerging variant BA.2.86 (simulated, ~2% true prevalence) | Wald: ",
                      n_wald_zero, "/", nrow(est),
                      " zero-width weeks (red) | Wilson: 0 zero-width"),
    x = "Epiweek", y = "BA.2.86 prevalence (%)",
    caption = paste0("Red dots: Wald CI = [0%, 0%] despite true prevalence > 0 | ",
                     "Wilson (1927) | Agresti & Coull (1998)")
  ) +
  theme_survinger(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    strip.text = ggplot2::element_text(size = 11, face = "bold")
  )

grDevices::png(file.path(fig_dir, "figC_ci_comparison.png"), width = 2800, height = 2000, res = 300)
print(pC); grDevices::dev.off()

file.copy(file.path(fig_dir, "figC_ci_comparison.png"),
          "vignettes/figures/figC_ci_comparison.png", overwrite = TRUE)

cat("\n=== figC: simulated rare variant, zero-width contrast visible ===\n")
