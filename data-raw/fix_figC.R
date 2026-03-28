.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all(quiet = TRUE)

cols <- .surv_colors()
fig_dir <- "man/figures"

# Use a RARE lineage where p-hat = 0 happens frequently
# This is where Wilson CI matters
coguk <- readRDS("data-raw/coguk_surveillance.rds")
d <- surv_design(
  data = coguk$sequences, strata = ~ region,
  sequencing_rate = coguk$population[c("region", "seq_rate")],
  population = coguk$population
)

# Find a rare lineage
lin_counts <- sort(table(coguk$sequences$lineage))
cat("Lineage counts:\n")
print(lin_counts)

# Use CH.1.1.1 (rare) — or the rarest non-Other
rare_lin <- names(lin_counts)[which(lin_counts > 50 & lin_counts < 3000)[1]]
if (is.na(rare_lin)) rare_lin <- "CH.1.1.1"
cat("Using rare lineage:", rare_lin, "with", lin_counts[rare_lin], "sequences\n")

w <- surv_lineage_prevalence(d, rare_lin, method = "hajek")
est <- w$estimates[!is.na(w$estimates$prevalence), ]

# Compute Wald CI for comparison
z <- stats::qnorm(0.975)
est$wald_lower <- pmax(0, est$prevalence - z * est$se)
est$wald_upper <- pmin(1, est$prevalence + z * est$se)
est$wald_width <- est$wald_upper - est$wald_lower
est$wilson_width <- est$ci_upper - est$ci_lower

# Count zero-width Wald
n_wald_zero <- sum(est$wald_width == 0)
n_wilson_zero <- sum(est$wilson_width == 0)
cat("Wald zero-width weeks:", n_wald_zero, "/", nrow(est), "\n")
cat("Wilson zero-width weeks:", n_wilson_zero, "/", nrow(est), "\n")

est$ts <- gsub("2023-", "", est$time)

# Highlight weeks where Wald fails
est$wald_fails <- est$wald_width == 0 & est$wilson_width > 0

show_labels <- sort(unique(est$ts))
show_labels <- show_labels[seq(1, length(show_labels), by = 3)]

# Two-panel: Wald (top) vs Wilson (bottom) showing CI as error bars
# This makes zero-width CIs OBVIOUS
est_long <- dplyr::bind_rows(
  est |> dplyr::transmute(
    ts = .data$ts, prevalence = .data$prevalence,
    ci_lo = .data$wald_lower, ci_hi = .data$wald_upper,
    Method = "Wald (standard)", fails = .data$wald_fails
  ),
  est |> dplyr::transmute(
    ts = .data$ts, prevalence = .data$prevalence,
    ci_lo = .data$ci_lower, ci_hi = .data$ci_upper,
    Method = "Wilson (survinger)", fails = FALSE
  )
)
est_long$Method <- factor(est_long$Method,
  levels = c("Wald (standard)", "Wilson (survinger)"))

pC <- ggplot2::ggplot(est_long, ggplot2::aes(x = .data$ts, y = .data$prevalence * 100)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$ci_lo * 100, ymax = .data$ci_hi * 100),
    width = 0.3, linewidth = 0.6, color = cols[["primary"]]
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = .data$fails), size = 2
  ) +
  ggplot2::scale_color_manual(values = c("FALSE" = cols[["primary"]], "TRUE" = "#CC0000"),
                              guide = "none") +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggplot2::facet_wrap(~ Method, ncol = 1) +
  ggplot2::labs(
    title = paste0("Wilson CI prevents zero-width intervals for rare lineage ", rare_lin),
    subtitle = paste0("COG-UK real data | Wald has ", n_wald_zero, "/", nrow(est),
                      " weeks with CI = [0,0] (red dots) | Wilson: 0 zero-width"),
    x = "Epiweek (2023)", y = paste0(rare_lin, " prevalence (%)"),
    caption = paste0("Red dots: Wald CI collapsed to zero (coverage failure) | ",
                     "Wilson (1927) JASA | Agresti & Coull (1998) | Data: COG-UK")
  ) +
  theme_survinger(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
    strip.text = ggplot2::element_text(size = 11, face = "bold")
  )

grDevices::png(file.path(fig_dir, "figC_ci_comparison.png"), width = 3000, height = 2000, res = 300)
print(pC); grDevices::dev.off()

file.copy(file.path(fig_dir, "figC_ci_comparison.png"),
          "vignettes/figures/figC_ci_comparison.png", overwrite = TRUE)

cat("\n=== figC: rare lineage, Wald zero-width visible ===\n")
