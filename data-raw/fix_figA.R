.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all(quiet = TRUE)

cols <- .surv_colors()
fig_dir <- "man/figures"

# ============================================================
# FIG A: Use HIGH-INEQUALITY scenario where difference MATTERS
# Simulate a realistic LMIC setting: capital city sequences 40%
# while rural areas sequence <1%
# ============================================================

set.seed(2024)

# 5 regions with extreme inequality (like real cross-country or LMIC within-country)
sim <- surv_simulate(
  n_regions = 5, n_weeks = 20,
  sequencing_rates = c(0.40, 0.10, 0.03, 0.008, 0.002),
  seed = 2024
)

# Make prevalence heterogeneous: high-seq region has higher prevalence
dat <- sim$sequences
regions <- sort(unique(dat$region))
target_prevs <- c(0.35, 0.25, 0.15, 0.08, 0.05)
names(target_prevs) <- regions

for (r in regions) {
  idx <- which(dat$region == r)
  p <- target_prevs[r]
  dat$lineage[idx] <- sample(c("BA.2.86", "Other"), length(idx),
                              replace = TRUE, prob = c(p, 1 - p))
}
sim$sequences <- dat

d <- surv_design(dat, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)

weighted <- surv_lineage_prevalence(d, "BA.2.86", method = "hajek")
naive <- surv_naive_prevalence(d, "BA.2.86")

w_df <- weighted$estimates[!is.na(weighted$estimates$prevalence), ]
n_df <- naive$estimates[!is.na(naive$estimates$prevalence), ]

merged <- dplyr::inner_join(
  w_df[c("time", "prevalence", "ci_lower", "ci_upper")],
  n_df[c("time", "prevalence")],
  by = "time", suffix = c("_corrected", "_naive")
)

# Real population-weighted truth
pop_shares <- sim$population$pop_total / sum(sim$population$pop_total)
true_prev <- sum(pop_shares * target_prevs[sim$population$region])
cat("True population prevalence:", round(true_prev * 100, 1), "%\n")
cat("Naive mean:", round(mean(merged$prevalence_naive) * 100, 1), "%\n")
cat("Weighted mean:", round(mean(merged$prevalence_corrected) * 100, 1), "%\n")
cat("Gini:", round(.gini_coefficient(d$strata_info$seq_rate), 3), "\n")

# Decision threshold at true prevalence level
threshold <- 0.15

merged$naive_triggers <- merged$prevalence_naive >= threshold
merged$corrected_triggers <- merged$prevalence_corrected >= threshold
n_disagree <- sum(merged$naive_triggers != merged$corrected_triggers)
cat("Weeks disagreeing on", threshold * 100, "% threshold:", n_disagree, "\n")

# Build panels
panel <- dplyr::bind_rows(
  merged |> dplyr::transmute(
    time = .data$time,
    prevalence = .data$prevalence_naive,
    panel = "Without survinger (naive)"
  ),
  merged |> dplyr::transmute(
    time = .data$time,
    prevalence = .data$prevalence_corrected,
    panel = "With survinger (design-weighted)"
  )
)
panel$panel <- factor(panel$panel,
  levels = c("Without survinger (naive)", "With survinger (design-weighted)"))

# Add truth line
truth_df <- tibble::tibble(
  time = unique(panel$time),
  prevalence = true_prev,
  panel = factor("Without survinger (naive)",
                 levels = levels(panel$panel))
)
truth_df2 <- truth_df
truth_df2$panel <- factor("With survinger (design-weighted)",
                          levels = levels(panel$panel))
truth_both <- dplyr::bind_rows(truth_df, truth_df2)

show_labels <- sort(unique(panel$time))
show_labels <- show_labels[seq(1, length(show_labels), by = 3)]

pA <- ggplot2::ggplot(panel, ggplot2::aes(x = .data$time, y = .data$prevalence, group = 1)) +
  # Truth reference
  ggplot2::geom_hline(yintercept = true_prev, linetype = "solid",
                      color = cols[["tertiary"]], linewidth = 0.6, alpha = 0.7) +
  # Action threshold
  ggplot2::geom_hline(yintercept = threshold, linetype = "dashed",
                      color = "#CC0000", linewidth = 0.5, alpha = 0.7) +
  # Estimate
  ggplot2::geom_line(ggplot2::aes(color = .data$panel), linewidth = 1) +
  ggplot2::geom_point(ggplot2::aes(color = .data$panel), size = 1.5) +
  ggplot2::scale_color_manual(values = c(cols[["secondary"]], cols[["primary"]]), guide = "none") +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                              limits = c(0, 0.45)) +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::facet_wrap(~ panel, ncol = 2) +
  ggplot2::annotate("text", x = show_labels[length(show_labels)], y = threshold + 0.015,
                    label = paste0("Action threshold (", threshold * 100, "%)"),
                    hjust = 1, size = 2.8, color = "#CC0000", fontface = "italic") +
  ggplot2::annotate("text", x = show_labels[length(show_labels)], y = true_prev - 0.015,
                    label = paste0("True prevalence (", round(true_prev * 100, 1), "%)"),
                    hjust = 1, size = 2.8, color = cols[["tertiary"]], fontface = "italic") +
  ggplot2::labs(
    title = "Naive estimates systematically overestimate prevalence in high-inequality settings",
    subtitle = paste0("Simulated: 5 regions, Gini = ",
                      round(.gini_coefficient(d$strata_info$seq_rate), 2),
                      ", sequencing rates 0.2\u201340% | ",
                      n_disagree, " weeks disagree on ", threshold * 100, "% threshold"),
    x = "Epiweek", y = "BA.2.86 prevalence",
    caption = paste0("Red dashed: action threshold | Green: true population prevalence | ",
                     "Naive overestimates by ", round((mean(merged$prevalence_naive) - true_prev) * 100, 1),
                     " pp on average")
  ) +
  theme_survinger(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
    strip.text = ggplot2::element_text(size = 11, face = "bold")
  )

grDevices::png(file.path(fig_dir, "figA_necessity.png"), width = 3200, height = 1600, res = 300)
print(pA); grDevices::dev.off()

file.copy(file.path(fig_dir, "figA_necessity.png"),
          "vignettes/figures/figA_necessity.png", overwrite = TRUE)

cat("\n=== figA FIXED: high-inequality scenario, visible difference ===\n")
