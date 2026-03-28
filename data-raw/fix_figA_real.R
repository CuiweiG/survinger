.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all(quiet = TRUE)

cols <- .surv_colors()
fig_dir <- "man/figures"

# Use ECDC real data — 5 countries, 40x inequality
ecdc <- readRDS("data-raw/ecdc_surveillance.rds")
cat("ECDC:", nrow(ecdc$sequences), "sequences, 5 countries\n")

d <- surv_design(
  data = ecdc$sequences, strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population, source_type = "source_type"
)

cat("Gini:", round(.gini_coefficient(d$strata_info$seq_rate), 3), "\n")
cat("Rate ratio:", round(max(d$strata_info$seq_rate) / min(d$strata_info$seq_rate)), "x\n")

lineages <- sort(table(ecdc$sequences$lineage), decreasing = TRUE)
target_lin <- names(lineages)[!grepl("Other", names(lineages))][1]
cat("Target lineage:", target_lin, "\n")

weighted <- surv_lineage_prevalence(d, target_lin, method = "hajek")
naive <- surv_naive_prevalence(d, target_lin)

w_df <- weighted$estimates[!is.na(weighted$estimates$prevalence), ]
n_df <- naive$estimates[!is.na(naive$estimates$prevalence), ]

merged <- dplyr::inner_join(
  w_df[c("time", "prevalence", "ci_lower", "ci_upper")],
  n_df[c("time", "prevalence")],
  by = "time", suffix = c("_corrected", "_naive")
)

mean_diff <- round(mean(abs(merged$prevalence_naive - merged$prevalence_corrected)) * 100, 1)
max_diff <- round(max(abs(merged$prevalence_naive - merged$prevalence_corrected)) * 100, 1)
cat("Mean |diff|:", mean_diff, "pp\n")
cat("Max |diff|:", max_diff, "pp\n")

# Find threshold where they disagree
# Use a threshold near the crossing point
mid_prev <- median(c(merged$prevalence_naive, merged$prevalence_corrected))
threshold <- round(mid_prev * 4) / 4  # round to nearest 0.25
# Try several thresholds and pick one with most disagreement
best_th <- 0.15
best_n <- 0
for (th in seq(0.05, 0.50, by = 0.05)) {
  n_dis <- sum((merged$prevalence_naive >= th) != (merged$prevalence_corrected >= th))
  if (n_dis > best_n) { best_n <- n_dis; best_th <- th }
}
threshold <- best_th
n_disagree <- best_n
cat("Best threshold:", threshold * 100, "% with", n_disagree, "disagreements\n")

merged$ts <- gsub("2023-", "", merged$time)

# Build panels
panel <- dplyr::bind_rows(
  merged |> dplyr::transmute(
    ts = .data$ts, prevalence = .data$prevalence_naive,
    panel = "Without survinger (naive)"
  ),
  merged |> dplyr::transmute(
    ts = .data$ts, prevalence = .data$prevalence_corrected,
    panel = "With survinger (design-weighted)"
  )
)
panel$panel <- factor(panel$panel,
  levels = c("Without survinger (naive)", "With survinger (design-weighted)"))

show_labels <- sort(unique(panel$ts))
show_labels <- show_labels[seq(1, length(show_labels), by = 4)]

pA <- ggplot2::ggplot(panel, ggplot2::aes(x = .data$ts, y = .data$prevalence, group = 1)) +
  ggplot2::geom_hline(yintercept = threshold, linetype = "dashed",
                      color = "#CC0000", linewidth = 0.5, alpha = 0.7) +
  ggplot2::geom_area(fill = cols[["light_blue"]], alpha = 0.2) +
  ggplot2::geom_line(ggplot2::aes(color = .data$panel), linewidth = 1) +
  ggplot2::geom_point(ggplot2::aes(color = .data$panel), size = 1.5) +
  ggplot2::scale_color_manual(values = c(cols[["secondary"]], cols[["primary"]]), guide = "none") +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                              limits = c(0, NA)) +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::facet_wrap(~ panel, ncol = 2) +
  ggplot2::annotate("text", x = show_labels[length(show_labels)],
                    y = threshold + 0.02,
                    label = paste0("Action threshold (", threshold * 100, "%)"),
                    hjust = 1, size = 2.8, color = "#CC0000", fontface = "italic") +
  ggplot2::labs(
    title = paste0(target_lin, ": naive vs design-weighted estimates drive different decisions"),
    subtitle = paste0("ECDC real data, 5 EU countries, 40-fold sequencing inequality (Gini = ",
                      round(.gini_coefficient(d$strata_info$seq_rate), 2), ") | ",
                      n_disagree, " weeks disagree on ", threshold * 100, "% threshold | ",
                      "Mean |bias|: ", mean_diff, " pp"),
    x = "Epiweek (2023)", y = paste0(target_lin, " prevalence"),
    caption = paste0("Red dashed: hypothetical action threshold | ",
                     "Data: ECDC Open Data (n = ", formatC(nrow(ecdc$sequences), big.mark = ","), ")")
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

cat("\n=== figA: ECDC REAL DATA, 40x inequality ===\n")
