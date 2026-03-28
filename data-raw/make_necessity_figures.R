.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all(quiet = TRUE)

cols <- .surv_colors()
fig_dir <- "man/figures"

coguk <- readRDS("data-raw/coguk_surveillance.rds")
design <- surv_design(
  data = coguk$sequences, strata = ~ region,
  sequencing_rate = coguk$population[c("region", "seq_rate")],
  population = coguk$population, source_type = "source_type"
)

target_lin <- "XBB.1.5"
weighted <- surv_lineage_prevalence(design, target_lin, method = "hajek")
naive <- surv_naive_prevalence(design, target_lin)

# ============================================================
# NEW FIG A: "What happens if you DON'T use survinger"
# Side-by-side faceted: Naive (wrong) vs Corrected (right)
# with annotations showing the CONSEQUENCE of each error
# ============================================================
cat("Generating necessity figures...\n")

w_df <- weighted$estimates[!is.na(weighted$estimates$prevalence), ]
n_df <- naive$estimates[!is.na(naive$estimates$prevalence), ]

# Calculate per-week decision impact
merged <- dplyr::inner_join(
  w_df[c("time", "prevalence", "ci_lower", "ci_upper")],
  n_df[c("time", "prevalence")],
  by = "time", suffix = c("_corrected", "_naive")
)
merged$bias <- merged$prevalence_naive - merged$prevalence_corrected
merged$ts <- gsub("2023-", "", merged$time)

# Scenario: decision threshold at 20% prevalence triggers public health action
threshold <- 0.20
merged$naive_triggers <- merged$prevalence_naive >= threshold
merged$corrected_triggers <- merged$prevalence_corrected >= threshold
merged$false_alarm <- merged$naive_triggers & !merged$corrected_triggers
merged$missed <- !merged$naive_triggers & merged$corrected_triggers

n_false <- sum(merged$false_alarm)
n_missed <- sum(merged$missed)
n_disagree <- sum(merged$naive_triggers != merged$corrected_triggers)

cat("Decision disagreements:", n_disagree, "weeks\n")
cat("False alarms (naive only):", n_false, "\n")
cat("Missed events (corrected only):", n_missed, "\n")

# Build panel data
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
  ggplot2::geom_area(fill = cols[["light_blue"]], alpha = 0.3) +
  ggplot2::geom_line(ggplot2::aes(color = .data$panel), linewidth = 1) +
  ggplot2::geom_point(ggplot2::aes(color = .data$panel), size = 1.5) +
  ggplot2::scale_color_manual(values = c(cols[["secondary"]], cols[["primary"]]), guide = "none") +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                              limits = c(0, NA)) +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::facet_wrap(~ panel, ncol = 2) +
  ggplot2::annotate("text", x = show_labels[length(show_labels)], y = threshold + 0.01,
                    label = "Action threshold (20%)", hjust = 1, size = 2.8,
                    color = "#CC0000", fontface = "italic") +
  ggplot2::labs(
    title = paste0("Policy impact: XBB.1.5 prevalence estimates drive different decisions"),
    subtitle = paste0("COG-UK 2023 H1 | ", n_disagree,
                      " weeks where naive and corrected estimates disagree on action threshold"),
    x = "Epiweek", y = paste0(target_lin, " prevalence"),
    caption = "Red dashed: hypothetical 20% action threshold | Data: COG-UK (n = 65,166)"
  ) +
  theme_survinger(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
    strip.text = ggplot2::element_text(size = 11, face = "bold")
  )

grDevices::png(file.path(fig_dir, "figA_necessity.png"), width = 3200, height = 1600, res = 300)
print(pA); grDevices::dev.off()
cat("OK figA\n")

# ============================================================
# NEW FIG B: Bias decomposition — WHERE does bias come from?
# Show per-region contribution to overall bias
# ============================================================

# Per-region prevalence vs population weight
region_stats <- design$data |>
  dplyr::group_by(.data$region) |>
  dplyr::summarise(
    n_seq = dplyr::n(),
    local_prev = mean(.data$lineage == target_lin),
    .groups = "drop"
  ) |>
  dplyr::left_join(design$population[c("region", "pop_total", "seq_rate")], by = "region") |>
  dplyr::mutate(
    pop_share = .data$pop_total / sum(.data$pop_total),
    seq_share = .data$n_seq / sum(.data$n_seq),
    overrepresentation = .data$seq_share / .data$pop_share,
    region_short = gsub("UK-", "", .data$region)
  )

global_naive <- mean(design$data$lineage == target_lin)
global_weighted <- mean(weighted$estimates$prevalence, na.rm = TRUE)

region_stats$contribution_to_bias <- (region_stats$seq_share - region_stats$pop_share) *
  (region_stats$local_prev - global_weighted)

cat("\nBias decomposition:\n")
print(region_stats[c("region_short", "pop_share", "seq_share", "overrepresentation",
                      "local_prev", "contribution_to_bias")])

pB <- ggplot2::ggplot(region_stats, ggplot2::aes(x = .data$overrepresentation,
                                                  y = .data$local_prev)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = cols[["neutral"]], linewidth = 0.4) +
  ggplot2::geom_hline(yintercept = global_weighted, linetype = "dotted",
                      color = cols[["primary"]], linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(size = .data$pop_share),
                      color = cols[["primary"]], alpha = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = .data$region_short),
                     vjust = -1.2, size = 3.5, fontface = "bold") +
  ggplot2::scale_size_continuous(range = c(3, 15), labels = function(x) paste0(round(x*100), "%"),
                                name = "Population share") +
  ggplot2::scale_x_continuous(labels = function(x) paste0(x, "x")) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x*100), "%")) +
  ggplot2::annotate("text", x = 0.7, y = max(region_stats$local_prev) * 0.95,
                    label = "Under-\nrepresented", size = 3, color = cols[["neutral"]], fontface = "italic") +
  ggplot2::annotate("text", x = 1.3, y = max(region_stats$local_prev) * 0.95,
                    label = "Over-\nrepresented", size = 3, color = cols[["neutral"]], fontface = "italic") +
  ggplot2::annotate("text", x = max(region_stats$overrepresentation) * 0.9,
                    y = global_weighted + 0.005,
                    label = paste0("Weighted mean: ", round(global_weighted * 100, 1), "%"),
                    size = 2.8, color = cols[["primary"]], fontface = "italic") +
  ggplot2::labs(
    title = paste0("Why bias occurs: sequencing representation vs local prevalence"),
    subtitle = paste0("Bubble size = population share | Dashed line = proportional representation (1x) | ",
                      target_lin),
    x = "Sequencing over-representation (seq_share / pop_share)",
    y = paste0("Local ", target_lin, " prevalence"),
    caption = "Points right of 1x are over-represented in the naive estimate | Data: COG-UK"
  ) +
  theme_survinger(base_size = 11)

grDevices::png(file.path(fig_dir, "figB_bias_source.png"), width = 2800, height = 1800, res = 300)
print(pB); grDevices::dev.off()
cat("OK figB\n")

# ============================================================
# NEW FIG C: CI width comparison — Wilson vs hypothetical Wald
# Shows the improvement from the methodological upgrade
# ============================================================

# Recalculate with both methods for comparison
z <- stats::qnorm(0.975)
comparison_ci <- w_df |>
  dplyr::mutate(
    wilson_width = .data$ci_upper - .data$ci_lower,
    wald_lower = pmax(0, .data$prevalence - z * .data$se),
    wald_upper = pmin(1, .data$prevalence + z * .data$se),
    wald_width = .data$wald_upper - .data$wald_lower,
    ts = gsub("2023-", "", .data$time)
  )

ci_long <- tidyr::pivot_longer(
  comparison_ci[c("ts", "wilson_width", "wald_width")],
  cols = c("wilson_width", "wald_width"),
  names_to = "method", values_to = "width"
)
ci_long$Method <- ifelse(ci_long$method == "wilson_width",
                         "Wilson (survinger)", "Wald (standard)")
ci_long$Method <- factor(ci_long$Method, levels = c("Wald (standard)", "Wilson (survinger)"))

show_ci_labels <- sort(unique(ci_long$ts))
show_ci_labels <- show_ci_labels[seq(1, length(show_ci_labels), by = 4)]

pC <- ggplot2::ggplot(ci_long, ggplot2::aes(x = .data$ts, y = .data$width * 100,
                                             fill = .data$Method, group = .data$Method)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(0.7), width = 0.6, alpha = 0.8) +
  ggplot2::scale_fill_manual(values = c(cols[["secondary"]], cols[["primary"]])) +
  ggplot2::scale_x_discrete(breaks = show_ci_labels) +
  ggplot2::labs(
    title = "Wilson intervals provide valid coverage at low prevalence",
    subtitle = "95% CI width comparison: survinger uses Wilson (never zero width)",
    x = "Epiweek", y = "95% CI width (percentage points)",
    fill = NULL,
    caption = "Wald CI = 0 when p-hat = 0 (coverage failure) | Wilson (1927) JASA"
  ) +
  theme_survinger(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "top"
  )

grDevices::png(file.path(fig_dir, "figC_ci_comparison.png"), width = 3000, height = 1600, res = 300)
print(pC); grDevices::dev.off()
cat("OK figC\n")

# Copy to vignettes
for (f in c("figA_necessity.png", "figB_bias_source.png", "figC_ci_comparison.png")) {
  file.copy(file.path(fig_dir, f), file.path("vignettes/figures", f), overwrite = TRUE)
}

cat("\n=== 3 NEW NECESSITY FIGURES DONE ===\n")
figs <- list.files(fig_dir, "\\.png$")
for (f in figs) cat(sprintf("  %s (%d KB)\n", f, round(file.info(file.path(fig_dir, f))$size / 1024)))
