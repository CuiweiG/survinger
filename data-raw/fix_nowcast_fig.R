.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all()

cols <- .surv_colors()
fig_dir <- "man/figures"

# ============================================================
# FIG 4 FIX: Use surv_simulate with real right-truncation
# The ECDC reconstructed data has no real truncation.
# Simulated data has proper collection->report delays.
# ============================================================

set.seed(2024)
sim <- surv_simulate(n_regions = 5, n_weeks = 20,
                     sequencing_rates = c(0.15, 0.08, 0.03, 0.01, 0.005),
                     delay_params = list(mu = 12, size = 3),
                     seed = 2024)

d_full <- surv_design(sim$sequences, ~ region,
                      sim$population[c("region", "seq_rate")],
                      sim$population)

# Use a reference date that's 2 weeks before max report_date
# This creates genuine right-truncation
max_report <- max(sim$sequences$report_date)
ref_early <- max_report - 14  # pretend we're looking 2 weeks ago

# Filter to only what would have been reported by ref_early
trunc_seqs <- sim$sequences[sim$sequences$report_date <= ref_early, ]
cat("Full sequences:", nrow(sim$sequences), "\n")
cat("Truncated (reported by", as.character(ref_early), "):", nrow(trunc_seqs), "\n")
cat("Removed:", nrow(sim$sequences) - nrow(trunc_seqs), "sequences\n")

d_trunc <- surv_design(trunc_seqs, ~ region,
                       sim$population[c("region", "seq_rate")],
                       sim$population)

delay_fit <- surv_estimate_delay(d_full)  # fit on full data (known delay)

# Nowcast the truncated data
nc <- surv_nowcast_lineage(d_trunc, delay_fit, "BA.2.86",
                           ref_date = ref_early)

# Get full (eventual) counts for BA.2.86
full_counts <- sim$sequences[sim$sequences$lineage == "BA.2.86", ] |>
  dplyr::group_by(.data$epiweek) |>
  dplyr::summarise(n_full = dplyr::n(), .groups = "drop")

nc_df <- nc$estimates
nc_df <- dplyr::left_join(nc_df, full_counts, by = c("time" = "epiweek"))
nc_df$n_full[is.na(nc_df$n_full)] <- 0L

# Calculate recovery accuracy for annotation
recent <- nc_df[nc_df$is_nowcast, ]
if (nrow(recent) > 0 && any(!is.na(recent$n_full) & recent$n_full > 0)) {
  valid <- recent[!is.na(recent$n_full) & recent$n_full > 0, ]
  recovery <- mean(valid$n_estimated / valid$n_full, na.rm = TRUE)
  undercount <- mean(valid$n_observed / valid$n_full, na.rm = TRUE)
  cat("Nowcasted weeks - observed/full:", round(undercount, 2),
      "| nowcast/full:", round(recovery, 2), "\n")
}

show_labels <- sort(unique(nc_df$time))
show_labels <- show_labels[seq(1, length(show_labels), by = 3)]

p4 <- ggplot2::ggplot(nc_df, ggplot2::aes(x = .data$time)) +
  # Full eventual counts (light grey background)
  ggplot2::geom_col(ggplot2::aes(y = .data$n_full),
                    fill = "#E0E0E0", width = 0.7) +
  # Truncated observed counts (darker)
  ggplot2::geom_col(ggplot2::aes(y = .data$n_observed),
                    fill = cols[["neutral"]], alpha = 0.8, width = 0.5) +
  # Nowcast line
  ggplot2::geom_line(ggplot2::aes(y = .data$n_estimated, group = 1),
                     color = cols[["quaternary"]], linewidth = 1.1) +
  ggplot2::geom_point(ggplot2::aes(y = .data$n_estimated,
                                   shape = ifelse(.data$is_nowcast, "Nowcasted", "Complete")),
                      color = cols[["quaternary"]], size = 2.5) +
  ggplot2::scale_shape_manual(values = c("Complete" = 16, "Nowcasted" = 17)) +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::annotate("segment",
                    x = nc_df$time[which(nc_df$is_nowcast)[1]],
                    xend = nc_df$time[which(nc_df$is_nowcast)[1]],
                    y = 0, yend = max(nc_df$n_full, na.rm = TRUE) * 0.95,
                    linetype = "dashed", color = "#999999") +
  ggplot2::annotate("text",
                    x = nc_df$time[which(nc_df$is_nowcast)[1]],
                    y = max(nc_df$n_full, na.rm = TRUE) * 0.98,
                    label = "Reporting\ncutoff", hjust = 1.1,
                    size = 3, color = "#666666", fontface = "italic") +
  ggplot2::labs(
    title = "Delay-adjusted nowcasting recovers right-truncated counts",
    subtitle = paste0("Simulated data | Reference date set 2 weeks before data end | ",
                      "NegBin(\u03bc=12, k=3) delay model"),
    x = "Epiweek", y = "BA.2.86 sequence count",
    shape = NULL,
    caption = paste0("Light grey: eventual full count | Dark grey: reported by cutoff | ",
                     "Orange: nowcast | \u25b2 = nowcasted weeks")
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    legend.position = c(0.15, 0.85)
  )

grDevices::png(file.path(fig_dir, "fig4_nowcast.png"),
               width = 3000, height = 1800, res = 300)
print(p4); grDevices::dev.off()
cat("OK fig4 (simulation-based, real truncation)\n")

# ============================================================
# FIG 7 NEW: Detection probability curve
# ============================================================
ecdc <- readRDS("data-raw/ecdc_surveillance.rds")
design <- surv_design(
  data = ecdc$sequences, strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population, source_type = "source_type"
)

prev_seq <- seq(0.001, 0.10, by = 0.001)
det_probs <- vapply(prev_seq, function(p) {
  surv_detection_probability(design, p)$overall
}, numeric(1))

det_df <- tibble::tibble(prevalence = prev_seq * 100, detection = det_probs * 100)

# Key thresholds
p95 <- prev_seq[which(det_probs >= 0.95)[1]] * 100
p50 <- prev_seq[which(det_probs >= 0.50)[1]] * 100

p7 <- ggplot2::ggplot(det_df, ggplot2::aes(x = .data$prevalence, y = .data$detection)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = .data$detection),
                       fill = cols[["light_blue"]], alpha = 0.3) +
  ggplot2::geom_line(color = cols[["primary"]], linewidth = 1.2) +
  ggplot2::geom_hline(yintercept = 95, linetype = "dashed", color = cols[["secondary"]], linewidth = 0.5) +
  ggplot2::geom_hline(yintercept = 50, linetype = "dotted", color = cols[["neutral"]], linewidth = 0.5) +
  ggplot2::annotate("text", x = 8, y = 97,
                    label = "95% detection threshold", size = 3.2,
                    color = cols[["secondary"]], fontface = "italic") +
  ggplot2::annotate("point", x = p95, y = 95, size = 4, color = cols[["secondary"]]) +
  ggplot2::annotate("label", x = p95 + 0.8, y = 88,
                    label = paste0(round(p95, 2), "% prevalence\nfor 95% detection"),
                    size = 3, fill = "#FFFFFFCC", label.size = 0.3) +
  ggplot2::scale_x_continuous(labels = function(x) paste0(x, "%")) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 102)) +
  ggplot2::labs(
    title = "Variant detection probability under current surveillance design",
    subtitle = "ECDC data: 5 countries, current sequencing allocation (n = 99,093)",
    x = "True variant prevalence",
    y = "Weekly probability of detecting \u22651 sequence",
    caption = "Based on per-stratum sequencing volumes | Assumes random sampling within strata"
  ) +
  theme_survinger(base_size = 12)

grDevices::png(file.path(fig_dir, "fig7_detection.png"),
               width = 2800, height = 1800, res = 300)
print(p7); grDevices::dev.off()
cat("OK fig7 detection curve\n")

# Copy to vignettes
file.copy(file.path(fig_dir, "fig4_nowcast.png"),
          "vignettes/figures/fig4_nowcast.png", overwrite = TRUE)
file.copy(file.path(fig_dir, "fig7_detection.png"),
          "vignettes/figures/fig7_detection.png", overwrite = TRUE)

cat("\n=== ALL FIXES DONE ===\n")
