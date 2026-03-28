.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all()

cols <- .surv_colors()
fig_dir <- "man/figures"

set.seed(2024)
sim <- surv_simulate(n_regions = 5, n_weeks = 20,
                     sequencing_rates = c(0.15, 0.08, 0.03, 0.01, 0.005),
                     delay_params = list(mu = 12, size = 3), seed = 2024)

# Reference date = max collection date (realistic: "today")
# This means sequences collected today but not yet reported are missing
ref_date <- max(sim$sequences$collection_date)
trunc_seqs <- sim$sequences[sim$sequences$report_date <= ref_date, ]
cat("Full:", nrow(sim$sequences), "Truncated:", nrow(trunc_seqs),
    "Removed:", nrow(sim$sequences) - nrow(trunc_seqs), "\n")

d_full <- surv_design(sim$sequences, ~ region,
                      sim$population[c("region", "seq_rate")], sim$population)
d_trunc <- surv_design(trunc_seqs, ~ region,
                       sim$population[c("region", "seq_rate")], sim$population)

delay_fit <- surv_estimate_delay(d_full)
nc <- surv_nowcast_lineage(d_trunc, delay_fit, "BA.2.86", ref_date = ref_date)

# Full counts for comparison
full_counts <- sim$sequences[sim$sequences$lineage == "BA.2.86", ] |>
  dplyr::group_by(.data$epiweek) |>
  dplyr::summarise(n_full = dplyr::n(), .groups = "drop")

nc_df <- dplyr::left_join(nc$estimates, full_counts, by = c("time" = "epiweek"))
nc_df$n_full[is.na(nc_df$n_full)] <- 0L

show_labels <- sort(unique(nc_df$time))
show_labels <- show_labels[seq(1, length(show_labels), by = 3)]

# Find the first nowcasted week for annotation
nc_start <- nc_df$time[which(nc_df$is_nowcast)[1]]

p4 <- ggplot2::ggplot(nc_df, ggplot2::aes(x = .data$time)) +
  ggplot2::geom_col(ggplot2::aes(y = .data$n_full), fill = "#E0E0E0", width = 0.7) +
  ggplot2::geom_col(ggplot2::aes(y = .data$n_observed), fill = cols[["neutral"]], alpha = 0.8, width = 0.5) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(0, .data$n_estimated - 1.96 * .data$se),
                                    ymax = .data$n_estimated + 1.96 * .data$se, group = 1),
                       fill = cols[["quaternary"]], alpha = 0.12) +
  ggplot2::geom_line(ggplot2::aes(y = .data$n_estimated, group = 1),
                     color = cols[["quaternary"]], linewidth = 1.1) +
  ggplot2::geom_point(ggplot2::aes(y = .data$n_estimated,
                                   shape = ifelse(.data$is_nowcast, "Nowcasted", "Complete")),
                      color = cols[["quaternary"]], size = 2.5) +
  ggplot2::scale_shape_manual(values = c("Complete" = 16, "Nowcasted" = 17)) +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  {if (!is.na(nc_start)) ggplot2::geom_vline(xintercept = nc_start,
                                              linetype = "dashed", color = "#999999")} +
  {if (!is.na(nc_start)) ggplot2::annotate("text", x = nc_start, y = max(nc_df$n_full, na.rm = TRUE) * 0.95,
                                           label = "Nowcast\nwindow \u2192", hjust = 1.1,
                                           size = 3, color = "#666666", fontface = "italic")} +
  ggplot2::labs(
    title = "Delay-adjusted nowcasting recovers incomplete recent counts",
    subtitle = paste0("Simulated surveillance | ref_date = max(collection_date) | ",
                      "NegBin(\u03bc=", round(delay_fit$parameters$mu[1], 1),
                      ", k=", round(delay_fit$parameters$size[1], 1), ") delay"),
    x = "Epiweek", y = "BA.2.86 sequence count",
    shape = NULL,
    caption = paste0("Light grey: eventual full count | Dark grey: reported by ref_date | ",
                     "Orange + CI: delay-corrected nowcast | \u25b2 = incomplete weeks")
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    legend.position = c(0.15, 0.85)
  )

grDevices::png(file.path(fig_dir, "fig4_nowcast.png"), width = 3000, height = 1800, res = 300)
print(p4); grDevices::dev.off()

file.copy(file.path(fig_dir, "fig4_nowcast.png"),
          "vignettes/figures/fig4_nowcast.png", overwrite = TRUE)

cat("=== fig4 nowcast FIXED (visible truncation) ===\n")
