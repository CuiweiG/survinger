.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all(quiet = TRUE)

cols <- .surv_colors()
fig_dir <- "man/figures"

# ECDC data has simulated report_date (NegBin delays added during processing)
# This is transparently documented. Use it for nowcast demo.
ecdc <- readRDS("data-raw/ecdc_surveillance.rds")

d <- surv_design(
  data = ecdc$sequences, strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population, source_type = "source_type"
)

delay <- surv_estimate_delay(d)
cat("Delay mu:", round(delay$parameters$mu[1], 1), "\n")

# Truncate: use ref_date = max(collection_date) to create real right-truncation
ref_date <- max(ecdc$sequences$collection_date)
trunc_seqs <- ecdc$sequences[ecdc$sequences$report_date <= ref_date, ]
cat("Full:", nrow(ecdc$sequences), "Truncated:", nrow(trunc_seqs),
    "Removed:", nrow(ecdc$sequences) - nrow(trunc_seqs), "\n")

d_trunc <- surv_design(trunc_seqs, ~ region,
                       ecdc$population[c("region", "seq_rate")], ecdc$population,
                       source_type = "source_type")

target_lin <- "XBB.1.5-like"
nc <- surv_nowcast_lineage(d_trunc, delay, target_lin, ref_date = ref_date)

# Full counts for comparison
full_counts <- ecdc$sequences[ecdc$sequences$lineage == target_lin, ] |>
  dplyr::group_by(.data$epiweek) |>
  dplyr::summarise(n_full = dplyr::n(), .groups = "drop")

nc_df <- dplyr::left_join(nc$estimates, full_counts, by = c("time" = "epiweek"))
nc_df$n_full[is.na(nc_df$n_full)] <- 0L
nc_df$ts <- gsub("2023-", "", nc_df$time)

show_labels <- sort(unique(nc_df$ts))
show_labels <- show_labels[seq(1, length(show_labels), by = 4)]
nc_start <- nc_df$ts[which(nc_df$is_nowcast)[1]]

p4 <- ggplot2::ggplot(nc_df, ggplot2::aes(x = .data$ts)) +
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
                                              linetype = "dashed", color = "#999")} +
  ggplot2::labs(
    title = paste0("Delay-adjusted nowcasting for ", target_lin),
    subtitle = paste0("ECDC 5 EU countries | ref_date = max(collection_date) | ",
                      "NegBin(\u03bc=", round(delay$parameters$mu[1], 1),
                      ") delay | n = ", formatC(nrow(ecdc$sequences), big.mark = ",")),
    x = "Epiweek (2023)", y = paste0(target_lin, " count"),
    shape = NULL,
    caption = paste0("Light grey: eventual full count | Dark grey: reported by ref_date | ",
                     "Orange: nowcast | Note: reporting delays modeled from ECDC data")
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

cat("=== fig4: ECDC data with modeled delays ===\n")
