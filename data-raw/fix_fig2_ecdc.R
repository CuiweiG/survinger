.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all(quiet = TRUE)

cols <- .surv_colors()
fig_dir <- "man/figures"

ecdc <- readRDS("data-raw/ecdc_surveillance.rds")
d <- surv_design(
  data = ecdc$sequences, strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population, source_type = "source_type"
)

lineages <- sort(table(ecdc$sequences$lineage), decreasing = TRUE)
target_lin <- names(lineages)[!grepl("Other", names(lineages))][1]

weighted <- surv_lineage_prevalence(d, target_lin, method = "hajek")
naive <- surv_naive_prevalence(d, target_lin)

w_df <- weighted$estimates[!is.na(weighted$estimates$prevalence), ]
n_df <- naive$estimates[!is.na(naive$estimates$prevalence), ]

diff_mean <- round(mean(abs(w_df$prevalence - n_df$prevalence), na.rm = TRUE) * 100, 1)
max_diff <- round(max(abs(w_df$prevalence - n_df$prevalence), na.rm = TRUE) * 100, 1)

cat("Mean |diff|:", diff_mean, "pp\n")
cat("Max |diff|:", max_diff, "pp\n")

merged <- dplyr::inner_join(
  w_df[c("time", "prevalence")], n_df[c("time", "prevalence")],
  by = "time", suffix = c("_w", "_n")
)
merged$ts <- gsub("2023-", "", merged$time)

both <- dplyr::bind_rows(
  w_df |> dplyr::transmute(ts = gsub("2023-","", .data$time),
    prevalence = .data$prevalence, ci_lo = .data$ci_lower, ci_hi = .data$ci_upper,
    Method = "Design-weighted (Hajek)"),
  n_df |> dplyr::transmute(ts = gsub("2023-","", .data$time),
    prevalence = .data$prevalence, ci_lo = .data$ci_lower, ci_hi = .data$ci_upper,
    Method = "Naive (unweighted)")
)

show_labels <- sort(unique(both$ts))
show_labels <- show_labels[seq(1, length(show_labels), by = 4)]

# Find max diff week for annotation
max_idx <- which.max(abs(merged$prevalence_w - merged$prevalence_n))
max_week <- merged$ts[max_idx]

p2 <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = merged,
    ggplot2::aes(x = .data$ts, ymin = .data$prevalence_w, ymax = .data$prevalence_n, group = 1),
    fill = cols[["secondary"]], alpha = 0.20
  ) +
  ggplot2::geom_ribbon(
    data = both, ggplot2::aes(x = .data$ts, ymin = .data$ci_lo, ymax = .data$ci_hi,
                               fill = .data$Method, group = .data$Method),
    alpha = 0.08
  ) +
  ggplot2::geom_line(
    data = both, ggplot2::aes(x = .data$ts, y = .data$prevalence,
                               color = .data$Method, group = .data$Method),
    linewidth = 1
  ) +
  ggplot2::geom_point(
    data = both, ggplot2::aes(x = .data$ts, y = .data$prevalence,
                               color = .data$Method, group = .data$Method),
    size = 1.5
  ) +
  ggplot2::scale_color_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_fill_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x*100), "%"), limits = c(0, NA)) +
  ggplot2::annotate("segment",
    x = max_week, xend = max_week,
    y = merged$prevalence_w[max_idx], yend = merged$prevalence_n[max_idx],
    color = cols[["secondary"]], linewidth = 1.5,
    arrow = ggplot2::arrow(length = ggplot2::unit(3, "pt"), ends = "both")) +
  ggplot2::annotate("label", x = max_week,
    y = (merged$prevalence_w[max_idx] + merged$prevalence_n[max_idx]) / 2,
    label = paste0(max_diff, " pp"), size = 3.5, fontface = "bold",
    fill = "#FFFFFFEE", color = cols[["secondary"]]) +
  ggplot2::annotate("label", x = show_labels[2],
    y = max(n_df$prevalence, na.rm = TRUE) * 0.95,
    label = paste0("Mean |bias|: ", diff_mean, " pp\nMax |bias|: ", max_diff, " pp\n",
                   "Gini: ", round(.gini_coefficient(d$strata_info$seq_rate), 2)),
    hjust = 0, size = 3, fill = "#FFFFFFDD", color = "#333333") +
  ggplot2::labs(
    title = paste0(target_lin, ": design weighting corrects up to ", max_diff, " pp bias"),
    subtitle = paste0("ECDC real data, 5 EU countries, 40-fold inequality (n = ",
                      formatC(nrow(ecdc$sequences), big.mark = ","), ")"),
    x = "Epiweek (2023)", y = "Estimated prevalence",
    color = NULL, fill = NULL,
    caption = "Red shading: bias area | 95% Wilson CI | Data: ECDC Open Data Portal"
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
    legend.position = c(0.80, 0.90),
    legend.background = ggplot2::element_rect(fill = "#FFFFFFEE", color = NA)
  )

grDevices::png(file.path(fig_dir, "fig2_compare.png"), width = 3000, height = 1800, res = 300)
print(p2); grDevices::dev.off()
file.copy(file.path(fig_dir, "fig2_compare.png"),
          "vignettes/figures/fig2_compare.png", overwrite = TRUE)

cat("=== fig2: ECDC, mean bias", diff_mean, "pp, max", max_diff, "pp ===\n")
