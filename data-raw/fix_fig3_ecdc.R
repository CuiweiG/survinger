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

per_rw <- d$data
per_rw$.target <- as.integer(per_rw$lineage == target_lin)
per_rw <- per_rw |>
  dplyr::group_by(.data$region, .data$epiweek) |>
  dplyr::summarise(naive_prev = mean(.data$.target), .groups = "drop")
gw <- weighted$estimates[c("time", "prevalence")]
names(gw) <- c("epiweek", "weighted_prev")
per_rw <- dplyr::left_join(per_rw, gw, by = "epiweek")
per_rw$bias <- per_rw$naive_prev - per_rw$weighted_prev
per_rw$es <- gsub("2023-", "", per_rw$epiweek)

ro <- d$strata_info[order(d$strata_info$seq_rate, decreasing = TRUE), ]
per_rw$region <- factor(per_rw$region, levels = ro$region)

show_wk <- sort(unique(per_rw$es))
show_wk <- show_wk[seq(1, length(show_wk), by = 3)]

cat("Bias range:", round(range(per_rw$bias, na.rm=TRUE) * 100, 1), "pp\n")

p3 <- ggplot2::ggplot(per_rw, ggplot2::aes(x = .data$es, y = .data$region, fill = .data$bias)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.5) +
  ggplot2::scale_fill_gradient2(
    low = cols[["primary"]], mid = "#FFFFFF", high = cols[["secondary"]],
    midpoint = 0, limits = c(-0.15, 0.15), oob = scales::squish,
    labels = function(x) paste0(ifelse(x > 0, "+", ""), round(x * 100, 1), " pp"),
    name = "Bias"
  ) +
  ggplot2::scale_x_discrete(breaks = show_wk) +
  ggplot2::labs(
    title = paste0("Per-country estimation bias for ", target_lin),
    subtitle = paste0("ECDC 5 EU countries, ordered by sequencing rate (highest top) | n = ",
                      formatC(nrow(ecdc$sequences), big.mark = ",")),
    x = "Epiweek (2023)", y = NULL,
    caption = "Blue: naive underestimates | Red: naive overestimates | Data: ECDC"
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "right", legend.key.height = ggplot2::unit(40, "pt")
  )

grDevices::png(file.path(fig_dir, "fig3_bias_heatmap.png"), width = 3200, height = 1400, res = 300)
print(p3); grDevices::dev.off()
file.copy(file.path(fig_dir, "fig3_bias_heatmap.png"),
          "vignettes/figures/fig3_bias_heatmap.png", overwrite = TRUE)

cat("=== fig3: ECDC 5-country heatmap ===\n")
