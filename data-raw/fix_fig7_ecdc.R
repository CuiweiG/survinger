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

prev_seq <- c(seq(0.0001, 0.001, 0.0001), seq(0.001, 0.01, 0.0005), seq(0.01, 0.02, 0.001))
prev_seq <- sort(unique(prev_seq))
det_probs <- vapply(prev_seq, function(p) surv_detection_probability(d, p)$overall, numeric(1))

idx_95 <- which(det_probs >= 0.95)[1]
idx_50 <- which(det_probs >= 0.50)[1]
p95 <- if (!is.na(idx_95)) prev_seq[idx_95] * 100 else NA
p50 <- if (!is.na(idx_50)) prev_seq[idx_50] * 100 else NA
cat("50% at:", p50, "% | 95% at:", p95, "%\n")

x_max <- min(max(prev_seq * 100), ifelse(!is.na(p95), p95 * 3, 1))
det_df <- tibble::tibble(prevalence = prev_seq * 100, detection = det_probs * 100)
ddc <- det_df[det_df$prevalence <= x_max, ]

p7 <- ggplot2::ggplot(ddc, ggplot2::aes(x = .data$prevalence, y = .data$detection)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = .data$detection),
                       fill = cols[["light_blue"]], alpha = 0.25) +
  ggplot2::geom_line(color = cols[["primary"]], linewidth = 1.2) +
  ggplot2::geom_hline(yintercept = 95, linetype = "dashed", color = cols[["secondary"]], linewidth = 0.5) +
  ggplot2::geom_hline(yintercept = 50, linetype = "dotted", color = cols[["neutral"]], linewidth = 0.4) +
  {if (!is.na(p95)) ggplot2::annotate("point", x = p95, y = 95, size = 4, color = cols[["secondary"]])} +
  {if (!is.na(p95)) ggplot2::annotate("label", x = p95, y = 87,
    label = paste0("95% detection\nat ", sprintf("%.2f", p95), "% prevalence"),
    size = 3, fill = "#FFFFFFDD", color = cols[["secondary"]])} +
  ggplot2::scale_x_continuous(labels = function(x) paste0(x, "%")) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = c(0, 25, 50, 80, 95, 100)) +
  ggplot2::labs(
    title = "Detection probability under current ECDC surveillance design",
    subtitle = paste0("ECDC 5 EU countries (n = ", formatC(nrow(ecdc$sequences), big.mark = ","),
                      " sequences, 46 weeks)"),
    x = "True variant prevalence", y = "Weekly detection probability",
    caption = "Based on per-stratum sequencing volumes | Data: ECDC Open Data Portal"
  ) +
  theme_survinger(base_size = 12)

grDevices::png(file.path(fig_dir, "fig7_detection.png"), width = 2800, height = 1800, res = 300)
print(p7); grDevices::dev.off()
file.copy(file.path(fig_dir, "fig7_detection.png"),
          "vignettes/figures/fig7_detection.png", overwrite = TRUE)
cat("=== fig7: ECDC data ===\n")
