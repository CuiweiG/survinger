.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all()

cols <- .surv_colors()
fig_dir <- "man/figures"

# USE COG-UK DATA, NOT ECDC
coguk <- readRDS("data-raw/coguk_surveillance.rds")
design <- surv_design(
  data = coguk$sequences, strata = ~ region,
  sequencing_rate = coguk$population[c("region", "seq_rate")],
  population = coguk$population, source_type = "source_type"
)

n_total <- nrow(coguk$sequences)

prev_seq <- c(seq(0.0001, 0.001, by = 0.0001),
              seq(0.001, 0.01, by = 0.0005),
              seq(0.01, 0.02, by = 0.001))
prev_seq <- sort(unique(prev_seq))

det_probs <- vapply(prev_seq, function(p) {
  surv_detection_probability(design, p)$overall
}, numeric(1))

det_df <- tibble::tibble(prevalence = prev_seq * 100, detection = det_probs * 100)

idx_95 <- which(det_probs >= 0.95)[1]
idx_80 <- which(det_probs >= 0.80)[1]
idx_50 <- which(det_probs >= 0.50)[1]

p95 <- if (!is.na(idx_95)) prev_seq[idx_95] * 100 else NA
p80 <- if (!is.na(idx_80)) prev_seq[idx_80] * 100 else NA
p50 <- if (!is.na(idx_50)) prev_seq[idx_50] * 100 else NA

cat("50% detection at:", p50, "%\n")
cat("80% detection at:", p80, "%\n")
cat("95% detection at:", p95, "%\n")

x_max <- min(max(prev_seq * 100), ifelse(!is.na(p95), p95 * 3, 1))
det_df_clip <- det_df[det_df$prevalence <= x_max, ]

p7 <- ggplot2::ggplot(det_df_clip, ggplot2::aes(x = .data$prevalence, y = .data$detection)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = .data$detection),
                       fill = cols[["light_blue"]], alpha = 0.25) +
  ggplot2::geom_line(color = cols[["primary"]], linewidth = 1.2) +
  ggplot2::geom_hline(yintercept = 95, linetype = "dashed",
                      color = cols[["secondary"]], linewidth = 0.5) +
  ggplot2::geom_hline(yintercept = 80, linetype = "dotted",
                      color = cols[["neutral"]], linewidth = 0.4) +
  ggplot2::geom_hline(yintercept = 50, linetype = "dotted",
                      color = cols[["neutral"]], linewidth = 0.4) +
  {if (!is.na(p95)) ggplot2::annotate("point", x = p95, y = 95,
                                      size = 4, color = cols[["secondary"]])} +
  {if (!is.na(p80)) ggplot2::annotate("point", x = p80, y = 80,
                                      size = 3, color = cols[["neutral"]])} +
  {if (!is.na(p50)) ggplot2::annotate("point", x = p50, y = 50,
                                      size = 3, color = cols[["neutral"]])} +
  {if (!is.na(p95)) ggplot2::annotate("label", x = p95, y = 87,
    label = paste0("95% detection\nat ", sprintf("%.2f", p95), "% prevalence"),
    size = 3, fill = "#FFFFFFDD", label.padding = ggplot2::unit(4, "pt"),
    color = cols[["secondary"]])} +
  {if (!is.na(p50)) ggplot2::annotate("label", x = p50, y = 42,
    label = paste0("50% at ", sprintf("%.3f", p50), "%"),
    size = 2.8, fill = "#FFFFFFDD", label.padding = ggplot2::unit(3, "pt"),
    color = cols[["neutral"]])} +
  ggplot2::scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    breaks = scales::pretty_breaks(n = 8)
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    breaks = c(0, 25, 50, 80, 95, 100)
  ) +
  ggplot2::labs(
    title = "Variant detection probability under current surveillance design",
    subtitle = paste0("COG-UK 4-nation design (n = ",
                      formatC(n_total, big.mark = ","),
                      " sequences, 26 weeks) | P(detect \u2265 1 per week)"),
    x = "True variant prevalence in population",
    y = "Weekly detection probability",
    caption = "Based on per-stratum sequencing volumes | Data: COG-UK Consortium"
  ) +
  theme_survinger(base_size = 12)

grDevices::png(file.path(fig_dir, "fig7_detection.png"),
               width = 2800, height = 1800, res = 300)
print(p7); grDevices::dev.off()

file.copy(file.path(fig_dir, "fig7_detection.png"),
          "vignettes/figures/fig7_detection.png", overwrite = TRUE)

cat("=== fig7 detection FIXED — now uses COG-UK data ===\n")
