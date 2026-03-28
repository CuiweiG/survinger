.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all(quiet = TRUE)

cols <- .surv_colors()
fig_dir <- "man/figures"

ecdc <- readRDS("data-raw/ecdc_surveillance.rds")

# Use Romania only — lowest sequencing country (0.3%, only 4795 seqs)
# and BQ.1 which declines to zero in late weeks
rom_seqs <- ecdc$sequences[ecdc$sequences$region == "Romania", ]
cat("Romania sequences:", nrow(rom_seqs), "\n")

rom_pop <- ecdc$population[ecdc$population$region == "Romania", ]
d_rom <- surv_design(rom_seqs, ~ region,
                     rom_pop[c("region", "seq_rate")], rom_pop)

# BQ.1 in Romania: declining lineage
target <- "BQ.1"
w <- surv_lineage_prevalence(d_rom, target, method = "hajek", min_obs = 1L)
est <- w$estimates[!is.na(w$estimates$prevalence), ]

cat("BQ.1 in Romania:\n")
cat("  Weeks:", nrow(est), "\n")
cat("  Weeks with prev=0:", sum(est$prevalence == 0), "\n")
cat("  Min n_obs:", min(est$n_obs), "\n")

# Compute Wald for comparison
z <- stats::qnorm(0.975)
est$wald_lower <- pmax(0, est$prevalence - z * est$se)
est$wald_upper <- pmin(1, est$prevalence + z * est$se)
est$wald_width <- est$wald_upper - est$wald_lower
est$wilson_width <- est$ci_upper - est$ci_lower

n_wald_zero <- sum(est$wald_width == 0)
cat("Wald zero-width:", n_wald_zero, "/", nrow(est), "\n")

est$wald_fails <- est$wald_width == 0 & est$wilson_width > 0

est_long <- dplyr::bind_rows(
  est |> dplyr::transmute(
    time = .data$time, prevalence = .data$prevalence,
    ci_lo = .data$wald_lower, ci_hi = .data$wald_upper,
    Method = "Wald (standard)", fails = .data$wald_fails
  ),
  est |> dplyr::transmute(
    time = .data$time, prevalence = .data$prevalence,
    ci_lo = .data$ci_lower, ci_hi = .data$ci_upper,
    Method = "Wilson (survinger)", fails = FALSE
  )
)
est_long$Method <- factor(est_long$Method,
  levels = c("Wald (standard)", "Wilson (survinger)"))
est_long$ts <- gsub("2023-", "", est_long$time)

show_labels <- sort(unique(est_long$ts))
show_labels <- show_labels[seq(1, length(show_labels), by = 3)]

pC <- ggplot2::ggplot(est_long, ggplot2::aes(x = .data$ts, y = .data$prevalence * 100)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$ci_lo * 100, ymax = .data$ci_hi * 100),
    width = 0.3, linewidth = 0.7, color = cols[["primary"]]
  ) +
  ggplot2::geom_point(ggplot2::aes(color = .data$fails), size = 2.5) +
  ggplot2::scale_color_manual(values = c("FALSE" = cols[["primary"]], "TRUE" = "#CC0000"),
                              guide = "none") +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggplot2::facet_wrap(~ Method, ncol = 1, scales = "free_y") +
  ggplot2::labs(
    title = paste0("Wilson CI prevents zero-width intervals for declining ", target),
    subtitle = paste0("ECDC Romania only (n = ", formatC(nrow(rom_seqs), big.mark = ","),
                      ", 0.3% sequencing rate) | Wald: ",
                      n_wald_zero, "/", nrow(est), " zero-width (red) | Wilson: 0"),
    x = "Epiweek (2023)", y = paste0(target, " prevalence (%)"),
    caption = "Red dots: Wald CI = [0%, 0%] (coverage failure) | Real ECDC data | Wilson (1927) JASA"
  ) +
  theme_survinger(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    strip.text = ggplot2::element_text(size = 11, face = "bold")
  )

grDevices::png(file.path(fig_dir, "figC_ci_comparison.png"), width = 2800, height = 2000, res = 300)
print(pC); grDevices::dev.off()
file.copy(file.path(fig_dir, "figC_ci_comparison.png"),
          "vignettes/figures/figC_ci_comparison.png", overwrite = TRUE)

cat("\n=== figC: ECDC Romania BQ.1 real data ===\n")
