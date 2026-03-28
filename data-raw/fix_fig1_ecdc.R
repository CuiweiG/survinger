.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all(quiet = TRUE)

cols <- .surv_colors()
fig_dir <- "man/figures"

ecdc <- readRDS("data-raw/ecdc_surveillance.rds")
d <- surv_design(
  data = ecdc$sequences, strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population
)

info <- d$strata_info
info$label <- info$region
info <- info[order(info$seq_rate), ]
info$label <- factor(info$label, levels = info$label)
ratio <- round(max(info$seq_rate) / min(info$seq_rate))

p1 <- ggplot2::ggplot(info, ggplot2::aes(x = .data$label, y = .data$seq_rate)) +
  ggplot2::geom_segment(ggplot2::aes(xend = .data$label, y = 0, yend = .data$seq_rate),
                        color = cols[["primary"]], linewidth = 1.5) +
  ggplot2::geom_point(size = 5, color = cols[["primary"]]) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(round(.data$seq_rate * 100, 1), "%")),
                     hjust = -0.3, size = 4, color = "#333333") +
  ggplot2::geom_hline(yintercept = mean(info$seq_rate), linetype = "dashed",
                      color = cols[["secondary"]], linewidth = 0.5) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100, 1), "%"),
                              expand = ggplot2::expansion(mult = c(0, 0.25))) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = paste0("Sequencing rates vary ", ratio, "-fold across European countries"),
    subtitle = paste0("ECDC variant surveillance, 2023 (n = ",
                      formatC(nrow(ecdc$sequences), big.mark = ","), " sequences, 5 EU countries)"),
    x = NULL, y = "Sequencing rate (% of confirmed cases sequenced)",
    caption = "Data: ECDC Open Data Portal | Dashed line: mean rate"
  ) +
  theme_survinger(base_size = 12)

grDevices::png(file.path(fig_dir, "fig1_inequality.png"), width = 2400, height = 1500, res = 300)
print(p1); grDevices::dev.off()
file.copy(file.path(fig_dir, "fig1_inequality.png"),
          "vignettes/figures/fig1_inequality.png", overwrite = TRUE)

cat("=== fig1: ECDC 5-country, 40x inequality ===\n")
