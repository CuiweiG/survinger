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

merged <- dplyr::inner_join(
  w_df[c("time", "prevalence")],
  n_df[c("time", "prevalence")],
  by = "time", suffix = c("_weighted", "_naive")
)
merged$ts <- gsub("2023-", "", merged$time)
merged$diff <- merged$prevalence_naive - merged$prevalence_weighted
merged$abs_diff <- abs(merged$diff)

# Both lines in same plot + shaded difference
both <- dplyr::bind_rows(
  merged |> dplyr::transmute(ts = .data$ts, prevalence = .data$prevalence_naive,
                              Method = "Naive (unweighted)"),
  merged |> dplyr::transmute(ts = .data$ts, prevalence = .data$prevalence_weighted,
                              Method = "Design-weighted (Hajek)")
)

show_labels <- sort(unique(merged$ts))
show_labels <- show_labels[seq(1, length(show_labels), by = 4)]

mean_diff <- round(mean(merged$abs_diff) * 100, 1)
max_diff <- round(max(merged$abs_diff) * 100, 1)
max_week <- merged$ts[which.max(merged$abs_diff)]

pA <- ggplot2::ggplot() +
  # Shaded area between the two lines (the bias)
  ggplot2::geom_ribbon(
    data = merged,
    ggplot2::aes(x = .data$ts, ymin = .data$prevalence_weighted,
                 ymax = .data$prevalence_naive, group = 1),
    fill = cols[["secondary"]], alpha = 0.20
  ) +
  # Both lines
  ggplot2::geom_line(
    data = both,
    ggplot2::aes(x = .data$ts, y = .data$prevalence,
                 color = .data$Method, group = .data$Method),
    linewidth = 1
  ) +
  ggplot2::geom_point(
    data = both,
    ggplot2::aes(x = .data$ts, y = .data$prevalence,
                 color = .data$Method, group = .data$Method),
    size = 1.5
  ) +
  ggplot2::scale_color_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                              limits = c(0, NA)) +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  # Annotate max difference
  ggplot2::annotate("segment",
    x = max_week, xend = max_week,
    y = merged$prevalence_weighted[which.max(merged$abs_diff)],
    yend = merged$prevalence_naive[which.max(merged$abs_diff)],
    color = cols[["secondary"]], linewidth = 1.5, arrow = ggplot2::arrow(length = ggplot2::unit(3, "pt"), ends = "both")
  ) +
  ggplot2::annotate("label",
    x = max_week, y = (merged$prevalence_naive[which.max(merged$abs_diff)] +
                        merged$prevalence_weighted[which.max(merged$abs_diff)]) / 2,
    label = paste0(max_diff, " pp"), size = 3.5, fontface = "bold",
    fill = "#FFFFFFEE", color = cols[["secondary"]]
  ) +
  # Summary box
  ggplot2::annotate("label",
    x = show_labels[2], y = max(merged$prevalence_naive) * 0.95,
    label = paste0("Mean |bias|: ", mean_diff, " pp\nMax |bias|: ", max_diff, " pp\n",
                   "Shaded area = estimation error\nfrom ignoring sequencing inequality"),
    hjust = 0, size = 3, fill = "#FFFFFFDD", color = "#333333"
  ) +
  ggplot2::labs(
    title = paste0("The cost of ignoring sequencing inequality: up to ",
                   max_diff, " pp error in ", target_lin, " prevalence"),
    subtitle = paste0("ECDC real data | 5 EU countries | 40-fold sequencing inequality (Gini = ",
                      round(.gini_coefficient(d$strata_info$seq_rate), 2),
                      ") | n = ", formatC(nrow(ecdc$sequences), big.mark = ",")),
    x = "Epiweek (2023)", y = paste0(target_lin, " prevalence"),
    color = NULL,
    caption = "Red shading: bias from naive estimation | Data: ECDC Open Data Portal"
  ) +
  theme_survinger(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    legend.position = c(0.80, 0.92),
    legend.background = ggplot2::element_rect(fill = "#FFFFFFEE", color = NA)
  )

grDevices::png(file.path(fig_dir, "figA_necessity.png"), width = 3000, height = 1800, res = 300)
print(pA); grDevices::dev.off()

file.copy(file.path(fig_dir, "figA_necessity.png"),
          "vignettes/figures/figA_necessity.png", overwrite = TRUE)

cat("=== figA: OVERLAY with shaded bias, max annotation ===\n")
cat("Mean bias:", mean_diff, "pp | Max:", max_diff, "pp\n")
