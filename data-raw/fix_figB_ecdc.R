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
global_weighted <- mean(weighted$estimates$prevalence, na.rm = TRUE)

region_stats <- d$data |>
  dplyr::group_by(.data$region) |>
  dplyr::summarise(
    n_seq = dplyr::n(),
    local_prev = mean(.data$lineage == target_lin),
    .groups = "drop"
  ) |>
  dplyr::left_join(d$population[c("region", "pop_total", "seq_rate")], by = "region") |>
  dplyr::mutate(
    pop_share = .data$pop_total / sum(.data$pop_total),
    seq_share = .data$n_seq / sum(.data$n_seq),
    overrepresentation = .data$seq_share / .data$pop_share
  )

cat("Region stats:\n")
print(region_stats[c("region", "seq_rate", "local_prev", "overrepresentation")])
cat("Prevalence range:", round(range(region_stats$local_prev) * 100, 1), "%\n")

pB <- ggplot2::ggplot(region_stats, ggplot2::aes(
  x = .data$overrepresentation, y = .data$local_prev)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = cols[["neutral"]], linewidth = 0.4) +
  ggplot2::geom_hline(yintercept = global_weighted, linetype = "dotted",
                      color = cols[["primary"]], linewidth = 0.4) +
  ggplot2::geom_point(ggplot2::aes(size = .data$pop_share),
                      color = cols[["primary"]], alpha = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = .data$region),
                     vjust = -1.5, size = 3.5, fontface = "bold") +
  ggplot2::scale_size_continuous(range = c(4, 18),
                                labels = function(x) paste0(round(x * 100), "%"),
                                name = "Population share") +
  ggplot2::scale_x_continuous(labels = function(x) paste0(x, "x")) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  ggplot2::annotate("text", x = 0.3, y = max(region_stats$local_prev) * 0.95,
                    label = "Under-\nrepresented", size = 3, color = cols[["neutral"]], fontface = "italic") +
  ggplot2::annotate("text", x = 2.5, y = max(region_stats$local_prev) * 0.95,
                    label = "Over-\nrepresented", size = 3, color = cols[["neutral"]], fontface = "italic") +
  ggplot2::annotate("text", x = max(region_stats$overrepresentation) * 0.85,
                    y = global_weighted + 0.01,
                    label = paste0("Weighted mean: ", round(global_weighted * 100, 1), "%"),
                    size = 2.8, color = cols[["primary"]], fontface = "italic") +
  ggplot2::labs(
    title = paste0("Why bias occurs: sequencing representation vs local ", target_lin, " prevalence"),
    subtitle = paste0("ECDC 5 EU countries | Bubble = population share | ",
                      "Dashed = proportional (1x) | n = ",
                      formatC(nrow(ecdc$sequences), big.mark = ",")),
    x = "Sequencing over-representation (seq_share / pop_share)",
    y = paste0("Local ", target_lin, " prevalence"),
    caption = "Points right of 1x are over-represented in the naive estimate | Data: ECDC"
  ) +
  theme_survinger(base_size = 11)

grDevices::png(file.path(fig_dir, "figB_bias_source.png"), width = 2800, height = 1800, res = 300)
print(pB); grDevices::dev.off()

file.copy(file.path(fig_dir, "figB_bias_source.png"),
          "vignettes/figures/figB_bias_source.png", overwrite = TRUE)

cat("=== figB: ECDC data, wider prevalence spread ===\n")
