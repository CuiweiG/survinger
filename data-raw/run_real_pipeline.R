.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all()

cols <- .surv_colors()
fig_dir <- "man/figures"

coguk <- readRDS("data-raw/coguk_surveillance.rds")
cat("=== FULL PIPELINE ON COG-UK REAL DATA ===\n")
cat("Sequences:", coguk$n_sequences, "\n")
cat("Source:", coguk$source, "\n\n")

# 1. Design
design <- surv_design(
  data = coguk$sequences, strata = ~ region,
  sequencing_rate = coguk$population[c("region", "seq_rate")],
  population = coguk$population,
  source_type = "source_type"
)
cat("--- Design ---\n")
print(design)

# 2. Report
cat("\n--- System Report ---\n")
report <- surv_report(design)

# 3. Target lineage
target_lin <- "XBB.1.5"
cat("\n--- Target:", target_lin, "---\n")
weighted <- surv_lineage_prevalence(design, target_lin, method = "hajek")
naive <- surv_naive_prevalence(design, target_lin)

w_mean <- round(mean(weighted$estimates$prevalence, na.rm = TRUE), 4)
n_mean <- round(mean(naive$estimates$prevalence, na.rm = TRUE), 4)
diff_mean <- round(mean(abs(weighted$estimates$prevalence -
                            naive$estimates$prevalence), na.rm = TRUE), 4)
cat("Weighted mean:", w_mean, "\n")
cat("Naive mean:", n_mean, "\n")
cat("Mean |diff|:", diff_mean, "\n")

# 4. Delay + Nowcast
delay <- surv_estimate_delay(design)
cat("\n--- Delay fit ---\n")
print(delay)

# 5. Adjusted
adj <- surv_adjusted_prevalence(design, delay, target_lin)

# 6. Allocation
alloc <- surv_optimize_allocation(design, "min_mse", total_capacity = 2000)

# 7. Detection
det <- surv_detection_probability(design, 0.01)
cat("\nDetection at 1%:", round(det$overall, 3), "\n")

# ============================================================
# REGENERATE ALL FIGURES WITH COG-UK REAL DATA
# ============================================================

cat("\n=== Generating figures ===\n")

# FIG 1: Inequality
info <- design$strata_info
info$label <- gsub("UK-", "", info$region)
info <- info[order(info$seq_rate), ]
info$label <- factor(info$label, levels = info$label)
ratio <- round(max(info$seq_rate) / min(info$seq_rate), 1)

p1 <- ggplot2::ggplot(info, ggplot2::aes(x = .data$label, y = .data$seq_rate)) +
  ggplot2::geom_segment(ggplot2::aes(xend = .data$label, y = 0, yend = .data$seq_rate),
                        color = cols[["primary"]], linewidth = 1.5) +
  ggplot2::geom_point(size = 5, color = cols[["primary"]]) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(round(.data$seq_rate * 100, 1), "%")),
                     hjust = -0.3, size = 4, color = "#333333") +
  ggplot2::geom_hline(yintercept = mean(info$seq_rate), linetype = "dashed",
                      color = cols[["secondary"]], linewidth = 0.5) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x*100,1), "%"),
                              expand = ggplot2::expansion(mult = c(0, 0.25))) +
  ggplot2::coord_flip() +
  ggplot2::labs(title = paste0("Sequencing rates vary ", ratio, "-fold across UK nations"),
                subtitle = paste0("COG-UK individual-level data, 2023 H1 (n = ",
                                  formatC(nrow(coguk$sequences), big.mark = ","), " sequences)"),
                x = NULL, y = "Sequencing rate (% of confirmed cases)",
                caption = "Data: COG-UK Consortium (CLIMB) | Dashed line: mean rate") +
  theme_survinger(base_size = 12)
grDevices::png(file.path(fig_dir, "fig1_inequality.png"), width = 2400, height = 1400, res = 300)
print(p1); grDevices::dev.off()
cat("OK fig1\n")

# FIG 2: Weighted vs Naive
w_df <- weighted$estimates[!is.na(weighted$estimates$prevalence), ]
w_df$Method <- "Design-weighted (Hajek)"
n_df <- naive$estimates[!is.na(naive$estimates$prevalence), ]
n_df$Method <- "Naive (unweighted)"
combined <- dplyr::bind_rows(w_df, n_df)
combined$time_short <- gsub("2023-", "", combined$time)

max_diff_val <- round(max(abs(w_df$prevalence - n_df$prevalence), na.rm = TRUE) * 100, 1)

show_labels <- sort(unique(combined$time_short))
show_labels <- show_labels[seq(1, length(show_labels), by = 3)]

p2 <- ggplot2::ggplot(combined, ggplot2::aes(
  x = .data$time_short, y = .data$prevalence,
  color = .data$Method, fill = .data$Method, group = .data$Method
)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
                       alpha = 0.10, color = NA) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::scale_color_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_fill_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_x_discrete(breaks = show_labels) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x*100), "%"), limits = c(0, NA)) +
  ggplot2::annotate("label", x = show_labels[2], y = max(w_df$prevalence, na.rm=TRUE) * 0.95,
                    label = paste0("Mean |bias|: ", round(diff_mean*100, 1), " pp\nMax |bias|: ", max_diff_val, " pp"),
                    hjust = 0, size = 3.2, fill = "#FFFFFFCC", color = "#333333") +
  ggplot2::labs(title = paste0(target_lin, ": design weighting corrects surveillance bias"),
                subtitle = paste0("COG-UK real data, 4 UK nations, 2023 H1 (n = ",
                                  formatC(nrow(coguk$sequences), big.mark=","), ")"),
                x = "Epiweek (2023)", y = "Estimated prevalence",
                color = NULL, fill = NULL,
                caption = "Shaded: 95% Wald CI | Data: COG-UK Consortium") +
  theme_survinger(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
                 legend.position = c(0.70, 0.90),
                 legend.background = ggplot2::element_rect(fill = "#FFFFFFDD", color = NA))
grDevices::png(file.path(fig_dir, "fig2_compare.png"), width = 3000, height = 1800, res = 300)
print(p2); grDevices::dev.off()
cat("OK fig2\n")

# FIG 3: Bias heatmap
per_rw <- design$data
per_rw$.target <- as.integer(per_rw$lineage == target_lin)
per_rw <- per_rw |>
  dplyr::group_by(.data$region, .data$epiweek) |>
  dplyr::summarise(naive_prev = mean(.data$.target), .groups = "drop")
global_w <- weighted$estimates[c("time", "prevalence")]
names(global_w) <- c("epiweek", "weighted_prev")
per_rw <- dplyr::left_join(per_rw, global_w, by = "epiweek")
per_rw$bias <- per_rw$naive_prev - per_rw$weighted_prev
per_rw$ew_short <- gsub("2023-", "", per_rw$epiweek)
per_rw$region_short <- gsub("UK-", "", per_rw$region)

rate_order <- design$strata_info[order(design$strata_info$seq_rate, decreasing = TRUE), ]
per_rw$region_short <- factor(gsub("UK-", "", per_rw$region),
                              levels = gsub("UK-", "", rate_order$region))

show_wk <- sort(unique(per_rw$ew_short))
show_wk <- show_wk[seq(1, length(show_wk), by = 3)]

p3 <- ggplot2::ggplot(per_rw, ggplot2::aes(
  x = .data$ew_short, y = .data$region_short, fill = .data$bias
)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.5) +
  ggplot2::scale_fill_gradient2(low = cols[["primary"]], mid = "#FFFFFF", high = cols[["secondary"]],
    midpoint = 0, limits = c(-0.15, 0.15), oob = scales::squish,
    labels = function(x) paste0(ifelse(x>0, "+", ""), round(x*100, 1), " pp"),
    name = "Bias") +
  ggplot2::scale_x_discrete(breaks = show_wk) +
  ggplot2::labs(title = paste0("Per-nation estimation bias for ", target_lin),
                subtitle = "Nations ordered by sequencing rate (highest top) | Positive = naive overestimates",
                x = "Epiweek (2023)", y = NULL,
                caption = "Blue: naive < weighted | Red: naive > weighted | Data: COG-UK") +
  theme_survinger(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                 legend.position = "right", legend.key.height = ggplot2::unit(40, "pt"))
grDevices::png(file.path(fig_dir, "fig3_bias_heatmap.png"), width = 3200, height = 1200, res = 300)
print(p3); grDevices::dev.off()
cat("OK fig3\n")

# FIG 5: Allocation
H <- design$n_strata
strats <- c("equal", "proportional", "min_mse")
alloc_details <- purrr::map_dfr(strats, function(s) {
  if (s == "equal") {
    n <- rep(2000L %/% H, H)
    n[1] <- n[1] + 2000L - sum(n)
    tibble::tibble(region = gsub("UK-", "", design$strata_info$region), n = n, strategy = s)
  } else if (s == "proportional") {
    pi_h <- design$population$pop_total / sum(design$population$pop_total)
    n <- as.integer(round(pi_h * 2000))
    n[1] <- n[1] + 2000L - sum(n)
    tibble::tibble(region = gsub("UK-", "", design$population$region), n = n, strategy = s)
  } else {
    res <- surv_optimize_allocation(design, s, total_capacity = 2000)
    tibble::tibble(region = gsub("UK-", "", res$allocation$region),
                   n = res$allocation$n_allocated, strategy = s)
  }
})
alloc_details$strategy <- factor(alloc_details$strategy,
  levels = c("equal", "proportional", "min_mse"),
  labels = c("Equal", "Population-proportional", "MSE-optimal (Neyman)"))

p5 <- ggplot2::ggplot(alloc_details, ggplot2::aes(x = .data$region, y = .data$n, fill = .data$strategy)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::scale_fill_manual(values = c(cols[["neutral"]], cols[["light_blue"]], cols[["primary"]])) +
  ggplot2::labs(title = "MSE-optimal allocation concentrates resources in England",
                subtitle = "Distribution of 2,000 weekly sequences across 4 UK nations",
                x = NULL, y = "Sequences allocated per week", fill = "Strategy",
                caption = "MSE-optimal: Neyman allocation | Data: COG-UK") +
  theme_survinger(base_size = 12)
grDevices::png(file.path(fig_dir, "fig5_allocation.png"), width = 2800, height = 1600, res = 300)
print(p5); grDevices::dev.off()
cat("OK fig5\n")

# Copy all to vignettes/figures
for (f in list.files(fig_dir, pattern = "\\.png$", full.names = TRUE)) {
  file.copy(f, file.path("vignettes/figures", basename(f)), overwrite = TRUE)
}

cat("\n=== ALL COG-UK FIGURES DONE ===\n")
cat("Key result: Mean |weighted - naive| =", round(diff_mean * 100, 1), "pp\n")
cat("Gini:", round(report$gini, 3), "\n")
cat("Rate ratio:", ratio, "x\n")
