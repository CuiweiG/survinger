.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all()

ecdc <- readRDS("data-raw/ecdc_surveillance.rds")
cols <- .surv_colors()

design <- surv_design(
  data = ecdc$sequences, strata = ~ region,
  sequencing_rate = ecdc$population[c("region", "seq_rate")],
  population = ecdc$population, source_type = "source_type"
)

lineages <- sort(table(ecdc$sequences$lineage), decreasing = TRUE)
target_lin <- names(lineages)[!grepl("Other", names(lineages))][1]

weighted <- surv_lineage_prevalence(design, target_lin, method = "hajek")
naive <- surv_naive_prevalence(design, target_lin)
delay <- surv_estimate_delay(design)
nc <- surv_nowcast_lineage(design, delay, target_lin)
adj <- surv_adjusted_prevalence(design, delay, target_lin)

fig_dir <- "vignettes/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# FIG 1: Inequality lollipop
info <- design$strata_info
info$label <- info$region
info <- info[order(info$seq_rate), ]
info$label <- factor(info$label, levels = info$label)
avg_rate <- mean(info$seq_rate)

p1 <- ggplot2::ggplot(info, ggplot2::aes(x = .data$label, y = .data$seq_rate)) +
  ggplot2::geom_segment(ggplot2::aes(xend = .data$label, y = 0, yend = .data$seq_rate),
                        color = cols[["primary"]], linewidth = 1.2) +
  ggplot2::geom_point(size = 4, color = cols[["primary"]]) +
  ggplot2::geom_hline(yintercept = avg_rate, linetype = "dashed",
                      color = cols[["secondary"]], linewidth = 0.6) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x*100,1), "%"),
                              expand = ggplot2::expansion(mult = c(0, 0.1))) +
  ggplot2::coord_flip() +
  ggplot2::labs(title = "Sequencing rates vary 40-fold across countries",
                subtitle = "ECDC COVID-19 variant surveillance, 5 European countries",
                x = NULL, y = "Sequencing rate (% of positive cases)",
                caption = "Data: ECDC Open Data Portal") +
  theme_survinger()

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

p2 <- ggplot2::ggplot(combined, ggplot2::aes(
  x = .data$time_short, y = .data$prevalence,
  color = .data$Method, fill = .data$Method, group = .data$Method
)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
                       alpha = 0.12, color = NA) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_fill_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x*100,1), "%"), limits = c(0, NA)) +
  ggplot2::labs(title = paste0("Design weighting changes ", target_lin, " prevalence estimates"),
                x = "Epiweek", y = "Estimated prevalence", color = NULL, fill = NULL,
                caption = "Shaded bands: 95% CI | Data: ECDC") +
  theme_survinger() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                 legend.position = "top")

grDevices::png(file.path(fig_dir, "fig2_compare.png"), width = 3000, height = 1800, res = 300)
print(p2); grDevices::dev.off()
cat("OK fig2\n")

# FIG 3: Allocation comparison
strats <- c("equal", "proportional", "min_mse")
alloc_details <- purrr::map_dfr(strats, function(s) {
  H <- design$n_strata
  if (s == "equal") {
    n <- rep(1000L %/% H, H)
    n[1] <- n[1] + 1000L - sum(n)
    tibble::tibble(region = design$strata_info$region, n = n, strategy = s)
  } else if (s == "proportional") {
    pi_h <- design$population$pop_total / sum(design$population$pop_total)
    n <- as.integer(round(pi_h * 1000))
    n[1] <- n[1] + 1000L - sum(n)
    tibble::tibble(region = design$population$region, n = n, strategy = s)
  } else {
    res <- surv_optimize_allocation(design, s, total_capacity = 1000)
    tibble::tibble(region = res$allocation$region, n = res$allocation$n_allocated, strategy = s)
  }
})
alloc_details$strategy <- factor(alloc_details$strategy,
  levels = c("equal", "proportional", "min_mse"),
  labels = c("Equal", "Proportional", "MSE-optimal"))

p3 <- ggplot2::ggplot(alloc_details, ggplot2::aes(x = .data$region, y = .data$n, fill = .data$strategy)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::scale_fill_manual(values = c(cols[["neutral"]], cols[["light_blue"]], cols[["primary"]])) +
  ggplot2::labs(title = "MSE-optimal allocation differs from equal or proportional",
                x = NULL, y = "Sequences allocated", fill = "Strategy",
                caption = "Total capacity: 1,000") +
  theme_survinger() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

grDevices::png(file.path(fig_dir, "fig3_allocation.png"), width = 2800, height = 1600, res = 300)
print(p3); grDevices::dev.off()
cat("OK fig3\n")

# FIG 4: Delay PMF + CDF
params <- delay$parameters[1, ]
d_seq <- 0:45
pmf_vals <- stats::dnbinom(d_seq, mu = params$mu, size = params$size)
cdf_vals <- stats::pnbinom(d_seq, mu = params$mu, size = params$size)
delay_df <- tibble::tibble(delay = d_seq, pmf = pmf_vals, cdf = cdf_vals)

p4 <- ggplot2::ggplot(delay_df) +
  ggplot2::geom_col(ggplot2::aes(x = .data$delay, y = .data$pmf),
                    fill = cols[["quaternary"]], alpha = 0.7) +
  ggplot2::geom_line(ggplot2::aes(x = .data$delay, y = .data$cdf * max(.data$pmf)),
                     color = cols[["primary"]], linewidth = 1) +
  ggplot2::scale_y_continuous(name = "Probability (PMF)",
    sec.axis = ggplot2::sec_axis(~ . / max(pmf_vals), name = "Cumulative (CDF)",
                                 labels = function(x) paste0(round(x*100), "%"))) +
  ggplot2::labs(title = "Reporting delay distribution (right-truncation corrected)",
                subtitle = paste0("NegBin fit: mu=", round(params$mu, 1), ", size=", round(params$size, 1)),
                x = "Delay (days)", caption = "Orange: PMF | Blue: CDF | Data: ECDC") +
  theme_survinger()

grDevices::png(file.path(fig_dir, "fig4_delay.png"), width = 2800, height = 1600, res = 300)
print(p4); grDevices::dev.off()
cat("OK fig4\n")

# FIG 5: Nowcast
nc_df <- nc$estimates
nc_df$time_short <- gsub("2023-", "", nc_df$time)
nc_df$type <- ifelse(nc_df$is_nowcast, "Nowcasted", "Complete")

p5 <- ggplot2::ggplot(nc_df, ggplot2::aes(x = .data$time_short)) +
  ggplot2::geom_col(ggplot2::aes(y = .data$n_observed), fill = cols[["neutral"]], alpha = 0.5, width = 0.6) +
  ggplot2::geom_ribbon(ggplot2::aes(y = .data$n_estimated, ymin = .data$ci_lower,
                                    ymax = .data$ci_upper, group = 1),
                       fill = cols[["quaternary"]], alpha = 0.15) +
  ggplot2::geom_line(ggplot2::aes(y = .data$n_estimated, group = 1),
                     color = cols[["quaternary"]], linewidth = 1) +
  ggplot2::geom_point(ggplot2::aes(y = .data$n_estimated, shape = .data$type),
                      color = cols[["quaternary"]], size = 2.5) +
  ggplot2::scale_shape_manual(values = c("Complete" = 16, "Nowcasted" = 17)) +
  ggplot2::labs(title = paste0("Nowcasting corrects right-truncation in ", target_lin),
                x = "Epiweek", y = "Sequence count", shape = "Status",
                caption = "Grey: reported | Orange: delay-corrected | Data: ECDC") +
  theme_survinger() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))

grDevices::png(file.path(fig_dir, "fig5_nowcast.png"), width = 3000, height = 1800, res = 300)
print(p5); grDevices::dev.off()
cat("OK fig5\n")

# FIG 6: Triple comparison
adj_df <- adj$estimates[!is.na(adj$estimates$prevalence), ]
adj_df$Method <- "Design + delay adjusted"
adj_df$time_short <- gsub("2023-", "", adj_df$time)
w_df2 <- weighted$estimates[!is.na(weighted$estimates$prevalence), ]
w_df2$Method <- "Design-weighted only"
w_df2$time_short <- gsub("2023-", "", w_df2$time)
n_df2 <- naive$estimates[!is.na(naive$estimates$prevalence), ]
n_df2$Method <- "Naive"
n_df2$time_short <- gsub("2023-", "", n_df2$time)

common <- c("time_short", "prevalence", "ci_lower", "ci_upper", "Method")
triple <- dplyr::bind_rows(
  adj_df[intersect(names(adj_df), common)],
  w_df2[intersect(names(w_df2), common)],
  n_df2[intersect(names(n_df2), common)]
)
triple$Method <- factor(triple$Method,
  levels = c("Naive", "Design-weighted only", "Design + delay adjusted"))

p6 <- ggplot2::ggplot(triple, ggplot2::aes(
  x = .data$time_short, y = .data$prevalence, color = .data$Method, group = .data$Method
)) +
  ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 1.8) +
  ggplot2::scale_color_manual(values = c(cols[["secondary"]], cols[["primary"]], cols[["tertiary"]])) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x*100,1), "%"), limits = c(0, NA)) +
  ggplot2::labs(title = paste0("Three levels of correction for ", target_lin),
                x = "Epiweek", y = "Prevalence", color = NULL,
                caption = "Red: naive | Blue: design-weighted | Green: + delay") +
  theme_survinger() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                 legend.position = "top")

grDevices::png(file.path(fig_dir, "fig6_triple.png"), width = 3000, height = 1800, res = 300)
print(p6); grDevices::dev.off()
cat("OK fig6\n")

# FIG 7: Bias heatmap
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

p7 <- ggplot2::ggplot(per_rw, ggplot2::aes(x = .data$ew_short, y = .data$region, fill = .data$bias)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::scale_fill_gradient2(low = cols[["primary"]], mid = "#FFFFFF", high = cols[["secondary"]],
    midpoint = 0, labels = function(x) paste0(ifelse(x>0,"+",""), round(x*100,1), "pp"),
    name = "Bias") +
  ggplot2::labs(title = "Estimation bias varies by country and time",
                x = "Epiweek", y = NULL, caption = "Blue: underestimate | Red: overestimate") +
  theme_survinger() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
                 legend.position = "right")

grDevices::png(file.path(fig_dir, "fig7_bias_heatmap.png"), width = 3200, height = 1400, res = 300)
print(p7); grDevices::dev.off()
cat("OK fig7\n")

# FIG 8: Dashboard
report <- surv_report(design, target_lineage = target_lin)
metrics <- tibble::tibble(
  metric = c("Gini\ncoefficient", "Effective\nsample ratio", "Detection\nprobability", "Mean bias\n(pp)"),
  display = c(sprintf("%.3f", report$gini), sprintf("%.1f%%", report$eff_ratio*100),
              sprintf("%.1f%%", report$detection_prob*100), sprintf("%.1f pp", report$mean_bias*100)),
  quality = c(ifelse(report$gini>0.4,"poor","moderate"), ifelse(report$eff_ratio<0.5,"poor","moderate"),
              ifelse(report$detection_prob<0.5,"poor","good"), ifelse(report$mean_bias>0.05,"poor","moderate"))
)
metrics$metric <- factor(metrics$metric, levels = rev(metrics$metric))
cmap <- c(good = cols[["tertiary"]], moderate = cols[["quaternary"]], poor = cols[["secondary"]])

p8 <- ggplot2::ggplot(metrics, ggplot2::aes(x = .data$metric, y = 1, fill = .data$quality)) +
  ggplot2::geom_tile(width = 0.85, height = 0.85, color = "white", linewidth = 2) +
  ggplot2::geom_text(ggplot2::aes(label = .data$display), size = 6, fontface = "bold", color = "white") +
  ggplot2::scale_fill_manual(values = cmap, guide = "none") +
  ggplot2::coord_flip() +
  ggplot2::labs(title = "Surveillance System Health Dashboard",
                subtitle = "Green = good | Orange = moderate | Red = concern") +
  theme_survinger() +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                 axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                 panel.border = ggplot2::element_blank())

grDevices::png(file.path(fig_dir, "fig8_dashboard.png"), width = 2400, height = 1200, res = 300)
print(p8); grDevices::dev.off()
cat("OK fig8\n")

figs <- list.files(fig_dir, pattern = "\\.png$")
cat(sprintf("\n=== %d publication-quality figures generated ===\n", length(figs)))
for (f in figs) {
  sz <- round(file.info(file.path(fig_dir, f))$size / 1024)
  cat(sprintf("  %s (%d KB)\n", f, sz))
}
