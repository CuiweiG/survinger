# ============================================================
# Diagnostic and comparison tools
# ============================================================

#' Compare weighted vs naive prevalence estimates
#'
#' Side-by-side plot showing the impact of design correction.
#'
#' @param weighted A `surv_prevalence` object (design-weighted).
#' @param naive A `surv_prevalence` object (unweighted).
#' @param title Character or `NULL`. Plot title.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' w <- surv_lineage_prevalence(d, "BA.2.86")
#' n <- surv_naive_prevalence(d, "BA.2.86")
#' surv_compare_estimates(w, n)
#'
#' @export
surv_compare_estimates <- function(weighted, naive, title = NULL) {
  .assert_surv_class(weighted, "surv_prevalence")
  .assert_surv_class(naive, "surv_prevalence")

  if (is.null(title)) title <- paste0("Prevalence: ", weighted$lineage)

  w_df <- weighted$estimates
  w_df$method <- paste0("Weighted (", weighted$method, ")")
  n_df <- naive$estimates
  n_df$method <- "Naive (unweighted)"
  combined <- dplyr::bind_rows(w_df, n_df)
  combined <- combined[!is.na(combined$prevalence), , drop = FALSE]

  ggplot2::ggplot(combined, ggplot2::aes(
    x = .data$time, y = .data$prevalence,
    color = .data$method, fill = .data$method, group = .data$method
  )) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      alpha = 0.15, color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(x * 100, 1), "%"),
      limits = c(0, NA)
    ) +
    ggplot2::scale_color_manual(values = c("#D4652F", "#2E75B6")) +
    ggplot2::scale_fill_manual(values = c("#D4652F", "#2E75B6")) +
    ggplot2::labs(title = title, x = "Time", y = "Prevalence",
                  color = "Method", fill = "Method") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Plot sequencing rate inequality across strata
#'
#' @param design A `surv_design` object.
#' @return A ggplot2 object.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 5, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' surv_plot_sequencing_rates(d)
#'
#' @export
surv_plot_sequencing_rates <- function(design) {
  .assert_surv_design(design)
  info <- design$strata_info
  info$label <- .make_strata_key(info, design$strata_vars)
  info <- info[order(info$seq_rate), , drop = FALSE]
  info$label <- factor(info$label, levels = info$label)

  ggplot2::ggplot(info, ggplot2::aes(x = .data$label, y = .data$seq_rate)) +
    ggplot2::geom_col(fill = "#2E75B6") +
    ggplot2::geom_hline(yintercept = mean(info$seq_rate),
                        linetype = "dashed", color = "#D4652F") +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(x * 100, 1), "%")
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Sequencing Rates by Stratum",
                  subtitle = "Dashed line = mean rate",
                  x = NULL, y = "Sequencing Rate") +
    ggplot2::theme_minimal(base_size = 12)
}


#' Plot allocation plan
#'
#' @param allocation A `surv_allocation` object.
#' @return A ggplot2 object.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 5, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' a <- surv_optimize_allocation(d, "min_mse", total_capacity = 500)
#' surv_plot_allocation(a)
#'
#' @export
surv_plot_allocation <- function(allocation) {
  .assert_surv_class(allocation, "surv_allocation")
  alloc <- allocation$allocation
  meta_cols <- setdiff(names(alloc), c("n_allocated", "proportion"))
  alloc$label <- do.call(paste, c(as.list(alloc[meta_cols]), list(sep = ":")))
  alloc <- alloc[order(alloc$n_allocated), , drop = FALSE]
  alloc$label <- factor(alloc$label, levels = alloc$label)

  ggplot2::ggplot(alloc, ggplot2::aes(x = .data$label, y = .data$n_allocated)) +
    ggplot2::geom_col(fill = "#5B8C5A") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Allocation (", allocation$objective, ")"),
      subtitle = paste0("Total capacity: ", allocation$total_capacity),
      x = NULL, y = "Sequences Allocated"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Compute design effect over time
#'
#' @param weighted A `surv_prevalence` object.
#' @param naive A `surv_prevalence` object.
#' @return A tibble with time, deff, and bias_correction columns.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' w <- surv_lineage_prevalence(d, "BA.2.86")
#' n <- surv_naive_prevalence(d, "BA.2.86")
#' surv_design_effect(w, n)
#'
#' @export
surv_design_effect <- function(weighted, naive) {
  .assert_surv_class(weighted, "surv_prevalence")
  .assert_surv_class(naive, "surv_prevalence")

  dplyr::inner_join(
    weighted$estimates[c("time", "prevalence", "se")] |>
      dplyr::rename(prev_w = "prevalence", se_w = "se"),
    naive$estimates[c("time", "prevalence", "se")] |>
      dplyr::rename(prev_n = "prevalence", se_n = "se"),
    by = "time"
  ) |>
    dplyr::mutate(
      deff = (.data$se_w / pmax(.data$se_n, 1e-10))^2,
      bias_correction = .data$prev_w - .data$prev_n
    )
}
