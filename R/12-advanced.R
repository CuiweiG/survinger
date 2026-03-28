# ============================================================
# Advanced features for ecosystem integration
# ============================================================

#' Estimate prevalence by subgroup
#'
#' Applies [surv_lineage_prevalence()] within subgroups defined by
#' a grouping variable. Analogous to [survey::svyby()] for stratified
#' survey analysis.
#'
#' @param design A `surv_design` object.
#' @param lineage Character. Target lineage.
#' @param by Character. Column name to group by (e.g., `"region"` or
#'   `"source_type"`).
#' @param time Character. Time aggregation. Default `"epiweek"`.
#' @param method Character. Estimation method. Default `"hajek"`.
#' @param conf_level Numeric. Default 0.95.
#'
#' @return A tibble with columns: group, time, lineage, prevalence, se,
#'   ci_lower, ci_upper, n_obs, effective_n.
#'
#' @seealso [surv_lineage_prevalence()], [surv_filter()]
#'
#' @examples
#' sim <- surv_simulate(n_regions = 4, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")],
#'                  sim$population, source_type = "source_type")
#' surv_prevalence_by(d, "BA.2.86", by = "region")
#'
#' @export
surv_prevalence_by <- function(design,
                               lineage,
                               by,
                               time = "epiweek",
                               method = "hajek",
                               conf_level = 0.95) {
  .assert_surv_design(design)
  checkmate::assert_string(lineage)
  checkmate::assert_string(by)
  .assert_columns_exist(design$data, by)

  groups <- sort(unique(design$data[[by]]))

  purrr::map_dfr(groups, function(g) {
    sub_design <- tryCatch(
      surv_filter(design, .data[[by]] == g),
      error = function(e) NULL
    )
    if (is.null(sub_design) || sub_design$n_obs < 5L) {
      return(tibble::tibble(
        group = g, time = NA_character_, lineage = lineage,
        prevalence = NA_real_, se = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        n_obs = 0L, effective_n = NA_real_
      ))
    }

    prev <- surv_lineage_prevalence(sub_design, lineage, time, method,
                                    conf_level, min_obs = 1L)
    est <- prev$estimates
    est$group <- g
    est[c("group", "time", "lineage", "prevalence", "se",
          "ci_lower", "ci_upper", "n_obs", "effective_n")]
  })
}


#' Sensitivity analysis across methods
#'
#' Runs all three prevalence estimators and delay/no-delay variants
#' on the same design, producing a comparison table. Essential for
#' robustness checks in publications.
#'
#' @param design A `surv_design` object.
#' @param lineage Character. Target lineage.
#' @param delay_fit Optional `surv_delay_fit` object. If provided,
#'   includes delay-adjusted estimates.
#' @param time Character. Default `"epiweek"`.
#' @param conf_level Numeric. Default 0.95.
#'
#' @return A tibble with one row per method-time combination, columns:
#'   method, time, prevalence, se, ci_lower, ci_upper.
#'
#' @seealso [surv_lineage_prevalence()], [surv_adjusted_prevalence()]
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' surv_sensitivity(d, "BA.2.86")
#'
#' @export
surv_sensitivity <- function(design, lineage,
                             delay_fit = NULL,
                             time = "epiweek",
                             conf_level = 0.95) {
  .assert_surv_design(design)
  checkmate::assert_string(lineage)

  methods <- c("naive", "horvitz_thompson", "hajek", "poststratified")
  results <- purrr::map_dfr(methods, function(m) {
    if (m == "naive") {
      est <- surv_naive_prevalence(design, lineage, time, conf_level)$estimates
    } else {
      est <- surv_lineage_prevalence(design, lineage, time, m,
                                     conf_level)$estimates
    }
    est$method <- m
    est[c("method", "time", "prevalence", "se", "ci_lower", "ci_upper",
          "n_obs", "effective_n")]
  })

  # Add delay-adjusted if delay_fit provided
  if (!is.null(delay_fit)) {
    adj <- surv_adjusted_prevalence(design, delay_fit, lineage, time,
                                   conf_level = conf_level)$estimates
    adj$method <- "adjusted"
    adj_cols <- intersect(names(results), names(adj))
    results <- dplyr::bind_rows(results, adj[adj_cols])
  }

  results$method <- factor(results$method,
    levels = c("naive", "horvitz_thompson", "hajek", "poststratified", "adjusted"))

  results
}


#' Compute power curve for detection across prevalence range
#'
#' Generates a detection probability curve that can be directly
#' plotted or included in publications. Answers: "At what prevalence
#' does our surveillance achieve X% detection?"
#'
#' @param design A `surv_design` object.
#' @param prevalence_range Numeric vector of prevalences to evaluate.
#'   Default `seq(0.001, 0.05, by = 0.001)`.
#' @param delay_fit Optional `surv_delay_fit`.
#' @param thresholds Numeric vector of detection thresholds to mark.
#'   Default `c(0.5, 0.8, 0.95)`.
#'
#' @return A list with:
#' \describe{
#'   \item{curve}{Tibble with prevalence and detection columns.}
#'   \item{thresholds}{Tibble with threshold, prevalence_needed columns.}
#' }
#'
#' @seealso [surv_detection_probability()], [surv_required_sequences()]
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' pc <- surv_power_curve(d)
#' pc$thresholds
#'
#' @export
surv_power_curve <- function(design,
                             prevalence_range = seq(0.001, 0.05, by = 0.001),
                             delay_fit = NULL,
                             thresholds = c(0.50, 0.80, 0.95)) {
  .assert_surv_design(design)

  det_probs <- vapply(prevalence_range, function(p) {
    surv_detection_probability(design, p, delay_fit = delay_fit)$overall
  }, numeric(1))

  curve <- tibble::tibble(
    prevalence = prevalence_range,
    detection = det_probs
  )

  # Find threshold crossings
  thresh_df <- purrr::map_dfr(thresholds, function(th) {
    idx <- which(det_probs >= th)[1]
    tibble::tibble(
      threshold = th,
      prevalence_needed = if (!is.na(idx)) prevalence_range[idx] else NA_real_
    )
  })

  structure(
    list(curve = curve, thresholds = thresh_df, design = design),
    class = "surv_power_curve"
  )
}


#' @param x A `surv_power_curve` object.
#' @param ... Additional arguments (unused).
#' @return A ggplot2 object.
#' @rdname surv_power_curve
#' @export
plot.surv_power_curve <- function(x, ...) {
  cols <- .surv_colors()
  curve <- x$curve
  thresh <- x$thresholds

  p <- ggplot2::ggplot(curve, ggplot2::aes(
    x = .data$prevalence * 100, y = .data$detection * 100)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = .data$detection * 100),
                         fill = cols[["light_blue"]], alpha = 0.25) +
    ggplot2::geom_line(color = cols[["primary"]], linewidth = 1.2)

  for (i in seq_len(nrow(thresh))) {
    th <- thresh$threshold[i] * 100
    pn <- thresh$prevalence_needed[i]
    if (!is.na(pn)) {
      p <- p +
        ggplot2::geom_hline(yintercept = th, linetype = "dashed",
                            color = cols[["neutral"]], linewidth = 0.4) +
        ggplot2::annotate("point", x = pn * 100, y = th, size = 3,
                          color = cols[["secondary"]]) +
        ggplot2::annotate("label", x = pn * 100, y = th - 5,
                          label = paste0(th, "% at ", round(pn * 100, 2), "%"),
                          size = 2.8, fill = "#FFFFFFDD",
                          color = cols[["secondary"]])
    }
  }

  p + ggplot2::scale_x_continuous(labels = function(v) paste0(v, "%")) +
    ggplot2::scale_y_continuous(labels = function(v) paste0(v, "%"),
                                breaks = c(0, 25, 50, 80, 95, 100)) +
    ggplot2::labs(title = "Variant detection power curve",
                  x = "True variant prevalence",
                  y = "Weekly detection probability") +
    theme_survinger(base_size = 12)
}

#' @export
print.surv_power_curve <- function(x, ...) {
  cli::cli_rule("Detection Power Curve")
  cli::cli_text("Prevalence range: {min(x$curve$prevalence)*100}% to {max(x$curve$prevalence)*100}%")
  cli::cli_text("")
  cli::cli_h3("Thresholds")
  for (i in seq_len(nrow(x$thresholds))) {
    th <- x$thresholds$threshold[i] * 100
    pn <- x$thresholds$prevalence_needed[i]
    if (!is.na(pn)) {
      pn_fmt <- round(pn * 100, 3)
      cli::cli_text("  {th}% detection requires {pn_fmt}% prevalence")
    }
  }
  invisible(x)
}
