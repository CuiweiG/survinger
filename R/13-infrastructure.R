# ============================================================
# Infrastructure features — making survinger a foundation package
# Inspired by: dplyr verbs, tidyverse pipes, knitr integration
# ============================================================

# ---- Pipe-friendly entry point ----

#' Pipe-friendly surveillance analysis
#'
#' Convenience wrapper that creates a design, estimates prevalence,
#' and optionally applies delay correction in a single pipe-friendly call.
#' Designed for rapid exploratory analysis in interactive sessions.
#'
#' @param data Data frame of sequence records.
#' @param strata One-sided formula for stratification.
#' @param sequencing_rate Sequencing rate specification (see [surv_design()]).
#' @param population Population data frame.
#' @param lineage Character. Target lineage to estimate.
#' @param correct_delay Logical. Apply delay correction? Default `FALSE`.
#' @param method Character. Prevalence method. Default `"hajek"`.
#' @param ... Additional arguments passed to [surv_design()].
#'
#' @return A `surv_prevalence` or `surv_adjusted` object.
#'
#' @seealso [surv_design()], [surv_lineage_prevalence()]
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' # One-liner analysis:
#' result <- surv_estimate(
#'   data = sim$sequences, strata = ~ region,
#'   sequencing_rate = sim$population[c("region", "seq_rate")],
#'   population = sim$population,
#'   lineage = "BA.2.86"
#' )
#' print(result)
#'
#' @export
surv_estimate <- function(data, strata, sequencing_rate, population,
                          lineage, correct_delay = FALSE,
                          method = "hajek", ...) {
  design <- surv_design(data = data, strata = strata,
                        sequencing_rate = sequencing_rate,
                        population = population, ...)

  if (correct_delay) {
    if (is.null(design$col_date_reported)) {
      cli::cli_warn("No report date column. Skipping delay correction.")
      return(surv_lineage_prevalence(design, lineage, method = method))
    }
    delay_fit <- surv_estimate_delay(design)
    surv_adjusted_prevalence(design, delay_fit, lineage,
                             prevalence_method = method)
  } else {
    surv_lineage_prevalence(design, lineage, method = method)
  }
}


# ---- Combine prevalence results ----

#' Combine multiple prevalence estimates
#'
#' Bind results from different lineages or methods into a single
#' tibble for comparison plots and tables.
#'
#' @param ... `surv_prevalence` objects to combine.
#'
#' @return A tibble with a `source` column identifying each input.
#'
#' @seealso [surv_lineage_prevalence()], [surv_sensitivity()]
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' p1 <- surv_lineage_prevalence(d, "BA.5")
#' p2 <- surv_lineage_prevalence(d, "XBB.1.5")
#' combined <- surv_bind(p1, p2)
#' head(combined)
#'
#' @export
surv_bind <- function(...) {
  objects <- list(...)
  purrr::imap_dfr(objects, function(obj, i) {
    if (inherits(obj, "surv_prevalence")) {
      est <- obj$estimates
      est$source <- paste0(obj$lineage, " (", obj$method, ")")
      est
    } else if (inherits(obj, "surv_adjusted")) {
      est <- obj$estimates
      est$source <- paste0(est$lineage[1], " (adjusted)")
      est
    } else {
      cli::cli_warn("Object {i} is not a surv_prevalence or surv_adjusted. Skipping.")
      tibble::tibble()
    }
  })
}


# ---- Automatic knitr printing ----

#' Format prevalence results for knitr tables
#'
#' Produces a publication-ready summary table from any survinger result.
#' Automatically called when objects are printed inside RMarkdown chunks
#' if `knitr` is loaded.
#'
#' @param x A survinger result object.
#' @param digits Integer. Decimal places for prevalence. Default 3.
#' @param percent Logical. Display as percentages? Default `TRUE`.
#' @param ... Additional arguments (unused).
#'
#' @return A tibble formatted for display.
#'
#' @seealso [surv_lineage_prevalence()], [tidy.surv_prevalence()]
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' prev <- surv_lineage_prevalence(d, "BA.2.86")
#' surv_table(prev)
#'
#' @export
surv_table <- function(x, digits = 3, percent = TRUE, ...) {
  if (inherits(x, "surv_prevalence") || inherits(x, "surv_adjusted")) {
    est <- if (inherits(x, "surv_prevalence")) x$estimates else x$estimates
    if (percent) {
      est$prevalence <- round(est$prevalence * 100, digits - 1)
      est$ci_lower <- round(est$ci_lower * 100, digits - 1)
      est$ci_upper <- round(est$ci_upper * 100, digits - 1)
      est$se <- round(est$se * 100, digits)
    } else {
      est$prevalence <- round(est$prevalence, digits)
      est$ci_lower <- round(est$ci_lower, digits)
      est$ci_upper <- round(est$ci_upper, digits)
      est$se <- round(est$se, digits)
    }
    # Add formatted CI column
    est$ci <- paste0("[", est$ci_lower, ", ", est$ci_upper, "]")
    # Select display columns
    display_cols <- intersect(c("time", "lineage", "prevalence", "ci", "se",
                                "n_obs", "effective_n"), names(est))
    est[display_cols]
  } else if (inherits(x, "surv_allocation")) {
    x$allocation
  } else if (inherits(x, "surv_delay_fit")) {
    x$parameters
  } else {
    cli::cli_abort("surv_table() does not support {.cls {class(x)[1]}} objects.")
  }
}


# ---- Design summary statistics ----

#' Compute surveillance quality metrics
#'
#' Returns a single-row tibble of design quality indicators suitable
#' for inclusion in manuscripts. Analogous to `broom::glance()` but
#' for the surveillance design itself.
#'
#' @param design A `surv_design` object.
#' @param target_lineage Character or `NULL`. Default `NULL` auto-selects.
#' @param target_prevalence Numeric. For detection calculation. Default 0.01.
#'
#' @return A single-row tibble with columns: n_obs, n_strata, gini,
#'   rate_ratio, effective_n, deff, detection_prob, mean_bias.
#'
#' @seealso [surv_report()], [surv_design()]
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' surv_quality(d)
#'
#' @export
surv_quality <- function(design, target_lineage = NULL,
                         target_prevalence = 0.01) {
  .assert_surv_design(design)
  info <- design$strata_info
  dat <- design$data

  if (is.null(target_lineage)) {
    lin_tab <- sort(table(dat[[design$col_lineage]]), decreasing = TRUE)
    non_other <- names(lin_tab)[!grepl("^other$", names(lin_tab), ignore.case = TRUE)]
    target_lineage <- if (length(non_other) > 0) non_other[1] else names(lin_tab)[1]
  }

  rate_range <- range(info$seq_rate)
  gini <- .gini_coefficient(info$seq_rate)

  w_dat <- .map_weights_to_obs(design)
  eff_n <- .kish_eff_n(w_dat$weight)

  det <- surv_detection_probability(design, target_prevalence)

  weighted <- surv_lineage_prevalence(design, target_lineage, method = "hajek")
  naive <- surv_naive_prevalence(design, target_lineage)
  bias <- mean(abs(weighted$estimates$prevalence - naive$estimates$prevalence),
               na.rm = TRUE)

  tibble::tibble(
    n_obs = design$n_obs,
    n_strata = design$n_strata,
    n_lineages = length(unique(dat[[design$col_lineage]])),
    gini = round(gini, 3),
    rate_min = round(rate_range[1], 4),
    rate_max = round(rate_range[2], 4),
    rate_ratio = round(rate_range[2] / rate_range[1], 1),
    effective_n = round(eff_n),
    deff = round(design$n_obs / eff_n, 3),
    detection_at_1pct = round(det$overall, 3),
    mean_abs_bias_pp = round(bias * 100, 2),
    target_lineage = target_lineage
  )
}
