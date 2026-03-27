# ============================================================
# S3 class constructors and print/summary methods
# ============================================================

# ---- surv_design (public constructor in 01-design.R) ----

#' @export
print.surv_design <- function(x, ...) {
  cli::cli_rule("Genomic Surveillance Design")
  n_obs_fmt <- formatC(x$n_obs, big.mark = ",")
  strata_label <- paste(x$strata_vars, collapse = " x ")
  cli::cli_text("Observations: {.strong {n_obs_fmt}}")
  cli::cli_text("Strata: {.strong {x$n_strata}} ({strata_label})")
  dr <- range(x$data[[x$col_date_collected]], na.rm = TRUE)
  cli::cli_text("Date range: {dr[1]} to {dr[2]}")
  n_lin <- length(unique(x$data[[x$col_lineage]]))
  cli::cli_text("Lineages: {.strong {n_lin}}")
  if (!is.null(x$col_source_type)) {
    src_label <- paste(sort(unique(x$data[[x$col_source_type]])), collapse = ", ")
    cli::cli_text("Sources: {src_label}")
  }
  wr <- range(x$weights$weight, na.rm = TRUE)
  wt_label <- paste0("[", round(wr[1], 3), ", ", round(wr[2], 3), "]")
  cli::cli_text("Weight range: {wt_label}")
  invisible(x)
}

#' @export
summary.surv_design <- function(object, ...) {
  out <- list(
    n_obs = object$n_obs,
    n_strata = object$n_strata,
    strata_vars = object$strata_vars,
    date_range = range(object$data[[object$col_date_collected]], na.rm = TRUE),
    lineage_counts = sort(table(object$data[[object$col_lineage]]),
                          decreasing = TRUE),
    weight_summary = summary(object$weights$weight)
  )
  structure(out, class = "summary.surv_design")
}

#' @export
print.summary.surv_design <- function(x, ...) {
  cli::cli_rule("Surveillance Design Summary")
  n_fmt <- formatC(x$n_obs, big.mark = ",")
  sv_label <- paste(x$strata_vars, collapse = " x ")
  cli::cli_text("Sequences: {.strong {n_fmt}}")
  cli::cli_text("Strata: {x$n_strata} ({sv_label})")
  cli::cli_text("Period: {x$date_range[1]} to {x$date_range[2]}")
  cli::cli_h3("Top lineages")
  top <- utils::head(x$lineage_counts, 5)
  for (i in seq_along(top)) cli::cli_text("  {names(top)[i]}: {top[i]}")
  cli::cli_h3("Sampling weights")
  print(x$weight_summary)
  invisible(x)
}

# ---- surv_allocation ----

#' Internal constructor for surv_allocation
#' @keywords internal
new_surv_allocation <- function(allocation, objective, total_capacity,
                                constraints, diagnostics) {
  structure(list(
    allocation = tibble::as_tibble(allocation),
    objective = objective,
    total_capacity = total_capacity,
    constraints = constraints,
    diagnostics = diagnostics
  ), class = "surv_allocation")
}

#' @export
print.surv_allocation <- function(x, ...) {
  cli::cli_rule("Optimal Sequencing Allocation")
  cli::cli_text("Objective: {.strong {x$objective}}")
  cli::cli_text("Total capacity: {.strong {x$total_capacity}} sequences")
  cli::cli_text("Strata: {nrow(x$allocation)}")
  cli::cli_text("")
  print(x$allocation, n = 10)
  invisible(x)
}

# ---- surv_prevalence ----

#' Internal constructor for surv_prevalence
#' @keywords internal
new_surv_prevalence <- function(estimates, design, method, lineage,
                                conf_level, time_unit) {
  structure(list(
    estimates = tibble::as_tibble(estimates),
    design = design, method = method, lineage = lineage,
    conf_level = conf_level, time_unit = time_unit
  ), class = "surv_prevalence")
}

#' @export
print.surv_prevalence <- function(x, ...) {
  cli::cli_rule("Lineage Prevalence Estimate")
  cli::cli_text("Lineage: {.val {x$lineage}}")
  cli::cli_text("Method: {.val {x$method}}")
  cli::cli_text("Confidence level: {x$conf_level}")
  cli::cli_text("Time periods: {nrow(x$estimates)}")
  cli::cli_text("")
  print(x$estimates, n = 10)
  invisible(x)
}

# ---- surv_delay_fit ----

#' Internal constructor for surv_delay_fit
#' @keywords internal
new_surv_delay_fit <- function(distribution, parameters, strata,
                               data_summary, diagnostics) {
  structure(list(
    distribution = distribution, parameters = parameters,
    strata = strata, data_summary = data_summary,
    diagnostics = diagnostics
  ), class = "surv_delay_fit")
}

#' @export
print.surv_delay_fit <- function(x, ...) {
  cli::cli_rule("Reporting Delay Distribution")
  cli::cli_text("Distribution: {.val {x$distribution}}")
  sl <- if (is.null(x$strata)) "none (pooled)" else paste(x$strata, collapse = " x ")
  cli::cli_text("Strata: {sl}")
  cli::cli_text("Observations: {x$data_summary$n}")
  cli::cli_text("Mean delay: {round(x$data_summary$mean_delay, 1)} days")
  cli::cli_text("")
  print(x$parameters)
  invisible(x)
}

# ---- surv_nowcast ----

#' Internal constructor for surv_nowcast
#' @keywords internal
new_surv_nowcast <- function(estimates, delay_fit, truncation_window,
                             method, lineage) {
  structure(list(
    estimates = tibble::as_tibble(estimates),
    delay_fit = delay_fit, truncation_window = truncation_window,
    method = method, lineage = lineage
  ), class = "surv_nowcast")
}

#' @export
print.surv_nowcast <- function(x, ...) {
  cli::cli_rule("Delay-Adjusted Nowcast")
  cli::cli_text("Method: {.val {x$method}}")
  cli::cli_text("Truncation window: {x$truncation_window} periods")
  if (!is.null(x$lineage)) cli::cli_text("Lineage: {.val {x$lineage}}")
  cli::cli_text("")
  print(x$estimates, n = 10)
  invisible(x)
}

# ---- surv_adjusted ----

#' Internal constructor for surv_adjusted
#' @keywords internal
new_surv_adjusted <- function(estimates, design_component, delay_component,
                              method) {
  structure(list(
    estimates = tibble::as_tibble(estimates),
    design_component = design_component,
    delay_component = delay_component,
    method = method
  ), class = "surv_adjusted")
}

#' @export
print.surv_adjusted <- function(x, ...) {
  cli::cli_rule("Design-Weighted Delay-Adjusted Prevalence")
  cli::cli_text("Correction: {.val {x$method}}")
  cli::cli_text("")
  print(x$estimates, n = 10)
  invisible(x)
}
