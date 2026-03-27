# ============================================================
# Sequencing resource allocation optimization
# ============================================================

#' Optimize sequencing allocation across strata
#'
#' Given fixed total sequencing capacity, finds the optimal allocation
#' across strata that minimizes a specified objective function.
#'
#' @param design A `surv_design` object.
#' @param objective Character. One of `"min_mse"`, `"max_detection"`,
#'   or `"min_imbalance"`.
#' @param total_capacity Integer. Total sequences available.
#' @param budget Numeric or `NULL`. Optional budget constraint.
#' @param min_per_stratum Integer. Minimum per stratum. Default 2.
#' @param target_lineage Character. Required for `"max_detection"`.
#' @param target_prevalence Numeric. Assumed prevalence for detection.
#'   Default 0.01.
#' @param cost_col Character or `NULL`. Column name for per-sequence cost.
#'
#' @return A `surv_allocation` object.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 4, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' a <- surv_optimize_allocation(d, "min_mse", total_capacity = 500)
#' print(a)
#'
#' @export
surv_optimize_allocation <- function(design,
                                     objective = c("min_mse",
                                                   "max_detection",
                                                   "min_imbalance"),
                                     total_capacity,
                                     budget = NULL,
                                     min_per_stratum = 2L,
                                     target_lineage = NULL,
                                     target_prevalence = 0.01,
                                     cost_col = NULL) {
  .assert_surv_design(design)
  objective <- match.arg(objective)
  checkmate::assert_count(total_capacity, positive = TRUE)
  checkmate::assert_count(min_per_stratum)

  H <- design$n_strata
  info <- design$strata_info
  strata_vars <- design$strata_vars

  if (total_capacity < H * min_per_stratum) {
    cli::cli_abort(
      "Total capacity ({total_capacity}) < strata ({H}) x min_per_stratum ({min_per_stratum})."
    )
  }

  # Population shares
  if ("pop_total" %in% names(info)) {
    pi_h <- info$pop_total / sum(info$pop_total)
  } else if ("n_positive" %in% names(info)) {
    pi_h <- info$n_positive / sum(info$n_positive)
  } else {
    pi_h <- rep(1 / H, H)
  }

  # Estimate p_h from data
  p_h <- .estimate_stratum_prevalence(design, target_lineage)

  # Per-stratum cost
  if (!is.null(cost_col) && cost_col %in% names(info)) {
    c_h <- info[[cost_col]]
  } else {
    c_h <- rep(1, H)
  }

  # Solve
  alloc <- switch(objective,
    "min_mse"        = .solve_min_mse(pi_h, p_h, c_h, total_capacity, min_per_stratum, budget),
    "max_detection"  = .solve_max_detection(p_h, target_prevalence, total_capacity, min_per_stratum),
    "min_imbalance"  = .solve_min_imbalance(pi_h, total_capacity, min_per_stratum)
  )

  alloc_tbl <- info[strata_vars]
  alloc_tbl$n_allocated <- alloc$n
  alloc_tbl$proportion <- alloc$n / total_capacity

  new_surv_allocation(
    allocation = alloc_tbl,
    objective = objective,
    total_capacity = total_capacity,
    constraints = list(budget = budget, min_per_stratum = min_per_stratum),
    diagnostics = alloc$diagnostics
  )
}


#' Compare multiple allocation strategies
#'
#' @param design A `surv_design` object.
#' @param strategies Character vector. Default includes all built-in.
#' @param total_capacity Integer. Total sequences.
#' @param target_prevalence Numeric. For detection objective.
#' @param ... Passed to [surv_optimize_allocation()].
#'
#' @return A tibble comparing strategies.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 4, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' surv_compare_allocations(d, total_capacity = 200)
#'
#' @export
surv_compare_allocations <- function(design,
                                     strategies = c("equal", "proportional",
                                                    "min_mse", "max_detection",
                                                    "min_imbalance"),
                                     total_capacity,
                                     target_prevalence = 0.01,
                                     ...) {
  .assert_surv_design(design)
  H <- design$n_strata
  info <- design$strata_info

  if ("pop_total" %in% names(info)) {
    pi_h <- info$pop_total / sum(info$pop_total)
  } else {
    pi_h <- rep(1 / H, H)
  }
  p_h <- .estimate_stratum_prevalence(design, NULL)

  results <- purrr::map_dfr(strategies, function(s) {
    n_alloc <- if (s == "equal") {
      .equal_allocation(H, total_capacity)
    } else if (s == "proportional") {
      .proportional_allocation(pi_h, total_capacity)
    } else {
      a <- surv_optimize_allocation(design, s, total_capacity,
                                    target_prevalence = target_prevalence, ...)
      a$allocation$n_allocated
    }
    tibble::tibble(
      strategy = s,
      total_mse = .compute_allocation_mse(pi_h, p_h, n_alloc),
      detection_prob = 1 - prod((1 - target_prevalence)^n_alloc),
      imbalance = sum((n_alloc / sum(n_alloc) - pi_h)^2)
    )
  })
  results
}


# ---- Internal solvers ----

#' @keywords internal
.solve_min_mse <- function(pi_h, p_h, c_h, N, min_n, budget) {
  H <- length(pi_h)
  # Neyman initial: n_h* proportional to pi_h * sqrt(p*(1-p)) / sqrt(c)
  sigma_h <- sqrt(p_h * (1 - p_h) + 0.001)
  raw <- pi_h * sigma_h / sqrt(c_h)
  n_init <- pmax(round(raw / sum(raw) * N), min_n)
  n_init <- .integerize_allocation(n_init, N, min_n)
  list(n = n_init, diagnostics = list(method = "neyman_rounded"))
}

#' @keywords internal
.solve_max_detection <- function(p_h, target_prev, N, min_n) {
  H <- length(p_h)
  # Allocate more to strata with lower prevalence (more "surprising")
  inv_p <- 1 / (p_h + target_prev + 0.001)
  raw <- inv_p / sum(inv_p) * N
  n_alloc <- pmax(round(raw), min_n)
  n_alloc <- .integerize_allocation(n_alloc, N, min_n)
  list(n = n_alloc, diagnostics = list(method = "detection_heuristic"))
}

#' @keywords internal
.solve_min_imbalance <- function(pi_h, N, min_n) {
  # Proportional allocation = minimum imbalance
  n_alloc <- pmax(round(pi_h * N), min_n)
  n_alloc <- .integerize_allocation(n_alloc, N, min_n)
  list(n = n_alloc, diagnostics = list(method = "proportional_rounded"))
}


# ---- Allocation helpers ----

#' @keywords internal
.equal_allocation <- function(H, N) {
  base <- rep(N %/% H, H)
  remainder <- N - sum(base)
  if (remainder > 0) base[seq_len(remainder)] <- base[seq_len(remainder)] + 1L
  as.integer(base)
}

#' @keywords internal
.proportional_allocation <- function(pi_h, N) {
  raw <- round(pi_h * N)
  .integerize_allocation(raw, N, min_n = 0L)
}

#' @keywords internal
.integerize_allocation <- function(n, total, min_n = 0L) {
  n <- pmax(as.integer(round(n)), as.integer(min_n))
  diff_val <- total - sum(n)
  if (diff_val > 0) {
    # Add to largest strata first
    idx <- order(n, decreasing = TRUE)
    for (i in seq_len(abs(diff_val))) {
      n[idx[(i - 1L) %% length(idx) + 1L]] <- n[idx[(i - 1L) %% length(idx) + 1L]] + 1L
    }
  } else if (diff_val < 0) {
    idx <- order(n, decreasing = TRUE)
    for (i in seq_len(abs(diff_val))) {
      pos <- idx[(i - 1L) %% length(idx) + 1L]
      if (n[pos] > min_n) n[pos] <- n[pos] - 1L
    }
  }
  as.integer(n)
}

#' @keywords internal
.estimate_stratum_prevalence <- function(design, target_lineage) {
  lin_col <- design$col_lineage
  strata_vars <- design$strata_vars
  dat <- design$data

  if (is.null(target_lineage)) {
    # Use most common non-"Other" lineage
    lin_tab <- sort(table(dat[[lin_col]]), decreasing = TRUE)
    target_lineage <- names(lin_tab)[1]
  }

  prev_by_stratum <- dat |>
    dplyr::group_by(dplyr::across(dplyr::all_of(strata_vars))) |>
    dplyr::summarise(
      p_h = mean(.data[[lin_col]] == target_lineage),
      .groups = "drop"
    )

  # Match to strata_info order
  info <- design$strata_info
  merged <- dplyr::left_join(info[strata_vars], prev_by_stratum, by = strata_vars)
  merged$p_h[is.na(merged$p_h)] <- 0.01
  merged$p_h
}

#' @keywords internal
.compute_allocation_mse <- function(pi_h, p_h, n_h) {
  sum(pi_h^2 * p_h * (1 - p_h) / pmax(n_h, 1))
}
