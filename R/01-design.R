# ============================================================
# Core data container: surv_design
# ============================================================

#' Create a genomic surveillance design object
#'
#' Constructs a survey design object tailored for pathogen genomic
#' surveillance, encoding stratification structure, sequencing rates,
#' and population information needed for design-weighted inference.
#'
#' @param data Data frame of individual sequence records.
#' @param strata One-sided formula specifying stratification variables
#'   (e.g., `~ region` or `~ region + source_type`).
#' @param sequencing_rate Either a one-sided formula
#'   (`~ n_sequenced / n_positive`), a named numeric vector, or a
#'   data frame with strata columns and `seq_rate`.
#' @param population Data frame with one row per stratum, containing
#'   stratification variables and population-level denominators.
#' @param date_collected Column name for collection date.
#'   Default `"collection_date"`.
#' @param date_reported Column name for report date.
#'   Default `"report_date"`. Set `NULL` if unavailable.
#' @param lineage Column name for lineage. Default `"lineage"`.
#' @param source_type Column name for sample source. Default `NULL`.
#' @param source_config Optional tibble of per-source characteristics.
#'
#' @return An object of class `surv_design`.
#'
#' @seealso [surv_simulate()], [surv_lineage_prevalence()],
#'   [surv_optimize_allocation()]
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 42)
#' design <- surv_design(
#'   data = sim$sequences,
#'   strata = ~ region,
#'   sequencing_rate = sim$population[c("region", "seq_rate")],
#'   population = sim$population
#' )
#' print(design)
#'
#' @export
surv_design <- function(data,
                        strata,
                        sequencing_rate,
                        population,
                        date_collected = "collection_date",
                        date_reported = "report_date",
                        lineage = "lineage",
                        source_type = NULL,
                        source_config = NULL) {

  checkmate::assert_data_frame(data, min.rows = 1)
  checkmate::assert_formula(strata)
  checkmate::assert_data_frame(population, min.rows = 1)
  checkmate::assert_string(date_collected)
  checkmate::assert_string(date_reported, null.ok = TRUE)
  checkmate::assert_string(lineage)
  checkmate::assert_string(source_type, null.ok = TRUE)

  strata_vars <- all.vars(strata)
  .assert_columns_exist(data, c(strata_vars, date_collected, lineage))
  .assert_columns_exist(population, strata_vars)

  if (!is.null(date_reported)) .assert_columns_exist(data, date_reported)
  if (!is.null(source_type)) .assert_columns_exist(data, source_type)

  data[[date_collected]] <- as.Date(data[[date_collected]])
  if (!is.null(date_reported)) data[[date_reported]] <- as.Date(data[[date_reported]])

  # --- Compute strata info ---
  obs_counts <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(strata_vars))) |>
    dplyr::summarise(n_sequenced_obs = dplyr::n(), .groups = "drop")

  info <- dplyr::left_join(population, obs_counts, by = strata_vars)
  info$n_sequenced_obs[is.na(info$n_sequenced_obs)] <- 0L

  # --- Parse sequencing_rate ---
  if (inherits(sequencing_rate, "formula")) {
    rate_expr <- rlang::f_rhs(sequencing_rate)
    info$seq_rate <- rlang::eval_tidy(rate_expr, data = info)
  } else if (is.numeric(sequencing_rate)) {
    strata_key <- .make_strata_key(info, strata_vars)
    checkmate::assert_names(names(sequencing_rate),
                            must.include = unique(strata_key))
    info$seq_rate <- sequencing_rate[strata_key]
  } else if (is.data.frame(sequencing_rate)) {
    .assert_columns_exist(sequencing_rate, c(strata_vars, "seq_rate"))
    info <- dplyr::left_join(
      info |> dplyr::select(-dplyr::any_of("seq_rate")),
      sequencing_rate[c(strata_vars, "seq_rate")],
      by = strata_vars
    )
  } else {
    cli::cli_abort(
      "{.arg sequencing_rate} must be a formula, named numeric vector, or data frame."
    )
  }

  info$seq_rate <- .validate_seq_rates(info$seq_rate)
  info$stratum_id <- seq_len(nrow(info))

  # --- Compute weights ---
  weights <- tibble::tibble(
    stratum_id = info$stratum_id,
    weight = 1.0 / info$seq_rate,
    weight_normalized = (1.0 / info$seq_rate) / mean(1.0 / info$seq_rate)
  )

  # --- Add epiweek ---
  data[["epiweek"]] <- .date_to_epiweek(data[[date_collected]])

  structure(
    list(
      data = tibble::as_tibble(data),
      strata_formula = strata,
      strata_vars = strata_vars,
      strata_info = tibble::as_tibble(info),
      weights = weights,
      population = tibble::as_tibble(population),
      source_config = source_config,
      col_date_collected = date_collected,
      col_date_reported = date_reported,
      col_lineage = lineage,
      col_source_type = source_type,
      n_obs = nrow(data),
      n_strata = nrow(info),
      created = Sys.time()
    ),
    class = "surv_design"
  )
}


#' Update sequencing rates in a surveillance design
#'
#' @param design A `surv_design` object.
#' @param new_rates Data frame with strata columns + `seq_rate`,
#'   or named numeric vector.
#' @return Updated `surv_design`.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' new_r <- sim$population[c("region", "seq_rate")]
#' new_r$seq_rate <- new_r$seq_rate * 2
#' d2 <- surv_update_rates(d, new_r)
#'
#' @export
surv_update_rates <- function(design, new_rates) {
  .assert_surv_design(design)
  sv <- design$strata_vars

  if (is.data.frame(new_rates)) {
    .assert_columns_exist(new_rates, c(sv, "seq_rate"))
    design$strata_info <- design$strata_info |>
      dplyr::select(-"seq_rate") |>
      dplyr::left_join(new_rates[c(sv, "seq_rate")], by = sv)
  } else if (is.numeric(new_rates)) {
    strata_key <- .make_strata_key(design$strata_info, sv)
    design$strata_info$seq_rate <- new_rates[strata_key]
  } else {
    cli::cli_abort(
      "{.arg new_rates} must be a data frame or named numeric vector."
    )
  }

  design$strata_info$seq_rate <- .validate_seq_rates(design$strata_info$seq_rate)
  design$weights <- tibble::tibble(
    stratum_id = design$strata_info$stratum_id,
    weight = 1.0 / design$strata_info$seq_rate,
    weight_normalized = (1.0 / design$strata_info$seq_rate) /
      mean(1.0 / design$strata_info$seq_rate)
  )
  design
}


#' Override design weights with custom values
#'
#' @param design A `surv_design` object.
#' @param weights Numeric vector of length equal to number of strata.
#' @return Updated `surv_design`.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' d2 <- surv_set_weights(d, rep(1.0, d$n_strata))
#'
#' @export
surv_set_weights <- function(design, weights) {
  .assert_surv_design(design)
  checkmate::assert_numeric(weights, len = design$n_strata, lower = 0)
  design$weights$weight <- weights
  design$weights$weight_normalized <- weights / mean(weights)
  design
}


#' Subset a surveillance design by filter criteria
#'
#' Creates a new `surv_design` object containing only sequences
#' matching the specified filter criteria.
#'
#' @param design A `surv_design` object.
#' @param ... Filter expressions passed to [dplyr::filter()].
#'
#' @return A new `surv_design` object with filtered data.
#'
#' @seealso [surv_design()]
#'
#' @examples
#' sim <- surv_simulate(n_regions = 5, n_weeks = 12, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' d_sub <- surv_filter(d, region %in% c("Region_A", "Region_B"))
#' print(d_sub)
#'
#' @export
surv_filter <- function(design, ...) {
  .assert_surv_design(design)
  new_data <- dplyr::filter(design$data, ...)
  if (nrow(new_data) == 0L) {
    cli::cli_warn("Filter resulted in 0 observations.")
  }
  surv_design(
    data = new_data,
    strata = design$strata_formula,
    sequencing_rate = design$strata_info[c(design$strata_vars, "seq_rate")],
    population = design$population,
    date_collected = design$col_date_collected,
    date_reported = design$col_date_reported,
    lineage = design$col_lineage,
    source_type = design$col_source_type,
    source_config = design$source_config
  )
}
