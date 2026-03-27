#' Create strata key string from data frame rows
#' @param df Data frame.
#' @param strata_vars Character vector of column names.
#' @return Character vector of strata identifiers.
#' @keywords internal
.make_strata_key <- function(df, strata_vars) {
  if (length(strata_vars) == 1L) {
    return(as.character(df[[strata_vars]]))
  }
  do.call(paste, c(as.list(df[strata_vars]), list(sep = ":")))
}

#' Map stratum-level weights to individual observations
#' @param design A surv_design object.
#' @return Tibble with all data columns plus stratum_id, seq_rate, weight.
#' @keywords internal
.map_weights_to_obs <- function(design) {
  strata_vars <- design$strata_vars
  weight_info <- design$strata_info[c(strata_vars, "stratum_id", "seq_rate")]
  dat <- dplyr::left_join(design$data, weight_info, by = strata_vars,
                          relationship = "many-to-one")
  dat[["weight"]] <- 1.0 / dat[["seq_rate"]]
  dat
}

#' Kish effective sample size
#' @param weights Numeric vector of positive weights.
#' @return Scalar effective sample size.
#' @keywords internal
.kish_eff_n <- function(weights) {
  (sum(weights))^2 / sum(weights^2)
}

#' Generate unequal shares summing to 1
#' @param n Number of shares.
#' @param shape Gamma shape parameter.
#' @return Numeric vector summing to 1.
#' @keywords internal
.generate_unequal_shares <- function(n, shape = 2) {
  shares <- stats::rgamma(n, shape = shape, rate = 1)
  shares / sum(shares)
}
