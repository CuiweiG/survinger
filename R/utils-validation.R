#' Assert required columns exist in a data frame
#' @param df Data frame to check.
#' @param cols Character vector of required column names.
#' @param arg_name Name of the argument for error messages.
#' @return Invisible TRUE if all columns exist.
#' @keywords internal
.assert_columns_exist <- function(df, cols, arg_name = "data") {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      "Column{?s} {.val {missing_cols}} not found in {.arg {arg_name}}.",
      call = rlang::caller_env()
    )
  }
  invisible(TRUE)
}

#' Assert object inherits from surv_design
#' @param x Object to check.
#' @param arg_name Name for error messages.
#' @return Invisible TRUE.
#' @keywords internal
.assert_surv_design <- function(x, arg_name = "design") {
  if (!inherits(x, "surv_design")) {
    cli::cli_abort(
      "{.arg {arg_name}} must be a {.cls surv_design} object, not {.cls {class(x)[1]}}.",
      call = rlang::caller_env()
    )
  }
  invisible(TRUE)
}

#' Assert object inherits from a specific surv_* class
#' @param x Object to check.
#' @param cls Expected class name.
#' @param arg_name Name for error messages.
#' @return Invisible TRUE.
#' @keywords internal
.assert_surv_class <- function(x, cls, arg_name = deparse(substitute(x))) {
  if (!inherits(x, cls)) {
    cli::cli_abort(
      "{.arg {arg_name}} must be a {.cls {cls}} object, not {.cls {class(x)[1]}}.",
      call = rlang::caller_env()
    )
  }
  invisible(TRUE)
}

#' Validate and clamp sequencing rates to (0, 1]
#' @param rates Numeric vector of sequencing rates.
#' @return Numeric vector clamped to 0.001--1.
#' @keywords internal
.validate_seq_rates <- function(rates) {
  if (any(is.na(rates))) {
    cli::cli_abort("Sequencing rates contain NA values.", call = rlang::caller_env())
  }
  if (any(rates <= 0, na.rm = TRUE)) {
    cli::cli_warn("Sequencing rates <= 0 detected. Clamping to 0.001.")
    rates <- pmax(rates, 0.001)
  }
  if (any(rates > 1, na.rm = TRUE)) {
    cli::cli_warn("Sequencing rates > 1 detected. Clamping to 1.0.")
    rates <- pmin(rates, 1.0)
  }
  rates
}
