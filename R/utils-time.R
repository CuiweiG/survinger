#' Convert Date vector to ISO epiweek strings (YYYY-Wnn)
#' @param dates Date vector.
#' @return Character vector of epiweek labels.
#' @keywords internal
.date_to_epiweek <- function(dates) {
  dates <- as.Date(dates)
  yr <- format(dates, "%G")
  wk <- format(dates, "%V")
  paste0(yr, "-W", wk)
}

#' Resolve user-specified time column
#'
#' Handles "epiweek", "month", "date", or a literal column name.
#'
#' @param dat Data frame.
#' @param time User specification string.
#' @param date_col Name of the date column in dat.
#' @return List with data (possibly modified) and col (column name).
#' @keywords internal
.resolve_time_column <- function(dat, time, date_col) {
  if (identical(time, "epiweek")) {
    if (!"epiweek" %in% names(dat)) {
      dat[["epiweek"]] <- .date_to_epiweek(dat[[date_col]])
    }
    return(list(data = dat, col = "epiweek"))
  }
  if (identical(time, "month")) {
    if (!"month" %in% names(dat)) {
      dat[["month"]] <- format(as.Date(dat[[date_col]]), "%Y-%m")
    }
    return(list(data = dat, col = "month"))
  }
  if (identical(time, "date")) {
    return(list(data = dat, col = date_col))
  }
  if (time %in% names(dat)) {
    return(list(data = dat, col = time))
  }
  cli::cli_abort(
    c("Time specification {.val {time}} not recognized.",
      "i" = "Use {.val epiweek}, {.val month}, {.val date}, or a column name."),
    call = rlang::caller_env()
  )
}
