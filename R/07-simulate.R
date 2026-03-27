# ============================================================
# Surveillance data simulator
# ============================================================

#' Simulate genomic surveillance data
#'
#' Generates synthetic surveillance datasets with realistic features:
#' multiple regions with unequal sequencing rates, multiple lineages
#' with time-varying prevalence, configurable reporting delays, and
#' multiple sample sources.
#'
#' @param n_regions Integer. Number of geographic regions. Default 5.
#' @param n_weeks Integer. Number of epiweeks. Default 26.
#' @param total_positive_per_week Integer. Mean total positive cases
#'   per week across all regions. Default 1000.
#' @param sequencing_rates Numeric vector of length `n_regions`.
#'   Per-region sequencing probability. If `NULL`, generated from
#'   a Beta distribution with realistic inequality. Default `NULL`.
#' @param lineage_dynamics Named list of functions, each taking a
#'   week number and returning a positive weight. If `NULL`, uses
#'   a default four-lineage scenario. Default `NULL`.
#' @param delay_params List with `mu` and `size` for negative binomial
#'   reporting delay. Default `list(mu = 10, size = 3)`.
#' @param sources Character vector of sample source types.
#'   Default `c("clinical", "wastewater", "sentinel")`.
#' @param source_weights Numeric vector (same length as `sources`).
#'   Default `c(0.7, 0.2, 0.1)`.
#' @param seed Integer or `NULL`. Random seed. Default `NULL`.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{sequences}{Tibble of individual sequence records.}
#'   \item{population}{Tibble with one row per region.}
#'   \item{truth}{Tibble of true lineage prevalence by region and week.}
#'   \item{parameters}{List of all input parameters.}
#' }
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 42)
#' head(sim$sequences)
#' sim$population
#'
#' @export
surv_simulate <- function(n_regions = 5L,
                          n_weeks = 26L,
                          total_positive_per_week = 1000L,
                          sequencing_rates = NULL,
                          lineage_dynamics = NULL,
                          delay_params = list(mu = 10, size = 3),
                          sources = c("clinical", "wastewater", "sentinel"),
                          source_weights = c(0.7, 0.2, 0.1),
                          seed = NULL) {

  checkmate::assert_count(n_regions, positive = TRUE)
  checkmate::assert_count(n_weeks, positive = TRUE)
  checkmate::assert_count(total_positive_per_week, positive = TRUE)
  checkmate::assert_list(delay_params)
  checkmate::assert_character(sources, min.len = 1L)
  checkmate::assert_numeric(source_weights, len = length(sources), lower = 0)

  if (!is.null(seed)) set.seed(seed)

  region_names <- paste0("Region_", LETTERS[seq_len(n_regions)])

  pop_shares <- stats::rgamma(n_regions, shape = 2, rate = 1)
  pop_shares <- pop_shares / sum(pop_shares)

  if (is.null(sequencing_rates)) {
    sequencing_rates <- stats::rbeta(n_regions, shape1 = 2, shape2 = 10)
    sequencing_rates <- pmax(sequencing_rates, 0.005)
    sequencing_rates <- pmin(sequencing_rates, 0.5)
  } else {
    checkmate::assert_numeric(sequencing_rates, len = n_regions,
                              lower = 0, upper = 1)
  }
  names(sequencing_rates) <- region_names

  if (is.null(lineage_dynamics)) {
    lineage_dynamics <- list(
      "BA.5" = function(w) pmax(0.01, 0.6 * exp(-0.05 * w)),
      "XBB.1.5" = function(w) pmax(0.01, 0.3 * exp(-0.03 * (w - 5))),
      "BA.2.86" = function(w) {
        if (w < 8) 0.01 else 0.01 + 0.5 * stats::plogis((w - 15) / 3)
      },
      "Other" = function(w) 0.10
    )
  }
  lineage_names <- names(lineage_dynamics)
  n_lineages <- length(lineage_names)

  all_records <- vector("list", n_regions * n_weeks)
  truth_records <- vector("list", n_regions * n_weeks)
  idx <- 0L

  for (w in seq_len(n_weeks)) {
    week_start <- as.Date("2024-01-01") + (w - 1L) * 7L
    ew_label <- paste0(format(week_start, "%G"), "-W", format(week_start, "%V"))

    global_prev <- vapply(lineage_dynamics, function(f) f(w), numeric(1))
    global_prev <- global_prev / sum(global_prev)

    for (r in seq_len(n_regions)) {
      idx <- idx + 1L
      region <- region_names[r]
      n_positive <- max(1L, stats::rpois(1L, total_positive_per_week * pop_shares[r]))
      n_sequenced <- stats::rbinom(1L, n_positive, sequencing_rates[r])

      region_prev <- global_prev + stats::rnorm(n_lineages, 0, 0.02)
      region_prev <- pmax(region_prev, 0.001)
      region_prev <- region_prev / sum(region_prev)

      truth_records[[idx]] <- tibble::tibble(
        epiweek = ew_label, region = region,
        lineage = lineage_names, true_prevalence = region_prev,
        n_positive = n_positive, n_sequenced = n_sequenced
      )

      if (n_sequenced == 0L) next

      coll_dates <- week_start + sample(0:6, n_sequenced, replace = TRUE)
      delays <- stats::rnbinom(n_sequenced, mu = delay_params$mu, size = delay_params$size)
      all_records[[idx]] <- tibble::tibble(
        sequence_id = paste0("seq_", w, "_", r, "_", seq_len(n_sequenced)),
        region = region,
        source_type = sample(sources, n_sequenced, replace = TRUE,
                             prob = source_weights / sum(source_weights)),
        lineage = sample(lineage_names, n_sequenced, replace = TRUE,
                         prob = region_prev),
        collection_date = coll_dates,
        report_date = coll_dates + delays,
        epiweek = ew_label
      )
    }
  }

  sequences <- dplyr::bind_rows(all_records)
  truth <- dplyr::bind_rows(truth_records)

  pop_summary <- truth |>
    dplyr::group_by(.data$region) |>
    dplyr::summarise(
      n_positive_total = sum(.data$n_positive),
      n_sequenced_total = sum(.data$n_sequenced),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      n_positive = as.integer(.data$n_positive_total / n_lineages),
      n_sequenced = as.integer(.data$n_sequenced_total / n_lineages),
      seq_rate = .data$n_sequenced / pmax(.data$n_positive, 1L),
      pop_total = as.integer(.data$n_positive / mean(sequencing_rates) * 10)
    ) |>
    dplyr::select("region", "n_positive", "n_sequenced", "seq_rate", "pop_total")

  list(
    sequences = sequences,
    population = pop_summary,
    truth = truth,
    parameters = list(
      n_regions = n_regions, n_weeks = n_weeks,
      total_positive_per_week = total_positive_per_week,
      sequencing_rates = sequencing_rates, delay_params = delay_params,
      sources = sources, source_weights = source_weights,
      seed = seed, lineage_names = lineage_names
    )
  )
}
