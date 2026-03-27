#' Example SARS-CoV-2 genomic surveillance data
#'
#' Simulated genomic surveillance dataset with 5 regions, 26 weeks,
#' highly unequal sequencing rates (15% to 0.5%), three sample sources,
#' and negative binomial reporting delays. Contains known ground truth
#' for benchmarking.
#'
#' @format A named list with four elements:
#' \describe{
#'   \item{sequences}{Tibble of sequence records: sequence_id, region,
#'     source_type, lineage, collection_date, report_date, epiweek.}
#'   \item{population}{Tibble with one row per region: region,
#'     n_positive, n_sequenced, seq_rate, pop_total.}
#'   \item{truth}{Tibble of true lineage prevalence by region and week.}
#'   \item{parameters}{List of simulation parameters.}
#' }
#'
#' @source Simulated using `surv_simulate(seed = 20240101)`.
#'
#' @examples
#' data(sarscov2_surveillance)
#' head(sarscov2_surveillance$sequences)
#' sarscov2_surveillance$population
"sarscov2_surveillance"
