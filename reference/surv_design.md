# Create a genomic surveillance design object

Constructs a survey design object tailored for pathogen genomic
surveillance, encoding stratification structure, sequencing rates, and
population information needed for design-weighted inference.

## Usage

``` r
# S3 method for class 'surv_design'
print(x, ...)

# S3 method for class 'surv_design'
summary(object, ...)

# S3 method for class 'summary.surv_design'
print(x, ...)

surv_design(
  data,
  strata,
  sequencing_rate,
  population,
  date_collected = "collection_date",
  date_reported = "report_date",
  lineage = "lineage",
  source_type = NULL,
  source_config = NULL
)
```

## Arguments

- x:

  Object to print or summarize.

- ...:

  Additional arguments (unused).

- object:

  A `surv_design` object to summarize.

- data:

  Data frame of individual sequence records.

- strata:

  One-sided formula specifying stratification variables (e.g.,
  `~ region` or `~ region + source_type`).

- sequencing_rate:

  Either a one-sided formula (`~ n_sequenced / n_positive`), a named
  numeric vector, or a data frame with strata columns and `seq_rate`.

- population:

  Data frame with one row per stratum, containing stratification
  variables and population-level denominators.

- date_collected:

  Column name for collection date. Default `"collection_date"`.

- date_reported:

  Column name for report date. Default `"report_date"`. Set `NULL` if
  unavailable.

- lineage:

  Column name for lineage. Default `"lineage"`.

- source_type:

  Column name for sample source. Default `NULL`.

- source_config:

  Optional tibble of per-source characteristics.

## Value

Invisibly returns the input object.

A summary list of class `summary.surv_design`.

Invisibly returns the input object.

An object of class `surv_design`.

## See also

[`surv_simulate()`](https://cuiweig.github.io/survinger/reference/surv_simulate.md),
[`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md),
[`surv_optimize_allocation()`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 42)
design <- surv_design(
  data = sim$sequences,
  strata = ~ region,
  sequencing_rate = sim$population[c("region", "seq_rate")],
  population = sim$population
)
print(design)
#> ── Genomic Surveillance Design ─────────────────────────────────────────────────
#> Observations: 775
#> Strata: 3 (region)
#> Date range: 2024-01-01 to 2024-02-25
#> Lineages: 4
#> Weight range: [3.716, 18.053]
```
