# Simulate genomic surveillance data

Generates synthetic surveillance datasets with realistic features:
multiple regions with unequal sequencing rates, multiple lineages with
time-varying prevalence, configurable reporting delays, and multiple
sample sources.

## Usage

``` r
surv_simulate(
  n_regions = 5L,
  n_weeks = 26L,
  total_positive_per_week = 1000L,
  sequencing_rates = NULL,
  lineage_dynamics = NULL,
  delay_params = list(mu = 10, size = 3),
  sources = c("clinical", "wastewater", "sentinel"),
  source_weights = c(0.7, 0.2, 0.1),
  seed = NULL
)
```

## Arguments

- n_regions:

  Integer. Number of geographic regions. Default 5.

- n_weeks:

  Integer. Number of epiweeks. Default 26.

- total_positive_per_week:

  Integer. Mean total positive cases per week across all regions.
  Default 1000.

- sequencing_rates:

  Numeric vector of length `n_regions`. Per-region sequencing
  probability. If `NULL`, generated from a Beta distribution with
  realistic inequality. Default `NULL`.

- lineage_dynamics:

  Named list of functions, each taking a week number and returning a
  positive weight. If `NULL`, uses a default four-lineage scenario.
  Default `NULL`.

- delay_params:

  List with `mu` and `size` for negative binomial reporting delay.
  Default `list(mu = 10, size = 3)`.

- sources:

  Character vector of sample source types. Default
  `c("clinical", "wastewater", "sentinel")`.

- source_weights:

  Numeric vector (same length as `sources`). Default `c(0.7, 0.2, 0.1)`.

- seed:

  Integer or `NULL`. Random seed. Default `NULL`.

## Value

A named list with elements:

- sequences:

  Tibble of individual sequence records.

- population:

  Tibble with one row per region.

- truth:

  Tibble of true lineage prevalence by region and week.

- parameters:

  List of all input parameters.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 42)
head(sim$sequences)
#> # A tibble: 6 × 7
#>   sequence_id region   source_type lineage collection_date report_date epiweek 
#>   <chr>       <chr>    <chr>       <chr>   <date>          <date>      <chr>   
#> 1 seq_1_1_1   Region_A wastewater  BA.5    2024-01-07      2024-01-10  2024-W01
#> 2 seq_1_1_2   Region_A clinical    BA.5    2024-01-05      2024-01-07  2024-W01
#> 3 seq_1_1_3   Region_A clinical    BA.5    2024-01-05      2024-01-12  2024-W01
#> 4 seq_1_1_4   Region_A clinical    BA.5    2024-01-04      2024-01-09  2024-W01
#> 5 seq_1_1_5   Region_A clinical    XBB.1.5 2024-01-02      2024-01-07  2024-W01
#> 6 seq_1_1_6   Region_A wastewater  XBB.1.5 2024-01-04      2024-01-04  2024-W01
sim$population
#> # A tibble: 3 × 5
#>   region   n_positive n_sequenced seq_rate pop_total
#>   <chr>         <int>       <int>    <dbl>     <int>
#> 1 Region_A       4784         265   0.0554    350844
#> 2 Region_B       1163         313   0.269      85291
#> 3 Region_C       2112         197   0.0933    154887
```
