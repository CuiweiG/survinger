# Example SARS-CoV-2 genomic surveillance data

Simulated genomic surveillance dataset with 5 regions, 26 weeks, highly
unequal sequencing rates (15% to 0.5%), three sample sources, and
negative binomial reporting delays. Contains known ground truth for
benchmarking.

## Usage

``` r
sarscov2_surveillance
```

## Format

A named list with four elements:

- sequences:

  Tibble of sequence records: sequence_id, region, source_type, lineage,
  collection_date, report_date, epiweek.

- population:

  Tibble with one row per region: region, n_positive, n_sequenced,
  seq_rate, pop_total.

- truth:

  Tibble of true lineage prevalence by region and week.

- parameters:

  List of simulation parameters.

## Source

Simulated using `surv_simulate(seed = 20240101)`.

## Examples

``` r
data(sarscov2_surveillance)
head(sarscov2_surveillance$sequences)
#> # A tibble: 6 × 7
#>   sequence_id region   source_type lineage collection_date report_date epiweek 
#>   <chr>       <chr>    <chr>       <chr>   <date>          <date>      <chr>   
#> 1 seq_1_1_1   Region_A wastewater  XBB.1.5 2024-01-02      2024-01-13  2024-W01
#> 2 seq_1_1_2   Region_A wastewater  BA.5    2024-01-07      2024-01-11  2024-W01
#> 3 seq_1_1_3   Region_A wastewater  BA.5    2024-01-01      2024-01-07  2024-W01
#> 4 seq_1_1_4   Region_A clinical    Other   2024-01-01      2024-01-12  2024-W01
#> 5 seq_1_1_5   Region_A wastewater  BA.5    2024-01-04      2024-01-12  2024-W01
#> 6 seq_1_1_6   Region_A clinical    Other   2024-01-04      2024-01-15  2024-W01
sarscov2_surveillance$population
#> # A tibble: 5 × 5
#>   region   n_positive n_sequenced seq_rate pop_total
#>   <chr>         <int>       <int>    <dbl>     <int>
#> 1 Region_A       6784         983  0.145     1233454
#> 2 Region_B       2117         183  0.0864     384909
#> 3 Region_C       2385          72  0.0302     433636
#> 4 Region_D       8712          84  0.00964   1584000
#> 5 Region_E       6095          27  0.00443   1108181
```
