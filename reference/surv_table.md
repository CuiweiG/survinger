# Format prevalence results for knitr tables

Produces a publication-ready summary table from any survinger result.
Automatically called when objects are printed inside RMarkdown chunks if
`knitr` is loaded.

## Usage

``` r
surv_table(x, digits = 3, percent = TRUE, ...)
```

## Arguments

- x:

  A survinger result object.

- digits:

  Integer. Decimal places for prevalence. Default 3.

- percent:

  Logical. Display as percentages? Default `TRUE`.

- ...:

  Additional arguments (unused).

## Value

A tibble formatted for display.

## See also

[`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md),
[`tidy.surv_prevalence()`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
prev <- surv_lineage_prevalence(d, "BA.2.86")
surv_table(prev)
#> # A tibble: 10 × 7
#>    time     lineage prevalence ci               se n_obs effective_n
#>    <chr>    <chr>        <dbl> <chr>         <dbl> <int>       <dbl>
#>  1 2024-W01 BA.2.86       4.27 [1.28, 13.26]  3.18    74        54.7
#>  2 2024-W02 BA.2.86       3.66 [0.91, 13.55]  3.66    58        46.0
#>  3 2024-W03 BA.2.86       0    [0, 6.27]      0       74        57.5
#>  4 2024-W04 BA.2.86       0    [0, 6.07]      0       78        59.4
#>  5 2024-W05 BA.2.86       0    [0, 6.51]      0       78        55.1
#>  6 2024-W06 BA.2.86       2.83 [0.71, 10.68]  2.44    82        59.5
#>  7 2024-W07 BA.2.86       2.18 [0.48, 9.37]   1.52    82        63.1
#>  8 2024-W08 BA.2.86       6.05 [2.26, 15.2]   2.99    78        59.2
#>  9 2024-W09 BA.2.86       9.07 [3.9, 19.72]   4.20    69        53.7
#> 10 2024-W10 BA.2.86      12.3  [6.38, 22.3]   4.03    86        66.2
```
