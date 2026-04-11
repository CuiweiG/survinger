# Compute power curve for detection across prevalence range

Generates a detection probability curve that can be directly plotted or
included in publications. Answers: "At what prevalence does our
surveillance achieve X% detection?"

## Usage

``` r
surv_power_curve(
  design,
  prevalence_range = seq(0.001, 0.05, by = 0.001),
  delay_fit = NULL,
  thresholds = c(0.5, 0.8, 0.95)
)

# S3 method for class 'surv_power_curve'
plot(x, ...)
```

## Arguments

- design:

  A `surv_design` object.

- prevalence_range:

  Numeric vector of prevalences to evaluate. Default
  `seq(0.001, 0.05, by = 0.001)`.

- delay_fit:

  Optional `surv_delay_fit`.

- thresholds:

  Numeric vector of detection thresholds to mark. Default
  `c(0.5, 0.8, 0.95)`.

- x:

  A `surv_power_curve` object.

- ...:

  Additional arguments (unused).

## Value

A list with:

- curve:

  Tibble with prevalence and detection columns.

- thresholds:

  Tibble with threshold, prevalence_needed columns.

A ggplot2 object.

## See also

[`surv_detection_probability()`](https://cuiweig.github.io/survinger/reference/surv_detection_probability.md),
[`surv_required_sequences()`](https://cuiweig.github.io/survinger/reference/surv_required_sequences.md)

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
pc <- surv_power_curve(d)
pc$thresholds
#> # A tibble: 3 × 2
#>   threshold prevalence_needed
#>       <dbl>             <dbl>
#> 1      0.5              0.01 
#> 2      0.8              0.021
#> 3      0.95             0.039
```
