# Required sequences for target detection probability

Required sequences for target detection probability

## Usage

``` r
surv_required_sequences(
  true_prevalence,
  target_probability = 0.95,
  n_periods = 1L,
  detection_threshold = 1L
)
```

## Arguments

- true_prevalence:

  Numeric.

- target_probability:

  Numeric. Default 0.95.

- n_periods:

  Integer. Default 1.

- detection_threshold:

  Integer. Default 1.

## Value

Integer.

## Examples

``` r
surv_required_sequences(0.01)
#> [1] 299
surv_required_sequences(0.05, target_probability = 0.99)
#> [1] 90
```
