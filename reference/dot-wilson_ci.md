# Wilson score interval for a proportion

Provides valid coverage even at p-hat = 0 or 1, unlike Wald. Reference:
Wilson (1927) JASA; Agresti & Coull (1998) TAS.

## Usage

``` r
.wilson_ci(p, n, z, side = c("lower", "upper"))
```

## Arguments

- p:

  Estimated proportion.

- n:

  Effective sample size.

- z:

  Critical value (e.g., qnorm(0.975) for 95% CI).

- side:

  Character: "lower" or "upper".

## Value

Scalar bound, clamped to 0–1.
