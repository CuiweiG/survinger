# Publication-quality ggplot2 theme

A clean, high-contrast theme suitable for journal publication. Follows
Nature/Lancet figure guidelines.

## Usage

``` r
theme_survinger(base_size = 11, base_family = "")
```

## Arguments

- base_size:

  Base font size. Default 11.

- base_family:

  Font family. Default `""`.

## Value

A ggplot2 theme object.

## See also

[`surv_compare_estimates()`](https://cuiweig.github.io/survinger/reference/surv_compare_estimates.md)

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_survinger()

```
