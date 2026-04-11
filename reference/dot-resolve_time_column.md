# Resolve user-specified time column

Handles "epiweek", "month", "date", or a literal column name.

## Usage

``` r
.resolve_time_column(dat, time, date_col)
```

## Arguments

- dat:

  Data frame.

- time:

  User specification string.

- date_col:

  Name of the date column in dat.

## Value

List with data (possibly modified) and col (column name).
