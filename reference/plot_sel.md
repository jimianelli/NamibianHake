# Plot Selectivity

This function creates a plot of selectivity across years and ages.

## Usage

``` r
plot_sel(
  Year = M$Year,
  sel = M[, 2:10],
  styr = 1964,
  fage = NULL,
  lage = NULL,
  alpha = 0.2,
  scale = 3.8,
  fill = "purple"
)
```

## Arguments

- Year:

  A vector of years. Default is `M$Year`.

- sel:

  A matrix of selectivity values. Default is `M[,2:10]`.

- styr:

  An integer specifying the start year for the plot. Default is 1964.

- fage:

  An integer specifying the first age to include in the plot. Default is
  NULL.

- lage:

  An integer specifying the last age to include in the plot. Default is
  NULL.

- alpha:

  A numeric value specifying the transparency level of the fill. Default
  is 0.2.

- scale:

  A numeric value specifying the scale for the density ridges. Default
  is 3.8.

- fill:

  A character string specifying the fill color for the density ridges.
  Default is "purple".

## Value

A ggplot object showing selectivity across years and ages.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_sel(Year = M$Year, sel = M[, 2:10], styr = 1970, fage = 1, lage = 8)
} # }
```
