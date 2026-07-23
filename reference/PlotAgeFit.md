# Plot Age Fit

This function creates a plot comparing observed and predicted age
compositions for a specified type.

## Usage

``` r
PlotAgeFit(x = bc, title = NULL, type = "fishery", fage = 2, lage = 7)
```

## Arguments

- x:

  A list containing the observed and predicted age compositions. Default
  is `bc`.

- title:

  A character string for the plot title. Default is NULL.

- type:

  A character string specifying the type of composition to plot. Default
  is "fishery".

- fage:

  An integer specifying the first age to include in the plot. Default is
  2.

- lage:

  An integer specifying the last age to include in the plot. Default is
  7.

## Value

A ggplot object comparing observed and predicted age compositions.

## Examples

``` r
if (FALSE) { # \dontrun{
PlotAgeFit(x = bc, title = "Age Composition", type = "survey", fage = 1, lage = 8)
} # }
```
