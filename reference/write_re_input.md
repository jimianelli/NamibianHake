# Write an input file for the RE model

Write an input file for the RE model

## Usage

``` r
write_re_input(
  file,
  years,
  biomass,
  CV,
  tag = NULL,
  sp = "species",
  survey = "survey",
  yr_start = NULL,
  yr_end = NULL
)
```

## Arguments

- file:

  A path to the output file, usually ending in .dat

- years:

  A vector of years for the observations

- biomass:

  A vector of biomasses for the observations

- CV:

  A vector of CV values for the observations

- tag:

  Optional character string to print to file for bookkeeping (e.g., "EBS
  shelf flathead sole")

- yr_start:

  Optional start year for predictions

- yr_end:

  Optional end year for predictions

## Value

Nothing, a data file is written to `file`

## See also

read_re_output
