# Read output file from the random effects (RE) model into a data.frame

Read output file from the random effects (RE) model into a data.frame

## Usage

``` r
read_re_output(file = "rwout.rep", skip_data = FALSE, use_names = FALSE)
```

## Arguments

- file:

  A character string of the file to be read in, defaulting to
  'rwout.rep' in the current working directory.

- skip_data:

  Whether to skip returning columns representing the data inputs. If
  TRUE these are included and will have NA values for any missing years

- use_names:

  Whether to use the names written by the RE model or to use more
  meaningful ones (default)

## Value

A dataframe with columns named based on the tags in the output file
