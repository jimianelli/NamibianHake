# Test the RE on a single species. This compares a current ADMB univariate run to a previous run. It tests *continuity* not correctness

Test the RE on a single species. This compares a current ADMB univariate
run to a previous run. It tests *continuity* not correctness

## Usage

``` r
test_re_species(d, test_tmb = TRUE)
```

## Arguments

- d:

  The .dat file name, e.g., bigskate_t.dat, found in tests/data folder

## Value

A data.frame with ADMB and TMB estimates for comparing

## Details

It will throw an error if the results have changed, otherwise results
are consistent.
