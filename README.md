
README
================
John Kathena
Jim Ianelli
2024-07-19

<!-- badges: start -->
<!-- badges: end -->

The goal of NamibianHake is to setup a package for the Namibian hake assessment.

## Installation

You can install the development version of NamibianHake from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jimianelli/NamibianHake")
```

## Reproducing the 2026 update

Run the annual workflow from the repository root in this order:

```sh
Rscript scripts/01_validate_inputs.R
Rscript scripts/02_run_assessment.R
Rscript scripts/03_check_models.R
Rscript scripts/04_render_report.R
```

The stages validate and compare the annual inputs, compile and run the ADMB
models, classify convergence, write compact CSV summaries, and render the HTML
and PDF reports. Review these colleague-facing outputs:

- `reports/input_changes_2025_to_2026.html`
- `reports/model_diagnostics_2026.csv`
- `reports/key_results_2026.csv`
- `vignettes/02-Namibian_hake_model_2026.html`
- `vignettes/02-Namibian_hake_model_2026.pdf`

The dome-shaped selectivity and steepness 0.5 sensitivities currently fail the
declared diagnostic criteria. Their estimates are retained for troubleshooting
but excluded from inference.

## Example model runs

Below are examples which show how models can be run based on the directory location.

``` r
library(NamibianHake)
## basic example code
library(here)
library(tidyverse)
library(ggridges)
theme_set(ggthemes::theme_few())

bc <- run_nh("bc")
m1 <- run_nh("m1",runit=FALSE)
m2 <- run_nh("m2",runit=FALSE)
m3 <- run_nh("m3",runit=FALSE)
mods <-  rbind(
    read_csv(here("mods", "bc", "nh_out.csv")) |> mutate(Model = "Base Case"),
    read_csv(here("mods", "m1", "nh_out.csv")) |> mutate(Model = "Model 1"),
    read_csv(here("mods", "m2", "nh_out.csv")) |> mutate(Model = "Model 2"),
    read_csv(here("mods", "m3", "nh_out.csv")) |> mutate(Model = "Model 3")
    )
  

df <- rbind(
    data.frame(SSB = bc$SSB, R = bc$Pred_Rec, Model = "Base Case"),
    data.frame(SSB = m1$SSB, R = m1$Pred_Rec, Model = "Model 1"),
    data.frame(SSB = m2$SSB, R = m2$Pred_Rec, Model = "Model 2"),
    data.frame(SSB = m3$SSB, R = m3$Pred_Rec, Model = "Model 3") 
  ) 

p1 <- df %>% ggplot(aes(x = SSB, y = R, color = Model, shape=Model)) +
  geom_point() +
  geom_line();p1
PlotAgeFit(x=bc, title="Base case",type="fishery",fage=2,lage=7)
PlotAgeFit(x=bc, title="Base case",type="survey1",fage=1,lage=7)
M<- data.frame(Year = 1964:2023,Sel = bc$S[1:60,], Model = "h=0.7, base case")
plot_sel()

```
