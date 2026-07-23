# Namibian hake stock assessments

<!-- badges: start -->
<!-- badges: end -->

This repository contains the annual Namibian hake stock assessment, its data
and control files, model diagnostics, and the code needed to reproduce the
results.

## Assessment documents

| Assessment | What it represents | Read online | Files |
|:--|:--|:--|:--|
| **2026 update** | Current assessment using data through 2026 | [2026 report](articles/02-Namibian_hake_model_2026.html) | [PDF](articles/02-Namibian_hake_model_2026.pdf) · [source](https://github.com/jimianelli/NamibianHake/blob/main/reports/02-Namibian_hake_model_2026.qmd) |
| **2025 assessment** | Previous-year baseline used for the bridge comparison | [2025 report](articles/01-Namibian_hake_model_2025.html) | [archived source](https://github.com/jimianelli/NamibianHake/blob/main/reports/archive/01-Namibian_hake_model_2025.qmd) |
| **2024 assessment** | Historical assessment | — | [PDF](articles/00-Namibian_hake_model_2024.pdf) · [archived source](https://github.com/jimianelli/NamibianHake/blob/main/reports/archive/00-Namibian_hake_model_2024.qmd) |

The [2025-to-2026 input-difference
report](https://github.com/jimianelli/NamibianHake/blob/main/reports/input_changes_2025_to_2026.md)
identifies every intended input change. Model diagnostics and key results are
available in the [`reports`
directory](https://github.com/jimianelli/NamibianHake/tree/main/reports).

## Reproduce the 2026 assessment

From the repository root, run the four numbered workflow scripts in order:

``` sh
Rscript scripts/01_validate_inputs.R
Rscript scripts/02_run_assessment.R
Rscript scripts/03_check_models.R
Rscript scripts/04_render_report.R
```

These stages validate and compare the annual inputs, compile and run the ADMB
models, classify convergence, write compact CSV summaries, and render the HTML
and PDF reports. The main machine-readable outputs are:

- `reports/input_changes_2025_to_2026.html`
- `reports/model_diagnostics_2026.csv`
- `reports/key_results_2026.csv`
- `vignettes/02-Namibian_hake_model_2026.html`
- `vignettes/02-Namibian_hake_model_2026.pdf`

The dome-shaped selectivity and steepness 0.5 sensitivities currently fail the
declared diagnostic criteria. Their estimates are retained for troubleshooting
but excluded from inference.

For the exact 2025 reference state, use the immutable
[`baseline-2025`](https://github.com/jimianelli/NamibianHake/releases/tag/baseline-2025)
tag.
