#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
repository <- dirname(dirname(normalizePath(script_file)))
setwd(repository)

source(file.path("R", "read-admb2.R"))
source(file.path("R", "assessment-validation.R"))

gradient_tolerance <- 0.001
models <- data.frame(
  Directory = c(
    "m25",
    "m26",
    "m26_survey_q_1",
    "m26_survey_q_0_5",
    "m26_estimated_natural_mortality",
    "m26_dome_selectivity",
    "m26_steepness_0_5",
    "m26_steepness_0_9"
  ),
  Model = c(
    "2025 baseline",
    "2026 update",
    "Survey q = 1.0",
    "Survey q = 0.5",
    "Estimated natural mortality",
    "Dome-shaped selectivity",
    "Steepness = 0.5",
    "Steepness = 0.9"
  ),
  Terminal_year = c(2025L, rep(2026L, 7L))
)

diagnostics <- lapply(seq_len(nrow(models)), function(index) {
  model <- models[index, ]
  directory <- file.path("mods", model$Directory)
  tryCatch(
    {
      convergence <- read_admb_convergence(directory)
      output <- read.csv(file.path(directory, "nh_out.csv"))
      output_years <- output$Year[is.finite(output$Year)]
      output_year <- if (length(output_years)) {
        max(output_years)
      } else {
        NA_integer_
      }
      converged <- is.finite(convergence$objective) &&
        is.finite(convergence$max_gradient) &&
        convergence$max_gradient <= gradient_tolerance &&
        convergence$positive_definite_hessian &&
        output_year == model$Terminal_year
      data.frame(
        Directory = model$Directory,
        Model = model$Model,
        Terminal_year = output_year,
        Objective = convergence$objective,
        Maximum_gradient = convergence$max_gradient,
        Status = if (converged) "Converged" else "Failed diagnostic criteria"
      )
    },
    error = function(e) {
      data.frame(
        Directory = model$Directory,
        Model = model$Model,
        Terminal_year = NA_integer_,
        Objective = NA_real_,
        Maximum_gradient = NA_real_,
        Status = paste("Failed:", conditionMessage(e))
      )
    }
  )
}) |>
  (\(x) do.call(rbind, x))()

key_results <- lapply(seq_len(nrow(diagnostics)), function(index) {
  if (diagnostics$Status[[index]] != "Converged") {
    return(NULL)
  }
  output <- read.csv(
    file.path("mods", diagnostics$Directory[[index]], "nh_out.csv")
  )
  output[
    output$Year == diagnostics$Terminal_year[[index]] &
      output$Variable %in%
        c("SSB", "Depletion", "B_Bmsy", "Catch", "RY", "Catch_RY"),
    c("Variable", "Year", "value", "ymin", "ymax")
  ] |>
    transform(
      Directory = diagnostics$Directory[[index]],
      Model = diagnostics$Model[[index]]
    )
}) |>
  (\(x) do.call(rbind, x))()

dir.create("reports", showWarnings = FALSE)
write.csv(
  diagnostics,
  file.path("reports", "model_diagnostics_2026.csv"),
  row.names = FALSE
)
write.csv(
  key_results,
  file.path("reports", "key_results_2026.csv"),
  row.names = FALSE
)
print(diagnostics, row.names = FALSE)

base_status <- diagnostics$Status[diagnostics$Directory == "m26"]
if (!identical(base_status, "Converged")) {
  stop("The 2026 base model did not pass diagnostic criteria.")
}
