#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
repository <- dirname(dirname(normalizePath(script_file)))
setwd(repository)

source(file.path("R", "assessment-validation.R"))
source(file.path("R", "assessment-sensitivities.R"))

executable <- file.path("source", if (.Platform$OS.type == "windows") {
  "nh.exe"
} else {
  "nh"
})
if (!file.exists(executable) && .Platform$OS.type != "windows") {
  status <- system2("make", c("--directory=source"))
  if (status != 0L) {
    stop("Failed to compile the ADMB executable.")
  }
}
if (!file.exists(executable)) {
  stop("ADMB executable not found: ", executable)
}
executable <- normalizePath(executable)

sensitivity_controls <- create_annual_sensitivity_controls(
  file.path("mods", "m26", "nh.dat")
)
controls <- c(
  m25 = file.path("mods", "m25", "nh.dat"),
  m26 = file.path("mods", "m26", "nh.dat"),
  sensitivity_controls
)
restart_models <- c("m26_dome_selectivity", "m26_steepness_0_5")

statuses <- lapply(names(controls), function(model) {
  message("Running ", model)
  tryCatch(
    {
      first_pass_error <- tryCatch(
        {
          run_assessment_model(
            controls[[model]],
            executable,
            arguments = c(
              "-nox", "-iprint", "250", "-maxfn", "5000", "-hbf", "1"
            )
          )
          NULL
        },
        error = identity
      )

      if (model %in% restart_models) {
        run_directory <- dirname(controls[[model]])
        file.copy(
          file.path(run_directory, "run.stdout"),
          file.path(run_directory, "first-pass.stdout"),
          overwrite = TRUE
        )
        file.copy(
          file.path(run_directory, "run.stderr"),
          file.path(run_directory, "first-pass.stderr"),
          overwrite = TRUE
        )
        run_assessment_model(
          controls[[model]],
          executable,
          arguments = c("-binp", "nh.bar", "-phase", "22")
        )
        file.copy(
          file.path(run_directory, "run.stdout"),
          file.path(run_directory, "second-pass.stdout"),
          overwrite = TRUE
        )
        file.copy(
          file.path(run_directory, "run.stderr"),
          file.path(run_directory, "second-pass.stderr"),
          overwrite = TRUE
        )
      } else if (inherits(first_pass_error, "error")) {
        stop(first_pass_error)
      }

      data.frame(
        Model = model,
        Exit_status = 0L,
        Message = if (model %in% restart_models) {
          "Completed with nh.bar restart"
        } else {
          "Completed"
        }
      )
    },
    error = function(e) {
      data.frame(
        Model = model,
        Exit_status = 1L,
        Message = conditionMessage(e)
      )
    }
  )
}) |>
  (\(x) do.call(rbind, x))()

dir.create("reports", showWarnings = FALSE)
write.csv(
  statuses,
  file.path("reports", "run_status_2026.csv"),
  row.names = FALSE
)
print(statuses, row.names = FALSE)

required <- statuses$Model %in% c("m25", "m26")
if (any(statuses$Exit_status[required] != 0L)) {
  stop("A required bridge model failed.")
}
