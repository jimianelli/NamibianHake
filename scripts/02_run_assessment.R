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

statuses <- lapply(names(controls), function(model) {
  message("Running ", model)
  tryCatch(
    {
      run_assessment_model(
        controls[[model]],
        executable,
        arguments = c(
          "-nox", "-iprint", "250", "-maxfn", "5000", "-hbf", "1"
        )
      )
      data.frame(Model = model, Exit_status = 0L, Message = "Completed")
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
