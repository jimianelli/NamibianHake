create_annual_sensitivity_controls <- function(
  base_control,
  output_root = dirname(dirname(base_control)),
  assessment_year = 2026L
) {
  base <- readLines(base_control, warn = FALSE)
  stopifnot(
    basename(base_control) == "nh.dat",
    as.integer(base[[7L]]) == assessment_year,
    trimws(base[[144L]]) == "54321"
  )

  variants <- list(
    m26_survey_q_1 = list(
      label = "#M26_q1",
      descriptor = "minus1,h7,asymp,sigmafixed,tvs,q1",
      changes = c(`88` = "1.0", `89` = "1.0")
    ),
    m26_survey_q_0_5 = list(
      label = "#M26_q05",
      descriptor = "minus1,h7,asymp,sigmafixed,tvs,q0.5",
      changes = c(`88` = "0.5", `89` = "0.5")
    ),
    m26_estimated_natural_mortality = list(
      label = "#M26_M_estimated",
      descriptor = "minus1,h7,asymp,sigmafixed,tvs,Mestimated",
      changes = c(`11` = "6")
    ),
    m26_dome_selectivity = list(
      label = "#M26_dome_selectivity",
      descriptor = "minus1,h7,dome,sigmafixed,tvs",
      changes = c(`41` = "4", `44` = "4")
    ),
    m26_steepness_0_5 = list(
      label = "#M26_h05",
      descriptor = "minus1,h5,asymp,sigmafixed,tvs",
      changes = c(`16` = "0.5")
    ),
    m26_steepness_0_9 = list(
      label = "#M26_h09",
      descriptor = "minus1,h9,asymp,sigmafixed,tvs",
      changes = c(`16` = "0.9")
    )
  )

  paths <- character(length(variants))
  names(paths) <- names(variants)
  for (name in names(variants)) {
    variant <- variants[[name]]
    lines <- base
    lines[[2L]] <- variant$label
    lines[[3L]] <- variant$descriptor
    lines[[5L]] <- "20"
    lines[as.integer(names(variant$changes))] <- unname(variant$changes)

    directory <- file.path(output_root, name)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    paths[[name]] <- file.path(directory, "nh.dat")
    writeLines(lines, paths[[name]], useBytes = TRUE)
  }
  paths
}

run_assessment_model <- function(control_file, executable,
                                 arguments = c(
                                   "-nox", "-iprint", "150",
                                   "-hbf", "1"
                                 )) {
  run_directory <- dirname(control_file)
  executable <- normalizePath(executable, mustWork = TRUE)
  old_directory <- getwd()
  on.exit(setwd(old_directory), add = TRUE)
  setwd(run_directory)
  status <- system2(
    executable,
    args = arguments,
    stdout = "run.stdout",
    stderr = "run.stderr"
  )
  if (status != 0L) {
    stop("ADMB model failed with exit status ", status, ": ", run_directory)
  }
  invisible(status)
}
