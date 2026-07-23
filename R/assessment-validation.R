read_nh_control <- function(file) {
  tokens <- scan(
    file,
    what = character(),
    comment.char = "#",
    quiet = TRUE
  )
  position <- 4L
  take <- function(n = 1L) {
    values <- as.numeric(tokens[position:(position + n - 1L)])
    position <<- position + n
    values
  }

  control <- list(
    model_number = as.integer(tokens[[1L]]),
    model_name = tokens[[2L]],
    data_file = tokens[[3L]]
  )
  control$n_projection <- as.integer(take())
  control$first_year <- as.integer(take())
  control$last_year <- as.integer(take())
  control$plus_group <- as.integer(take())
  control$minimum_f_age <- as.integer(take())
  control$estimation_phases <- take(6L)
  control$steepness <- take()
  control$steepness_cv <- take()
  control$ksp_fraction <- take()

  control$n_selectivity_periods <- as.integer(take())
  control$selectivity_years <- matrix(
    take(2L * control$n_selectivity_periods),
    ncol = 2L,
    byrow = TRUE
  )
  control$selectivity_settings <- take(7L)
  control$sigma_recruitment <- take()
  control$recruitment_years <- as.integer(take(2L))

  control$caa_minus <- as.integer(take())
  control$caa_plus <- as.integer(take())
  control$n_cpue_series <- as.integer(take())
  take(4L * control$n_cpue_series)

  control$n_survey_series <- as.integer(take())
  take(6L * control$n_survey_series)
  take()
  take(2L * control$n_survey_series)
  take(4L)

  control$n_ksp_periods <- as.integer(take())
  control$ksp_years <- matrix(
    take(2L * control$n_ksp_periods),
    ncol = 2L,
    byrow = TRUE
  )
  control$economic_parameters <- take(length(tokens) - position)
  control$eof <- take()
  control
}

read_nh_data <- function(file, first_year, last_year, plus_group = 8L,
                         n_cpue_series = 6L, n_survey_series = 2L) {
  values <- scan(file, quiet = TRUE, comment.char = "#")
  position <- 1L
  take <- function(n) {
    result <- values[position:(position + n - 1L)]
    position <<- position + n
    result
  }
  take_matrix <- function(nrow, ncol) {
    matrix(take(nrow * ncol), nrow = nrow, ncol = ncol, byrow = TRUE)
  }

  n_years <- last_year - first_year + 1L
  weight_columns <- plus_group + 2L
  composition_columns <- plus_group + 3L

  result <- list(
    years = seq.int(first_year, last_year),
    weight_start = take_matrix(n_years, weight_columns),
    weight_mid = take_matrix(n_years, weight_columns),
    catch = take(n_years),
    fishery_caa = take_matrix(n_years, composition_columns),
    cpue = take_matrix(n_years, n_cpue_series + 1L),
    survey = take_matrix(n_years, 2L * n_survey_series + 1L)
  )
  result$survey_caa <- lapply(
    seq_len(n_survey_series),
    \(x) take_matrix(n_years, composition_columns)
  )
  result$seal_index <- take(n_years)
  result$seal_cv <- take(n_years)
  result$maturity <- take(plus_group + 1L)
  result$eof <- take(1L)
  result$remaining_values <- length(values) - position + 1L
  result
}

validate_nh_inputs <- function(control_file, gradient_tolerance = 0.001) {
  control <- read_nh_control(control_file)
  data_file <- file.path(dirname(control_file), control$data_file)
  data <- read_nh_data(
    data_file,
    first_year = control$first_year,
    last_year = control$last_year,
    plus_group = control$plus_group,
    n_cpue_series = control$n_cpue_series,
    n_survey_series = control$n_survey_series
  )

  expected_years <- seq.int(control$first_year, control$last_year)
  annual_years <- c(
    list(
      weight_start = data$weight_start[, 1L],
      weight_mid = data$weight_mid[, 1L],
      fishery_caa = data$fishery_caa[, 2L],
      cpue = data$cpue[, 1L],
      survey = data$survey[, 1L]
    ),
    setNames(
      lapply(data$survey_caa, \(x) x[, 2L]),
      paste0("survey_caa_", seq_along(data$survey_caa))
    )
  )

  mismatched_years <- names(annual_years)[
    !vapply(annual_years, identical, logical(1), as.numeric(expected_years))
  ]
  if (length(mismatched_years)) {
    stop(
      "Non-continuous or mismatched years in: ",
      paste(mismatched_years, collapse = ", ")
    )
  }
  if (control$selectivity_years[nrow(control$selectivity_years), 2L] !=
      control$last_year) {
    stop("Final selectivity period does not end in the terminal year.")
  }
  if (control$recruitment_years[[2L]] != control$last_year) {
    stop("Recruitment estimation does not end in the terminal year.")
  }
  if (control$ksp_years[nrow(control$ksp_years), 2L] != control$last_year) {
    stop("Final Ksp period does not end in the terminal year.")
  }
  if (control$eof != 54321) {
    stop("Control file EOF sentinel is not 54321.")
  }
  if (data$eof != 12345) {
    stop("Data file EOF sentinel is not 12345.")
  }
  if (data$remaining_values != 0L) {
    stop("Unexpected values remain after parsing the data file.")
  }

  invisible(list(
    control = control,
    data = data,
    gradient_tolerance = gradient_tolerance
  ))
}

read_admb_convergence <- function(run_directory) {
  par_file <- file.path(run_directory, "nh.par")
  header <- readLines(par_file, n = 1L, warn = FALSE)
  objective <- as.numeric(sub(
    ".*Objective function value = ([^ ]+).*",
    "\\1",
    header
  ))
  max_gradient <- as.numeric(sub(
    ".*Maximum gradient component = ([^ ]+).*",
    "\\1",
    header
  ))

  hessian_file <- file.path(run_directory, "admodel.hes")
  if (!file.exists(hessian_file)) {
    correlation_file <- file.path(run_directory, "nh.cor")
    if (!file.exists(correlation_file)) {
      stop("Neither admodel.hes nor nh.cor is available.")
    }
    correlation_header <- readLines(
      correlation_file,
      n = 1L,
      warn = FALSE
    )
    log_determinant <- as.numeric(sub(".*= *", "", correlation_header))
    if (!is.finite(log_determinant)) {
      stop("The archived Hessian log determinant is not finite.")
    }
    return(list(
      objective = objective,
      max_gradient = max_gradient,
      n_parameters = as.integer(sub(
        ".*Number of parameters = ([0-9]+).*",
        "\\1",
        header
      )),
      minimum_hessian_eigenvalue = NA_real_,
      positive_definite_hessian = TRUE,
      hessian_source = "Archived ADMB covariance output"
    ))
  }

  connection <- file(hessian_file, "rb")
  on.exit(close(connection))
  n_parameters <- readBin(connection, integer(), 1L, size = 4L)
  hessian <- matrix(
    readBin(connection, numeric(), n_parameters^2L, size = 8L),
    nrow = n_parameters,
    ncol = n_parameters
  )
  hessian <- (hessian + t(hessian)) / 2
  hessian_eigenvalues <- eigen(
    hessian,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  list(
    objective = objective,
    max_gradient = max_gradient,
    n_parameters = n_parameters,
    minimum_hessian_eigenvalue = min(hessian_eigenvalues),
    positive_definite_hessian = min(hessian_eigenvalues) > 0,
    hessian_source = "Binary ADMB Hessian"
  )
}

validate_nh_outputs <- function(run_directory, terminal_year,
                                gradient_tolerance = 0.001) {
  convergence <- read_admb_convergence(run_directory)
  if (!is.finite(convergence$objective)) {
    stop("The objective function is not finite.")
  }
  if (!is.finite(convergence$max_gradient) ||
      convergence$max_gradient > gradient_tolerance) {
    stop("The maximum gradient exceeds the convergence tolerance.")
  }
  if (!convergence$positive_definite_hessian) {
    stop("The ADMB Hessian is not positive definite.")
  }

  report <- read_rep(file.path(run_directory, "nh_R.rep"))
  report_years <- range(as.integer(report$Year))
  if (report_years[[2L]] != terminal_year) {
    stop("nh_R.rep does not end in the expected terminal year.")
  }

  output <- readr::read_csv(
    file.path(run_directory, "nh_out.csv"),
    show_col_types = FALSE
  )
  if (max(output$Year, na.rm = TRUE) != terminal_year) {
    stop("nh_out.csv does not end in the expected terminal year.")
  }

  invisible(list(
    convergence = convergence,
    report_years = report_years,
    output = output
  ))
}
