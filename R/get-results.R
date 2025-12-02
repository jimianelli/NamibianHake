#' Fetch Model Results in Parallel
#'
#' This function fetches model results based on specified model names and directories.
#'
#' @param mod_names. A character vector of model names. Default is 'mod_names'.
#' @param rundir The main sub directory path for the models. Default is 'runs'
#' @param moddir The main directory path for the models. Default is 'mod_dir'
#' @param run_on_mac boolean if running on a mac
#' @param do_estimation boolean run estimation
#' @return A list containing model results.
#' @export
# mod_names. = mod_label; rundir = "mods"; moddir = mod_dir;
get_results <- function(mod_names. = mod_names,
                        rundir = "mods",
                        moddir = mod_dir,
                        run_on_mac = TRUE,
                        do_estimation = FALSE) {
  # Build absolute paths so calls work even when the working directory
  # is not the package root (e.g., when pkgdown evaluates vignettes).
  override_path <- getOption("NamibianHake.mods_path", NULL)
  if (is.null(override_path)) {
    env_override <- Sys.getenv("NAMIBIANHAKE_MODS_PATH", unset = NA)
    if (!is.na(env_override) && nzchar(env_override)) {
      override_path <- env_override
    }
  }
  if (!is.null(override_path)) {
    rundir_full <- override_path
  } else {
  is_absolute_path <- function(x) grepl("^(?:[A-Za-z]:)?[\\\\/]", x)
  find_description_root <- function(start = getwd()) {
    cur <- normalizePath(start, winslash = "/", mustWork = TRUE)
    repeat {
      if (file.exists(file.path(cur, "DESCRIPTION"))) {
        return(cur)
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
    NULL
  }

  if (is_absolute_path(rundir)) {
    rundir_full <- rundir
  } else {
    root <- find_description_root()
    if (is.null(root)) root <- getwd()
    rundir_full <- file.path(root, rundir)
  }
  }

  if (!run_on_mac) {
    # Normalize Windows-style paths while still allowing \\\\ separators
    rundir_full <- normalizePath(rundir_full, winslash = "\\\\", mustWork = FALSE)
  }
  fn <- file.path(rundir_full, moddir, "nh_R.rep")
  missing_rep <- fn[!file.exists(fn)]
  if (length(missing_rep)) {
    stop("Model output not found. Expected files: ", paste(missing_rep, collapse = ", "))
  }
  nmods <- length(mod_names.)

  num_cores <- parallel::detectCores()
  if (is.na(num_cores)) num_cores <- 1
  num_cores <- max(1, min(num_cores - 2, nmods))
  cl <- NULL
  if (num_cores > 1) {
    cl <- parallel::makeCluster(num_cores)
  }
  on.exit({
    if (!is.null(cl)) parallel::stopCluster(cl)
  })

  # Export necessary functions and objects to the cluster
  if (!is.null(cl)) {
    parallel::clusterExport(cl, c("read_fit", "read_admb", "read_rep", "run_nh"),
      envir = environment()
    )
  }

  run_parallel <- function(x, fun) {
    if (!is.null(cl)) {
      tryCatch(
        parallel::parLapply(cl = cl, X = x, fun = fun),
        error = function(e) {
          warning("Parallel execution failed (", conditionMessage(e), "); falling back to sequential.")
          parallel::stopCluster(cl)
          cl <<- NULL
          lapply(x, fun)
        }
      )
    } else {
      lapply(x, fun)
    }
  }

  # Run model parallel
  if (do_estimation) {
    system.time(modlst <- run_parallel(fn, run_nh))
  }
  # Fetch model results in parallel
  system.time(modlst <- run_parallel(fn, read_rep))
  fn <- file.path(rundir_full, moddir, "nh_out.csv")
  system.time(moddiag <- run_parallel(fn, readr::read_csv))
  names(modlst) <- mod_names.
  names(moddiag) <- mod_names.
  res <- list(modlst, moddiag)
  names(res) <- c("modlst", "moddiag")
  return(res)
}
