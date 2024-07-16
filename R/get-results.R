#' Fetch Model Results in Parallel
#'
#' This function fetches model results based on specified model names and directories.
#'
#' @param mod_names. A character vector of model names. Default is 'mod_names'.
#' @param rundir The main sub directory path for the models. Default is 'runs'
#' @param moddir The main directory path for the models. Default is 'mod_dir'
#' @return A list containing model results.
#' @export
#mod_names. = mod_label; rundir = "mods"; moddir = mod_dir;
get_results <- function(mod_names. = mod_names,
                        rundir = "mods",
                        moddir = mod_dir,
                        run_on_mac = TRUE) {
  if (run_on_mac) {
    fn        <- paste0(rundir, "/", moddir, "/nh_R.rep")
  } else {
    #run on windows
    fn <- paste0(rundir, "\\", moddir, "\\nh_R.rep") #rundir is "C:\\GitProjects\\EBSpollock\\2023_runs\\"
  }
  nmods <- length(mod_names.)
  num_cores <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl))  # Ensure the cluster stops after function execution

  # Export necessary functions and objects to the cluster
  #source("R/read-admb2.R")
  parallel::clusterExport(cl, c("read_fit", "read_admb", "read_rep"), envir =
                            environment())

  # Fetch model results in parallel
  system.time(modlst <- parallel::parLapply(cl = cl, X = fn, fun = read_rep))
  fn <- paste0("mods/", moddir, "/nh_out.csv")
  system.time(moddiag <- parallel::parLapply(cl = cl, X = fn, fun = readr::read_csv))
  names(modlst) <- mod_names.
  names(moddiag) <- mod_names.
  res <- list(modlst, moddiag)
  names(res) <- c("modlst", "moddiag")
  return(res)
}
