#' Run NH Model
#'
#' This function runs the NH model with specified parameters and returns the output.
#'
#' @param model A character string specifying the model directory. Default is "m1".
#' @param steepness A numeric value for the steepness parameter. Default is 0.7.
#' @param runit A logical value indicating whether to run the model. Default is TRUE.
#'
#' @return The function returns the output read from the `nh_R.rep` file in the specified model directory.
#'
#' @examples
#' \dontrun{
#' m2 <- run_nh("m2")
#' }
#'
#' @export
model <- "m1"
steepness <- .7

run_nh <- function(model = "m1", runit = TRUE) {
  if (runit) {
    run_dir <- here("mods", model)
    arg <- " -iprint 400 "
    setwd(run_dir)
    if (.Platform$OS.type == "windows") {
      run <- paste0("..\\..\\source\\nh", arg)
      shell(run)
    } else {
      run <- paste0("../../source/nh", arg)
      system(run)
    }
    setwd(here())
  }
  out <- read_rep(here("mods", model, "nh_R.rep"))
  return(out)
}
