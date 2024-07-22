#' Get Relative Biomass
#'
#' This function retrieves and processes biomass data for a specified species and end year.
#'
#' @param spp A character string specifying the species. Default is "paradoxus".
#' @param endyr An integer specifying the end year for the data. Default is 2025.
#'
#' @return A tidy data frame with the relative biomass for the specified species and end year.
#'
#' @examples
#' get_relative_biomass(spp = "paradoxus", endyr = 2025)
#'
#' @import dplyr
#' @import readr
#' @import rema
#' @import here
#'
#' @export
get_relative_biomass <- function(spp = "paradoxus", endyr = 2025) {
  data <- rema::prepare_rema_input(
    model_name = spp,
    multi_survey = 0,
    admb_re = NULL,
    biomass_dat = read_csv(here("mods", "data", paste0(spp, ".csv"))) |>
      mutate(strata = spp),
    cpue_dat = NULL, sum_cpue_index = FALSE,
    start_year = NULL, end_year = endyr,
    wt_biomass = NULL, wt_cpue = NULL,
    PE_options = NULL, q_options = NULL,
    zeros = NULL, extra_biomass_cv = NULL, extra_cpue_cv = NULL
  )
  m <- rema::fit_rema(data)
  output <- rema::tidy_rema(rema_model = m)
  return(output)
}
