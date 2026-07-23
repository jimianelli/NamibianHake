#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
repository <- dirname(dirname(normalizePath(script_file)))
setwd(repository)

source(file.path("R", "read-admb2.R"))
source(file.path("R", "assessment-validation.R"))
source(file.path("R", "input-differences.R"))

for (model in c("m25", "m26")) {
  control <- file.path("mods", model, "nh.dat")
  result <- validate_nh_inputs(control)
  message(
    model,
    ": valid annual inputs through ",
    result$control$last_year
  )
}

cape_hakes <- read.csv(file.path("mods", "data", "CapeHakes.csv"))
by_species <- split(cape_hakes$year, cape_hakes$strata)
for (species in names(by_species)) {
  years <- by_species[[species]]
  stopifnot(
    identical(years, sort(unique(years))),
    identical(tail(years, 2L), 2025:2026)
  )
  message(species, ": survey inputs extend through 2026")
}

comparison <- compare_nh_inputs(
  file.path("mods", "m25", "nh.dat"),
  file.path("mods", "m26", "nh.dat")
)
report <- write_input_difference_report(
  comparison,
  file.path("mods", "data", "capehakes25.csv"),
  file.path("mods", "data", "CapeHakes.csv"),
  file.path("reports", "input_changes_2025_to_2026.html")
)
message("Wrote ", report)
