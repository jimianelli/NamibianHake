#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", arguments[grepl("^--file=", arguments)])
repository <- dirname(dirname(normalizePath(script_file)))
setwd(repository)

quarto <- Sys.which("quarto")
if (!nzchar(quarto)) {
  stop("Quarto is not installed or is not on PATH.")
}

report <- file.path(
  "vignettes",
  "02-Namibian_hake_model_2026.qmd"
)
for (format in c("html", "pdf")) {
  message("Rendering ", format)
  status <- system2(quarto, c("render", report, "--to", format))
  if (status != 0L) {
    stop("Failed to render the ", format, " report.")
  }
}
