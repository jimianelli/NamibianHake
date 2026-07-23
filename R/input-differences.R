compare_nh_inputs <- function(previous_control, updated_control) {
  previous <- validate_nh_inputs(previous_control)
  updated <- validate_nh_inputs(updated_control)
  previous_data <- previous$data
  updated_data <- updated$data
  common_rows <- seq_along(previous_data$years)

  compare_matrix <- function(previous_matrix, updated_matrix, block,
                             columns, year_column) {
    updated_matrix <- updated_matrix[common_rows, , drop = FALSE]
    changed <- which(previous_matrix != updated_matrix, arr.ind = TRUE)
    if (!nrow(changed)) {
      return(data.frame())
    }
    data.frame(
      Block = block,
      Year = previous_data$years[changed[, "row"]],
      Field = columns[changed[, "col"]],
      Previous = previous_matrix[changed],
      Updated = updated_matrix[changed],
      stringsAsFactors = FALSE
    ) |>
      subset(Field != year_column)
  }

  differences <- list(
    compare_matrix(
      previous_data$weight_start,
      updated_data$weight_start,
      "Start-year weight at age",
      c("Year", paste0("Age ", 0:8)),
      "Year"
    ),
    compare_matrix(
      previous_data$weight_mid,
      updated_data$weight_mid,
      "Mid-year weight at age",
      c("Year", paste0("Age ", 0:8)),
      "Year"
    ),
    data.frame(
      Block = "Catch",
      Year = previous_data$years,
      Field = "Catch",
      Previous = previous_data$catch,
      Updated = updated_data$catch[common_rows]
    ) |>
      subset(Previous != Updated),
    compare_matrix(
      previous_data$fishery_caa,
      updated_data$fishery_caa,
      "Fishery catch at age",
      c("Sample size", "Year", paste0("Age ", 0:8)),
      "Year"
    ),
    compare_matrix(
      previous_data$cpue,
      updated_data$cpue,
      "CPUE",
      c("Year", paste0("Series ", 1:6)),
      "Year"
    ),
    compare_matrix(
      previous_data$survey,
      updated_data$survey,
      "Survey",
      c("Year", "Survey 1", "Survey 1 CV", "Survey 2", "Survey 2 CV"),
      "Year"
    ),
    compare_matrix(
      previous_data$survey_caa[[1L]],
      updated_data$survey_caa[[1L]],
      "Survey 1 catch at age",
      c("Sample size", "Year", paste0("Age ", 0:8)),
      "Year"
    ),
    compare_matrix(
      previous_data$survey_caa[[2L]],
      updated_data$survey_caa[[2L]],
      "Survey 2 catch at age",
      c("Sample size", "Year", paste0("Age ", 0:8)),
      "Year"
    )
  )

  historical <- do.call(rbind, differences)
  rownames(historical) <- NULL
  new_year <- data.frame(
    Block = c(
      "Catch",
      paste0("CPUE series ", 1:6),
      "Survey 1",
      "Survey 1 CV",
      "Survey 2",
      "Survey 2 CV"
    ),
    Year = updated$control$last_year,
    Value = c(
      tail(updated_data$catch, 1L),
      tail(updated_data$cpue[, -1L], 1L),
      tail(updated_data$survey[, -1L], 1L)
    )
  )
  list(
    historical = historical,
    new_year = new_year,
    previous_year = previous$control$last_year,
    updated_year = updated$control$last_year
  )
}

write_input_difference_report <- function(comparison, previous_cape_file,
                                          updated_cape_file, output_file) {
  previous_cape <- read.csv(previous_cape_file)
  updated_cape <- read.csv(updated_cape_file)
  cape <- merge(
    previous_cape,
    updated_cape,
    by = c("strata", "year"),
    all = TRUE,
    suffixes = c("_previous", "_updated")
  )
  cape_changes <- cape[
    is.na(cape$biomass_previous) |
      is.na(cape$biomass_updated) |
      cape$biomass_previous != cape$biomass_updated |
      cape$cv_previous != cape$cv_updated, ,
    drop = FALSE
  ]

  summary <- aggregate(
    Field ~ Block,
    comparison$historical,
    length
  )
  names(summary)[[2L]] <- "Changed cells"
  table_html <- function(x, digits = 6L) {
    if (!nrow(x)) {
      return(htmltools::tags$p("No changes."))
    }
    htmltools::HTML(knitr::kable(
      x,
      format = "html",
      digits = digits,
      row.names = FALSE
    ))
  }

  document <- htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$title("Namibian hake input changes"),
      htmltools::tags$style(htmltools::HTML(
        "body{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem}table{border-collapse:collapse;width:100%;margin-bottom:2rem}th,td{border:1px solid #ccc;padding:.35rem;text-align:right}th:first-child,td:first-child{text-align:left}h1,h2{color:#17365d}"
      ))
    ),
    htmltools::tags$body(
      htmltools::tags$h1(
        paste0(
          "Input changes: ",
          comparison$previous_year,
          " to ",
          comparison$updated_year
        )
      ),
      htmltools::tags$p(
        "Generated from the parsed ADMB inputs. Historical changes are",
        "separated from values newly added for the update year."
      ),
      htmltools::tags$h2("Summary of historical revisions"),
      table_html(summary, 0L),
      htmltools::tags$h2("Historical values changed"),
      table_html(comparison$historical),
      htmltools::tags$h2(
        paste("Values added for", comparison$updated_year)
      ),
      table_html(comparison$new_year),
      htmltools::tags$h2("Species-specific survey changes"),
      table_html(cape_changes)
    )
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  htmltools::save_html(document, output_file, background = "white")
  invisible(output_file)
}
