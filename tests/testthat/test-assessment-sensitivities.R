test_that("2026 sensitivity controls preserve valid terminal years", {
  skip_if_not(file.exists(test_path("..", "..", "mods", "m26", "nh.dat")))

  root <- withr::local_tempdir()
  paths <- create_annual_sensitivity_controls(
    test_path("..", "..", "mods", "m26", "nh.dat"),
    output_root = root
  )

  expect_named(
    paths,
    c(
      "m26_survey_q_1",
      "m26_survey_q_0_5",
      "m26_estimated_natural_mortality",
      "m26_dome_selectivity",
      "m26_steepness_0_5",
      "m26_steepness_0_9"
    )
  )
  for (path in paths) {
    control <- read_nh_control(path)
    expect_equal(control$last_year, 2026)
    expect_equal(control$selectivity_years[nrow(control$selectivity_years), 2L],
                 2026)
    expect_equal(control$recruitment_years[[2L]], 2026)
    expect_equal(control$ksp_years[nrow(control$ksp_years), 2L], 2026)
    expect_equal(control$eof, 54321)
  }
})
