test_that("2025 and 2026 inputs have valid dimensions and terminal years", {
  skip_if_not(file.exists(test_path("..", "..", "mods", "m25", "nh.dat")))

  for (model in c("m25", "m26")) {
    result <- validate_nh_inputs(
      test_path("..", "..", "mods", model, "nh.dat")
    )
    expected_year <- if (model == "m25") 2025 else 2026

    expect_equal(result$control$last_year, expected_year)
    expect_equal(result$control$eof, 54321)
    expect_equal(result$data$eof, 12345)
    expect_equal(result$data$remaining_values, 0L)
    expect_equal(nrow(result$data$weight_start), expected_year - 1963L)
    expect_equal(ncol(result$data$weight_start), 10L)
    expect_equal(ncol(result$data$fishery_caa), 11L)
    expect_equal(ncol(result$data$cpue), 7L)
    expect_equal(ncol(result$data$survey), 5L)
  }
})

test_that("Cape hake biomass inputs are ordered and extend through 2026", {
  skip_if_not(file.exists(
    test_path("..", "..", "mods", "data", "CapeHakes.csv")
  ))

  biomass <- read.csv(
    test_path("..", "..", "mods", "data", "CapeHakes.csv")
  )

  by_species <- split(biomass$year, biomass$strata)
  for (years in by_species) {
    expect_equal(years, sort(unique(years)))
    expect_equal(tail(years, 2L), 2025:2026)
  }
})

test_that("available assessment outputs satisfy convergence checks", {
  skip_if_not(file.exists(test_path("..", "..", "mods", "m25")))

  for (model in c("m25", "m26")) {
    run_directory <- test_path("..", "..", "mods", model)
    if (!file.exists(file.path(run_directory, "admodel.hes"))) {
      next
    }
    terminal_year <- if (model == "m25") 2025 else 2026

    output <- validate_nh_outputs(run_directory, terminal_year)
    expect_lt(output$convergence$max_gradient, 0.001)
    expect_gt(output$convergence$minimum_hessian_eigenvalue, 0)
  }
})
