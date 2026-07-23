test_that("input comparison identifies the new assessment year", {
  skip_if_not(file.exists(test_path("..", "..", "mods", "m25", "nh.dat")))

  comparison <- compare_nh_inputs(
    test_path("..", "..", "mods", "m25", "nh.dat"),
    test_path("..", "..", "mods", "m26", "nh.dat")
  )

  expect_equal(comparison$previous_year, 2025)
  expect_equal(comparison$updated_year, 2026)
  expect_equal(unique(comparison$new_year$Year), 2026)
  expect_gt(nrow(comparison$historical), 0L)
})
