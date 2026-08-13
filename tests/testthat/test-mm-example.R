test_that("bundled MM example reproduces the frozen publication ranking", {
  files <- example_mm_files()
  expect_true(all(file.exists(files)))
  result <- run_mm_example()
  expect_s3_class(result, "proret_result")
  expect_equal(nrow(result$ranking), 1728L)
  expect_gt(result$disease_projection$coverage, 0.99)
  expect_equal(result$ranking$drug[[1L]], "as703026")
  expect_equal(result$ranking$reversal_score[[1L]], 0.459959205071598,
               tolerance = 1e-12)
  expect_equal(result$settings$cell_panel,
               sort(c("A375", "HA1E", "HELA", "HT29", "MCF7", "PC3",
                      "YAPC")))
})
