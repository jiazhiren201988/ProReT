test_that("reference manifest is complete and versioned", {
  manifest <- proret_reference_manifest()
  expect_equal(nrow(manifest), 7L)
  expect_true(all(grepl("/v0.1.4/", manifest$url, fixed = TRUE)))
  expect_true(all(nchar(manifest$sha256) == 64L))
  expect_equal(anyDuplicated(manifest$filename), 0L)
})

test_that("official LINCS manifest uses existing dated archives", {
  manifest <- lincs_manifest()
  expect_true(all(grepl("[.]gz$", manifest$archive_filename)))
  expect_true(all(grepl("GSE70138", manifest$url, fixed = TRUE)))
  expect_true(all(nchar(manifest$archive_sha512) == 128L))
})
