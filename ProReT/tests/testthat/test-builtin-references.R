test_that("all four bundled bases and matching references load", {
  catalogue <- available_builtin_references()
  expect_equal(nrow(catalogue), 4L)
  expected <- list(
    K562_common = c(5927L, 50L, 5927L),
    K562_complete = c(8139L, 50L, 5927L),
    KOLF_common = c(5927L, 30L, 5927L),
    KOLF_complete = c(24781L, 30L, 10129L)
  )
  for (key in names(expected)) {
    parts <- strsplit(key, "_", fixed = TRUE)[[1L]]
    reference <- load_builtin_reference(parts[[1L]], parts[[2L]])
    dims <- expected[[key]]
    expect_equal(nrow(reference$basis$weights), dims[[1L]], info = key)
    expect_equal(ncol(reference$basis$weights), dims[[2L]], info = key)
    expect_equal(length(reference$drug_reference$core_genes), dims[[3L]],
                 info = key)
    expect_identical(reference$basis$audit$checksum,
                     reference$drug_reference$basis_checksum, info = key)
  }
})

test_that("K562 common and complete interfaces preserve the frozen ranking", {
  files <- example_mm_files()
  ds <- data.table::fread(files[["disease_signature"]])
  signature <- stats::setNames(ds$moderated_t, as.character(ds$entrez_id))
  common <- load_builtin_reference("K562", "common")
  complete <- load_builtin_reference("K562", "complete")
  p1 <- project_disease_signature(
    signature, common$basis, common$drug_reference$core_genes
  )
  p2 <- project_disease_signature(
    signature, complete$basis, complete$drug_reference$core_genes
  )
  expect_equal(p1$program_signature, p2$program_signature, tolerance = 1e-12)
  r1 <- rank_drugs(p1, common$drug_reference)
  r2 <- rank_drugs(p2, complete$drug_reference)
  expect_identical(r1$ranking$drug, r2$ranking$drug)
  expect_equal(r1$ranking$reversal_score, r2$ranking$reversal_score,
               tolerance = 1e-12)
})
