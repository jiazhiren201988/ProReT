test_that("all bases and checksum-matched references are correctly declared", {
  skip_on_cran()
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
    dims <- expected[[key]]
    basis <- load_builtin_basis(parts[[1L]], parts[[2L]])
    expect_equal(nrow(basis$weights), dims[[1L]], info = key)
    expect_equal(ncol(basis$weights), dims[[2L]], info = key)
    if (key == "K562_complete") {
      expect_error(
        load_builtin_reference(parts[[1L]], parts[[2L]]),
        "No checksum-matched"
      )
      next
    }
    reference <- load_builtin_reference(parts[[1L]], parts[[2L]])
    expect_equal(length(reference$drug_reference$core_genes), dims[[3L]],
                 info = key)
    expect_identical(
      basis_projection_checksum(reference$basis,
                                reference$drug_reference$core_genes),
      reference$drug_reference$projection_checksum, info = key
    )
  }
})

test_that("altered drug-program payloads fail integrity checks", {
  skip_on_cran()
  reference <- load_builtin_reference("K562", "common")
  files <- example_mm_files()
  ds <- data.table::fread(files[["disease_signature"]])
  signature <- stats::setNames(ds$moderated_t, as.character(ds$entrez_id))
  projection <- project_disease_signature(
    signature, reference$basis, reference$drug_reference$core_genes
  )
  altered <- reference$drug_reference
  program <- altered$program_names[[1L]]
  altered$signatures[[program]][[1L]] <-
    altered$signatures[[program]][[1L]] + 0.01
  expect_error(score_reversal(projection, altered), "payload integrity")
})

