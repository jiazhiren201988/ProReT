test_that("basis audit rejects invalid declarations and values", {
  w <- matrix(abs(rnorm(120)), 40, 3,
              dimnames = list(as.character(1:40), paste0("P", 1:3)))
  w[1, 1] <- -1
  expect_error(validate_program_basis(w, "nonnegative_spectra"), "Negative")
  w[1, 1] <- NA_real_
  expect_error(validate_program_basis(w), "NA/NaN/Inf")
  w[1, 1] <- 0
  w[, 1] <- 0
  expect_warning(
    expect_error(validate_program_basis(w), "All-zero"),
    "contains no negative"
  )
})

test_that("gene row order does not alter projections", {
  b1 <- toy_basis()
  ord <- sample(rownames(b1$weights))
  b2 <- toy_basis(ord)
  sig <- setNames(seq_len(40), as.character(seq_len(40)))
  p1 <- project_disease_signature(sig, b1, as.character(1:40))
  p2 <- project_disease_signature(sig, b2, as.character(1:40))
  expect_equal(p1$program_signature, p2$program_signature, tolerance = 1e-12)
})
