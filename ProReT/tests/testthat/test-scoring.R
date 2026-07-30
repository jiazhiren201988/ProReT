test_that("perfect reversal is first and same direction is last", {
  b <- toy_basis()
  sig <- setNames(scale(seq_len(40))[, 1], as.character(1:40))
  dp <- project_disease_signature(sig, b, as.character(1:40))
  drugs <- toy_drugs(b, sig)
  result <- rank_drugs(dp, drugs, min_cell_lines = 2)
  expect_equal(result$ranking$drug[[1]], "reverse")
  expect_equal(tail(result$ranking$drug, 1), "same")
})

test_that("basis and core mismatches stop scoring", {
  b <- toy_basis()
  sig <- setNames(rnorm(40), as.character(1:40))
  dp <- project_disease_signature(sig, b, as.character(1:40))
  drugs <- toy_drugs(b, sig)
  drugs$basis_checksum <- "not-the-same"
  expect_error(score_reversal(dp, drugs), "different basis")
  drugs <- toy_drugs(b, sig)
  drugs$core_genes <- rev(drugs$core_genes)
  expect_error(score_reversal(dp, drugs), "different core")
})

test_that("fixed toy data are deterministic", {
  b <- toy_basis()
  sig <- setNames(sin(seq_len(40)), as.character(1:40))
  dp1 <- project_disease_signature(sig, b, as.character(1:40))
  dp2 <- project_disease_signature(sig, b, as.character(1:40))
  expect_identical(dp1$program_signature, dp2$program_signature)
})

test_that("HQ filtering is strict rather than silently disabled", {
  b <- toy_basis()
  sig <- setNames(rnorm(40), as.character(1:40))
  x <- matrix(sig, ncol = 1,
              dimnames = list(names(sig), "S1"))
  md <- data.frame(
    sig_id = "S1", pert_id = "A", pert_iname = "drug",
    cell_iname = "CELL1", pert_type = "trt_cp", pert_itime = "24 h"
  )
  expect_error(
    project_lincs_matrix(x, md, b, require_hq = TRUE),
    "no quality field"
  )
})
