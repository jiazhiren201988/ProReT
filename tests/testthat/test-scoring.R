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
  drugs$projection_checksum <- "not-the-same"
  expect_error(score_reversal(dp, drugs), "different core-basis")
  drugs <- toy_drugs(b, sig)
  drugs$core_genes <- rev(drugs$core_genes)
  expect_error(score_reversal(dp, drugs), "different core-basis|payload")
})

test_that("raw expression is rejected by the in-memory projection API", {
  b <- toy_basis()
  sig <- setNames(rnorm(40), as.character(1:40))
  x <- matrix(sig, ncol = 1, dimnames = list(names(sig), "S1"))
  md <- data.frame(
    sig_id = "S1", pert_id = "A", pert_iname = "drug",
    cell_iname = "CELL1", pert_type = "trt_cp", pert_itime = "24 h",
    pert_idose = "1 uM", is_hiq = 1
  )
  expect_error(
    project_lincs_matrix(x, md, b, input_level = "raw_expression"),
    "raw expression"
  )
})

test_that("dose, time and quality summaries are retained", {
  b <- toy_basis()
  sig <- setNames(rnorm(40), as.character(1:40))
  x <- cbind(S1 = sig, S2 = sig * 0.9)
  md <- data.frame(
    sig_id = c("S1", "S2"), pert_id = "A", pert_iname = "drug",
    cell_iname = "CELL1", pert_type = "trt_cp", pert_itime = "24 h",
    pert_idose = c("1 uM", "10 uM"), is_hiq = c(1, 0)
  )
  object <- project_lincs_matrix(x, md, b)
  expect_equal(object$signatures$n_doses, 2L)
  expect_match(object$signatures$doses, "1 uM")
  expect_equal(object$signatures$n_hq, 1L)
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
