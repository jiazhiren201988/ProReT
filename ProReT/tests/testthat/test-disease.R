test_that("moderated disease signature has the expected direction", {
  skip_if_not_installed("limma")
  set.seed(9)
  x <- matrix(rnorm(40 * 8), 40, 8,
              dimnames = list(as.character(1:40), paste0("S", 1:8)))
  x[1:5, 5:8] <- x[1:5, 5:8] + 3
  md <- data.frame(sample = colnames(x),
                   group = rep(c("Control", "Disease"), each = 4))
  fit <- fit_disease_signature(
    x, md, target = "Disease", reference = "Control",
    normalization = "none"
  )
  expect_s3_class(fit, "proret_disease")
  expect_gt(mean(fit$signature[as.character(1:5)]), 0)
})

test_that("NA is never interpreted as zero", {
  b <- toy_basis()
  sig <- setNames(rnorm(40), as.character(1:40))
  sig[1] <- NA_real_
  expect_error(project_disease_signature(sig, b, as.character(1:40)), "NA")
})

test_that("coverage gate is enforced", {
  b <- toy_basis()
  sig <- setNames(rnorm(25), as.character(1:25))
  expect_error(project_disease_signature(sig, b, as.character(1:40)),
               "coverage")
})
