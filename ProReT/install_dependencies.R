cran <- c("data.table", "digest", "curl", "testthat")
bioc <- c(
  "AnnotationDbi", "Biobase", "GEOquery", "limma", "org.Hs.eg.db",
  "hgu133a.db", "cmapR"
)

install.packages(setdiff(cran, rownames(installed.packages())))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(
  setdiff(bioc, rownames(installed.packages())),
  ask = FALSE, update = FALSE
)
