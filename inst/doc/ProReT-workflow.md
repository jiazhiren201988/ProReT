# ProReT workflow

## 1. Audit the GEO phenotype

```r
library(ProReT)
geo <- read_geo_series("GSEXXXXX_series_matrix.txt.gz")
colnames(geo$phenotype)
head(geo$phenotype)
```

Never infer case/control status from sample order. Define prespecified,
non-overlapping regular expressions.

## 2. Audit the program basis

```r
basis <- load_program_basis(
  "gene_spectra_score.tsv",
  gene_column = "gene",
  gene_id_type = "SYMBOL",
  basis_type = "signed_template"
)
basis
basis$audit
```

Declare a signed score matrix as `signed_template`. A genuinely non-negative
cNMF spectrum should be declared `nonnegative_spectra`.

## 3. Prepare the disease

```r
disease <- prepare_disease_from_series(
  "GSEXXXXX_series_matrix.txt.gz",
  group_column = "characteristics_ch1",
  groups = c(Control = "^healthy", Disease = "^disease"),
  target = "Disease", reference = "Control",
  gene_id_type = "PROBEID",
  annotation_package = "hgu133plus2.db"
)
```

For paired samples, set `paired = TRUE` and provide `block_column`. Covariates
must be named explicitly.

## 4. Project LINCS and rank drugs

```r
paths <- download_lincs()
drugs <- project_lincs_gctx(
  paths[["gctx"]], paths[["sig_info"]], paths[["gene_info"]],
  basis, gene_space = "BING", pert_type = "trt_cp", time_hours = 24
)
disease_program <- project_disease_signature(
  disease, basis, drugs$core_genes, min_coverage = 0.80
)
result <- rank_drugs(
  disease_program, drugs,
  geometry = "orthogonal_subspace",
  ridge_fraction = 1e-3,
  min_cell_lines = 3,
  aggregator = "mean"
)
head(result$ranking, 20)
```

The expensive drug projection is cached by basis checksum, LINCS checksums and
filter settings. A different basis necessarily generates a different cache.
