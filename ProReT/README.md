# ProReT

**ProReT** (Program-space Reversal of Transcription) ranks compounds by the
inverse alignment between a disease signature and LINCS perturbation
signatures in a shared, user-supplied gene-program representation.

ProReT is the reusable model implementation accompanying a Bioinformatics
Original Paper. Paper-only benchmarking (baseline/SOTA comparisons, RepoDB
validation, phenotype permutations and external cohorts) is intentionally
kept outside the package so that the public API remains focused and stable.

## Model

For a signed gene-level signature \(z\) and a gene-by-program basis \(W\),
ProReT first restricts both disease and drug data to the same frozen gene
universe and L2-normalizes each basis column:

\[
s = W^\top z.
\]

Correlated templates are adjusted using the regularized Gram geometry

\[
\widetilde{s}=(W^\top W+\lambda I)^{-1/2}s,
\]

and reversal is the negative cosine between disease and drug coordinates.
LINCS instances are summarized by the median within compound-cell conditions,
then reversal scores are averaged across at least three cell lines.

Because the supported cNMF `gene_spectra_score` matrices are signed, ProReT
calls these quantities **signed program-template alignments**, not cNMF usage
or non-negative program activity.

## Inputs

You provide:

1. a local GEO `*_series_matrix.txt` or `.txt.gz` file;
2. a tab-delimited gene-by-program basis produced by cNMF or another method;
3. explicit case/control regular expressions and, when necessary, pairing and
   covariates.

ProReT obtains the frozen LINCS Level 5 reference automatically. The GCTX file
is about 5.6 GB, so it is cached outside the installed package.

## Installation

```r
install.packages("BiocManager")
BiocManager::install(c(
  "AnnotationDbi", "Biobase", "GEOquery", "limma",
  "org.Hs.eg.db", "cmapR"
))
install.packages(c("data.table", "digest", "matrixStats", "testthat"))

# After publishing the repository:
# install.packages("remotes")
# remotes::install_github("YOUR_GITHUB_ACCOUNT/ProReT")
```

## Quick start: real MM data without the 5.6 GB download

The repository includes a compact, real-data example based on multiple
myeloma GSE6477, the K562 50-program basis and the precomputed 36-chunk LINCS
projection. It runs the actual projection, Gram correction, reversal scoring
and cross-cell aggregation, but does not need the original GCTX:

```r
library(ProReT)

mm <- run_mm_example("mm_quick_results")
head(mm$ranking, 20)
mm$settings$cell_panel
```

The installed paths are available for inspection:

```r
example_mm_files()
```

The bundled moderated-t signature makes the quick example independent of GEO
annotation packages. The original GSE6477 Series Matrix is also included so
users can reproduce disease-side preprocessing after installing `GEOquery`,
`limma` and `hgu133a.db`. See
`system.file("examples", "run_mm_from_series.R", package = "ProReT")`.

## Download LINCS once

```r
library(ProReT)
paths <- download_lincs()
paths
```

Interrupted downloads remain as `.part` files and are never treated as valid
data. File sizes are checked before analysis. For institutional storage:

```r
paths <- download_lincs("D:/references/ProReT_LINCS")
```

## Complete workflow

```r
library(ProReT)

result <- run_proret(
  series_matrix = "GSEXXXXX_series_matrix.txt.gz",
  basis_file = "gene_spectra_score.tsv",
  output_dir = "results/GSEXXXXX",

  # Inspect pData(read_geo_series(...)) before setting these fields.
  group_column = "characteristics_ch1",
  groups = c(Control = "^healthy", Disease = "^multiple myeloma"),
  target = "Disease",
  reference = "Control",

  basis_gene_column = "gene",
  basis_gene_id_type = "SYMBOL",
  basis_type = "signed_template",

  disease_gene_id_type = "PROBEID",
  annotation_package = "hgu133plus2.db",
  paired = FALSE,
  covariates = c("age", "sex"),

  lincs_cache = "D:/references/ProReT_LINCS",
  gene_space = "BING",
  pert_type = "trt_cp",
  time_hours = 24,
  require_hq = FALSE,

  geometry = "orthogonal_subspace",
  ridge_fraction = 1e-3,
  min_coverage = 0.80,
  min_cell_lines = 3,
  aggregator = "mean"
)

head(result$ranking)
```

`require_hq = TRUE` is strict: if the selected metadata lacks a quality field,
the analysis stops rather than silently disabling quality control.

## Output

`run_proret()` writes:

- `disease/disease_DEG.tsv`: limma moderated statistics;
- `disease_program_signature.tsv`: signed program alignment;
- `drug_cell_scores.tsv`: cell-specific reversal scores;
- `drug_ranking.tsv`: cross-cell compound ranking;
- `run_manifest.rds`: input checksums, basis checksum and all model settings;
- `proret_result.rds`: complete result object.

## Reproducibility contract

- Disease and drugs must use the same basis checksum and ordered core genes.
- Missing disease genes are set to zero only after the observed t statistics
  have been standardized; NA values are never interpreted as zero.
- The default coverage gate is 80%.
- Basis columns are normalized on the actual shared core universe.
- Projection caches include the basis contents, LINCS checksums and filter
  settings.
- Changing gene row order does not change the result.
- Cross-cell aggregation uses a disease-independent, coverage-derived common
  cell panel by default (80% compound coverage; each retained drug must cover
  at least 60% of the panel).

## Scope

The package produces drug rankings; it does not claim clinical efficacy.
Benchmarking and validation belong to the accompanying reproducibility
workflow, not the model API.

## License

MIT.
