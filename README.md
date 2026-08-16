# ProReT

[![R-CMD-check](https://github.com/jiazhiren201988/ProReT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jiazhiren201988/ProReT/actions/workflows/R-CMD-check.yaml)
[![GitHub release](https://img.shields.io/github/v/release/jiazhiren201988/ProReT)](https://github.com/jiazhiren201988/ProReT/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21943239.svg)](https://doi.org/10.5281/zenodo.21943239)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)

**ProReT** (Program-space Reversal of Transcription) ranks compounds by the
inverse alignment between a disease signature and LINCS perturbation
signatures in a shared, user-supplied gene-program representation.

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

install.packages("remotes")
remotes::install_github("jiazhiren201988/ProReT")
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

## Built-in K562 and KOLF program systems

K562 and KOLF2.1J are provided as equal, complementary examples of the
basis-agnostic interface. The compact MM reference is installed with the
package; the complete versioned reference collection is obtained once with:

```r
download_references()
```

| System | Variant | Genes | Programs | Precomputed LINCS scope |
|---|---:|---:|---:|---|
| K562 | common | 5,927 | 50 | 5,927 shared genes |
| K562 | complete | 8,139 | 50 | rebuild required from LINCS |
| KOLF2.1J | common | 5,927 | 30 | 5,927 shared genes |
| KOLF2.1J | complete | 24,781 | 30 | 10,129 KOLF-LINCS BING genes |

```r
available_builtin_references()

k562 <- load_builtin_reference("K562", "common")
kolf  <- load_builtin_reference("KOLF", "common")

# Each object contains a checksum-matched pair:
k562$basis
k562$drug_reference
```

To analyse the same named disease signature with either reference:

```r
kolf_disease <- project_disease_signature(
  disease_signature,
  kolf$basis,
  kolf$drug_reference$core_genes
)
kolf_result <- rank_drugs(kolf_disease, kolf$drug_reference)
```

Use `variant = "common"` when a matched gene universe is scientifically
important. Use `variant = "complete"` to deploy each trained basis with its
native LINCS BING intersection. A precomputed K562-complete drug reference is
not supplied because its weights on the frozen core differ from the K562-common
basis; it must be rebuilt with `project_lincs_gctx()`. K562 and KOLF are examples of the general
program-space method, not claims that either cell line is universally optimal.

## Download LINCS once

```r
library(ProReT)
paths <- download_lincs()
paths
```

Interrupted LINCS downloads remain as `.part` files and are never treated as
valid data. The official NCBI file size and SHA-512 checksum are verified
before decompression. For institutional storage:

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



## Data and licensing

The MIT license covers ProReT source code. GEO, LINCS and program-reference
data retain the terms of their source repositories. File versions, source
URLs and checksums are documented in `inst/extdata/DATA_PROVENANCE.md` and in
the versioned release manifest.

## Citation and support

Use `citation("ProReT")` for the current citation record. The archived v0.1.6
release is available at https://doi.org/10.5281/zenodo.21943239. Reproducible
bug reports can be submitted at
https://github.com/jiazhiren201988/ProReT/issues. Contribution and security
policies are provided in `CONTRIBUTING.md` and `SECURITY.md`.

## License

MIT.
