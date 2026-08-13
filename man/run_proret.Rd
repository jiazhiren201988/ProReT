\name{run_proret}
\alias{run_proret}
\title{Run the complete ProReT workflow}
\usage{
run_proret(series_matrix, basis_file, output_dir, group_column, groups,
  target, reference, basis_gene_column = 1L,
  basis_gene_id_type = "SYMBOL", basis_type = "signed_template",
  disease_gene_id_type = "PROBEID", annotation_package = NULL,
  paired = FALSE, block_column = NULL, covariates = character(),
  normalization = "microarray_quantile",
  lincs_cache = default_lincs_cache(),
  download_missing = TRUE, gene_space = "BING", pert_type = "trt_cp",
  time_hours = 24, cell_lines = NULL, require_hq = FALSE,
  geometry = "orthogonal_subspace", ridge_fraction = 1e-3,
  min_coverage = 0.8, min_cell_lines = 3L, aggregator = "mean",
  cell_panel_policy = "coverage", min_cell_panel_coverage = 0.8,
  min_drug_panel_fraction = 0.6)
}
\arguments{
\item{series_matrix, basis_file, output_dir}{Input and output paths.}
\item{group_column, groups, target, reference}{Disease contrast definition.}
\item{basis_gene_column, basis_gene_id_type, basis_type}{Basis settings.}
\item{disease_gene_id_type, annotation_package}{Disease gene mapping.}
\item{paired, block_column, covariates, normalization}{Disease design settings.}
\item{lincs_cache, download_missing}{LINCS cache behavior.}
\item{gene_space, pert_type, time_hours, cell_lines, require_hq}{LINCS filters.}
\item{geometry, ridge_fraction, min_coverage, min_cell_lines, aggregator}{
Projection, scoring and aggregation settings.}
\item{cell_panel_policy, min_cell_panel_coverage, min_drug_panel_fraction}{
Disease-independent common-cell-panel settings.}
}
\value{A \code{proret_result}; files and provenance are written to
\code{output_dir}.}
\description{Runs disease modelling, LINCS projection and compound ranking.}
