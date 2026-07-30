\name{scoring}
\alias{project_disease_signature}
\alias{score_reversal}
\alias{aggregate_drug_scores}
\alias{rank_drugs}
\title{Project disease signatures and rank reversal compounds}
\usage{
project_disease_signature(disease, basis, core_genes,
  min_coverage = 0.8, min_observed = 20L)

score_reversal(disease_projection, drug_programs,
  geometry = c("orthogonal_subspace", "template"),
  ridge_fraction = 1e-3)

aggregate_drug_scores(cell_scores, min_cell_lines = 3L,
  cell_panel = NULL, aggregator = c("mean", "median"),
  cell_panel_policy = c("coverage", "all", "configured"),
  min_cell_panel_coverage = 0.8, min_drug_panel_fraction = 0.6)

rank_drugs(disease_projection, drug_programs,
  geometry = "orthogonal_subspace", ridge_fraction = 1e-3,
  min_cell_lines = 3L, cell_panel = NULL, aggregator = "mean",
  cell_panel_policy = "coverage", min_cell_panel_coverage = 0.8,
  min_drug_panel_fraction = 0.6)
}
\arguments{
\item{disease}{A \code{proret_disease} or named moderated-statistic vector.}
\item{basis}{A \code{proret_basis}.}
\item{core_genes}{Frozen shared gene universe.}
\item{min_coverage, min_observed}{Disease coverage gates.}
\item{disease_projection}{A projected disease object.}
\item{drug_programs}{A projected perturbation reference.}
\item{geometry, ridge_fraction}{Program-space geometry settings.}
\item{cell_scores}{Cell-specific reversal scores.}
\item{min_cell_lines}{Minimum cells per compound.}
\item{cell_panel}{Optional fixed eligible cell panel.}
\item{aggregator}{Cross-cell mean or median.}
\item{cell_panel_policy}{Coverage-derived, all-cell or configured panel.}
\item{min_cell_panel_coverage}{Minimum drug coverage used to select a cell.}
\item{min_drug_panel_fraction}{Minimum selected-panel fraction required per
drug.}
}
\value{A disease projection, cell-score table, ranking table or
\code{proret_result}.}
\description{
Computes signed template alignment and negative-cosine reversal in a shared,
checksum-guarded coordinate system.
}
