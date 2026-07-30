\name{basis}
\alias{validate_program_basis}
\alias{load_program_basis}
\alias{build_program_geometry}
\title{Load, validate and transform a program basis}
\usage{
validate_program_basis(basis,
  basis_type = c("signed_template", "nonnegative_spectra",
    "pathway_membership"), min_programs = 2L, min_genes = 20L)

load_program_basis(file, gene_column = 1L, gene_id_type = "SYMBOL",
  basis_type = "signed_template", annotation_package = NULL, sep = "\t")

build_program_geometry(basis, ridge_fraction = 1e-3,
  mode = c("orthogonal_subspace", "template"))
}
\arguments{
\item{basis}{A numeric gene-by-program matrix or a \code{proret_basis}.}
\item{basis_type}{Declared semantics of the basis weights.}
\item{min_programs, min_genes}{Minimum acceptable dimensions.}
\item{file}{Delimited basis file.}
\item{gene_column}{Gene identifier column name or position.}
\item{gene_id_type}{Input gene identifier type.}
\item{annotation_package}{Probe annotation package for PROBEID input.}
\item{sep}{Field separator.}
\item{ridge_fraction}{Ridge relative to the mean Gram diagonal.}
\item{mode}{Correlation-adjusted or raw-template geometry.}
}
\value{A validated matrix, \code{proret_basis}, or geometry diagnostics.}
\description{Audits a user basis and constructs its regularized geometry.}
