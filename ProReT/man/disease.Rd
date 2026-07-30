\name{disease}
\alias{read_geo_series}
\alias{fit_disease_signature}
\alias{prepare_disease_from_series}
\title{Construct a moderated disease signature}
\usage{
read_geo_series(series_matrix)

fit_disease_signature(expression, sample_metadata, sample_column = "sample",
  group_column = "group", target, reference, gene_id_type = "ENTREZID",
  annotation_package = NULL, paired = FALSE, block_column = NULL,
  covariates = character(), normalization = "microarray_quantile",
  force_log2 = NA)

prepare_disease_from_series(series_matrix, group_column, groups, target,
  reference, gene_id_type = "PROBEID", annotation_package = NULL,
  paired = FALSE, block_column = NULL, covariates = character(),
  normalization = "microarray_quantile", force_log2 = NA,
  output_dir = NULL)
}
\arguments{
\item{series_matrix}{Local GEO Series Matrix text file, optionally gzipped.}
\item{expression}{Numeric gene-by-sample matrix.}
\item{sample_metadata}{Data frame with one row per sample.}
\item{sample_column, group_column}{Sample and group metadata columns.}
\item{groups}{Named, mutually exclusive regular expressions.}
\item{target, reference}{Disease and control group labels.}
\item{gene_id_type}{Expression-row identifier type.}
\item{annotation_package}{Annotation package for probe mapping.}
\item{paired}{Whether the design is paired.}
\item{block_column}{Pairing identifier column.}
\item{covariates}{Explicit nuisance covariate columns.}
\item{normalization}{One of \code{microarray_quantile}, \code{log2_only},
or \code{none}.}
\item{force_log2}{Override automatic log-scale detection.}
\item{output_dir}{Optional output directory.}
}
\value{A GEO input list or a \code{proret_disease} object.}
\description{Fits an explicit limma design and returns moderated t statistics.}
