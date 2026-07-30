\name{lincs-projection}
\alias{project_lincs_matrix}
\alias{project_lincs_gctx}
\title{Project perturbation signatures through a program basis}
\usage{
project_lincs_matrix(expression, metadata, basis, gene_space_ids = NULL,
  pert_type = "trt_cp", time_hours = 24, cell_lines = NULL,
  require_hq = FALSE)

project_lincs_gctx(gctx_file, sig_info_file, gene_info_file, basis,
  gene_space = "BING", pert_type = "trt_cp", time_hours = 24,
  cell_lines = NULL, require_hq = FALSE, chunk_size = 3000L,
  cache_dir = file.path(default_lincs_cache(),
    "projections"), overwrite = FALSE)
}
\arguments{
\item{expression}{Numeric gene-by-signature matrix.}
\item{metadata}{Signature metadata.}
\item{basis}{A \code{proret_basis}.}
\item{gene_space_ids}{Optional permitted gene IDs.}
\item{pert_type, time_hours, cell_lines, require_hq}{Perturbation filters.}
\item{gctx_file, sig_info_file, gene_info_file}{Frozen LINCS files.}
\item{gene_space}{\code{BING} or \code{ALL}.}
\item{chunk_size}{Signatures parsed per GCTX block.}
\item{cache_dir}{Projection cache directory.}
\item{overwrite}{Whether to replace an existing projection cache.}
}
\value{A \code{proret_drug_programs} object.}
\description{
Uses a single frozen core-gene universe, normalizes basis columns on that
universe and summarizes instances within compound-cell conditions by median.
}
