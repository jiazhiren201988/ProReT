\name{lincs-reference}
\alias{lincs_manifest}
\alias{check_lincs_files}
\alias{download_lincs}
\title{Manage the frozen LINCS Level 5 reference}
\usage{
lincs_manifest()

check_lincs_files(cache_dir = default_lincs_cache(),
  strict_size = TRUE)

download_lincs(cache_dir = default_lincs_cache(),
  overwrite = FALSE, verify = TRUE)
}
\arguments{
\item{cache_dir}{Directory holding frozen LINCS files.}
\item{strict_size, verify}{Whether expected file sizes must match.}
\item{overwrite}{Whether existing files are replaced.}
}
\value{A manifest data frame or named local file paths.}
\description{Downloads, caches and audits the external LINCS data resource.}
