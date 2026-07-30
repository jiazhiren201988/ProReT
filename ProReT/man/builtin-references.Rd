\name{builtin-references}
\alias{available_builtin_references}
\alias{load_builtin_basis}
\alias{load_builtin_reference}
\title{Bundled K562 and KOLF program references}
\usage{
available_builtin_references()

load_builtin_basis(cell_line = c("K562", "KOLF"),
  variant = c("common", "complete"))

load_builtin_reference(cell_line = c("K562", "KOLF"),
  variant = c("common", "complete"))
}
\arguments{
\item{cell_line}{\code{K562} or \code{KOLF}.}
\item{variant}{The matched 5,927-gene \code{common} basis or the mapped
\code{complete} trained basis.}
}
\value{
A reference catalogue, a \code{proret_basis}, or a list containing the basis
and checksum-matched precomputed LINCS drug reference.
}
\description{
Provides direct access to the two illustrative program systems without
requiring users to reconstruct their bases or download the original GCTX.
}
