\name{mm-example}
\alias{example_mm_files}
\alias{run_mm_example}
\title{Bundled real-data multiple-myeloma example}
\usage{
example_mm_files()

run_mm_example(output_dir = NULL)
}
\arguments{
\item{output_dir}{Optional directory for ranking and provenance outputs.}
}
\value{
\code{example_mm_files()} returns named installed paths.
\code{run_mm_example()} returns a \code{proret_result}.
}
\description{
Runs GSE6477 multiple myeloma against a precomputed K562/LINCS program
reference without downloading the full 5.6 GB GCTX matrix.
}
