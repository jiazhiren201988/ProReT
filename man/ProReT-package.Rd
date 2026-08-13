\name{ProReT-package}
\alias{ProReT-package}
\alias{ProReT}
\title{Program-space Reversal of Transcription}
\description{
ProReT estimates moderated disease signatures from GEO Series Matrix files,
projects disease and LINCS drug perturbations through the same user-supplied
gene-by-program basis, and ranks compounds by inverse program alignment.
}
\details{
Signed cNMF score matrices are treated as program templates rather than
non-negative program usage. The default geometry applies regularized Gram
whitening before negative-cosine scoring.
}
\docType{package}
