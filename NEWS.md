# ProReT 0.1.4

- Bound file-loaded program bases and precomputed drug references to the exact
  source-file SHA-256 plus their gene and program coordinates, providing a
  platform-independent integrity contract for distributed references.

# ProReT 0.1.3

- Bound projection references to canonical source basis weights rather than
  platform-computed normalized weights, eliminating BLAS-dependent checksum
  differences without changing the projection or scoring calculations.

# ProReT 0.1.2

- Canonicalized floating-point basis matrices before hashing so integrity
  checks are stable across BLAS implementations while remaining sensitive to
  substantive changes in program weights.

# ProReT 0.1.1

- Made content-derived reference checksums portable across operating systems
  and supported R releases by using the version-2 XDR serialization format.

# ProReT 0.1.0

- Initial basis-agnostic implementation.
- GEO Series Matrix disease modelling with moderated statistics.
- Automatic frozen LINCS Level 5 download and content-aware projection cache.
- Signed program-template projection with regularized Gram geometry.
- Drug-cell and cross-cell reversal ranking.
- Strict input, coordinate and checksum guards.
- Versioned complete K562 and KOLF2.1J reference downloads with SHA-256
  verification.
- Corrected NCBI LINCS archive URLs with official SHA-512 validation,
  resumable transfers and atomic decompression.
- Core-basis and drug-payload integrity checks that reject incompatible or
  modified precomputed references.
- Multi-platform GEO selection, explicit Level 5 input semantics,
  normalization provenance, and dose/time/quality summaries.
- Cross-platform continuous integration and expanded release documentation.
