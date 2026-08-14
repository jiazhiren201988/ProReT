# GitHub release checklist

Before publishing a release:

1. Run the complete test suite and `R CMD check` on the built source package.
2. Confirm that GitHub Actions pass on Windows, macOS and Linux, including the
   optional `cmapR` dependency.
3. Verify every file in `release-assets/v0.1.4` against
   `proret_reference_manifest()`.
4. Create the `v0.1.4` release and upload the seven frozen reference assets.
5. Confirm that `download_references()` succeeds in a clean cache.
6. Connect the repository to Zenodo and archive the tagged release.
7. Add the final software and article DOIs when they become available.

The original LINCS GCTX file is not uploaded to GitHub. ProReT downloads the
frozen Level 5 archive directly from NCBI GEO and validates the official
SHA-512 checksum before decompression.
