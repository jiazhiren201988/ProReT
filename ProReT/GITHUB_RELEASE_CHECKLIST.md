# GitHub release checklist

Before making this repository public:

1. Replace `Author Name`, `author@example.org`, `USERNAME` and the placeholder
   citation with the final author and repository details.
2. Add the manuscript DOI when available.
3. Create a clean Git repository and commit the frozen `0.1.0` source.
4. Enable the included R-CMD-check GitHub Action.
5. Confirm that the Action passes with all suggested Bioconductor packages,
   including `cmapR`.
6. Create a signed `v0.1.0` tag and archive the release through Zenodo.
7. Record the release DOI in the manuscript Code Availability statement.
8. Do not upload the 5.6 GB LINCS GCTX file to GitHub; the package downloads it
   from GEO and caches it locally.

Local checks:

```r
devtools::test()
devtools::check()
```

The source tree also includes deterministic toy tests for score direction,
coordinate checksums, row-order invariance, missing values, coverage and strict
HQ-field handling.
