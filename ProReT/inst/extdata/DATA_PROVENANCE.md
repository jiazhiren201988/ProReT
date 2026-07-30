# Bundled MM example provenance

The bundled example is provided for executable software demonstration and
reproducibility. It is not a new validation cohort.

- Disease: multiple myeloma versus normal donor plasma cells.
- GEO accession: GSE6477.
- Platform: GPL96 (Affymetrix Human Genome U133A Array).
- Disease statistic: limma moderated t, MM versus normal donor.
- Program representation: K562 signed `gene_spectra_score`, 50 programs,
  restricted to the frozen 5,927-gene K562/LINCS BING universe.
- Drug reference: GSE70138 LINCS Level 5, `trt_cp`, 24 h, all eligible cell
  lines; median within compound-cell conditions.
- Program geometry: regularized Gram whitening, ridge fraction 0.001.
- Reversal metric: negative cosine.
- Cell aggregation: coverage-derived common panel (cell coverage >= 0.80);
  retained compounds cover >= 0.60 of the panel and at least three cells.

Source SHA-256:

- GSE6477 Series Matrix:
  `255756A266D2740664D11580F76DB2BA64C6FB0C0C184DD31C0D535F3C7FF896`
- Uncompressed K562 common basis:
  `3DE9398CF0F2B90521A595A1522EB4D0608C4CB5206E42905E66183CCE51F37A`
- Original 36-chunk drug-program cache:
  `BFA7A03E66C5EA09C1A0ECA49F7CCBD0DD7042130AB1D07A7167872B97DD823D`
- GSE6477 moderated-t table:
  `BE70965053299EAD9F92ADE4E00D4015092607890213EA116643799DA91E3D3B`
