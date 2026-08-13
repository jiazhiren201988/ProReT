# Bundled MM example provenance

The MIT license applies to ProReT source code only. Third-party GEO, LINCS and
program-reference data retain their original terms and citation requirements.
Users are responsible for confirming that redistribution and downstream use
comply with the source repositories' current policies.

The bundled example is provided for executable software demonstration and
reproducibility. It is not a new validation cohort.

- Disease: multiple myeloma versus normal donor plasma cells.
- GEO accession: GSE6477.
- Source URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6477
- Platform: GPL96 (Affymetrix Human Genome U133A Array).
- Disease statistic: limma moderated t, MM versus normal donor.
- Program representation: K562 signed `gene_spectra_score`, 50 programs,
  restricted to the frozen 5,927-gene K562/LINCS BING universe.
- Drug reference: GSE70138 LINCS Level 5, `trt_cp`, 24 h, all eligible cell
  lines; median within compound-cell conditions.
- LINCS source URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138
- Access date for the frozen source records: 13 August 2026.
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

## Bundled program systems

K562 is a signed 50-program dictionary; KOLF2.1J is a signed 30-program
dictionary derived from human induced pluripotent stem cells. Both are
interpreted as signed program templates rather than non-negative activity.

Four basis files are supplied:

- K562 complete: 8,139 mapped Entrez genes.
- K562 common: 5,927 genes shared with the frozen comparison universe.
- KOLF2.1J complete: 24,781 mapped Entrez genes.
- KOLF2.1J common: the same 5,927-gene comparison universe.

The K562 LINCS payload uses its 5,927-gene BING intersection and is compatible
with both K562 interfaces after content validation. KOLF provides separate
5,927-gene common and 10,129-gene native-BING projections. Every payload uses
GSE70138 Level 5, `trt_cp`, 24 h, and median aggregation within drug-cell
conditions.

The complete frozen basis and drug-program files are distributed as versioned
GitHub Release assets and are verified against the package manifest before
use. K562 and KOLF2.1J are complementary illustrative references for the same
software interface; their inclusion does not imply a universal ordering of
cell lines or bases.

Additional source SHA-256:

- K562 complete program model:
  `C3338D886D6E2DAA3A706C50EF8F76299BDB5102F8FF4BC7CBAA913C587C5D05`
- KOLF2.1J complete program model:
  `E25DAFB43ECB6C4F478034760AC936807F03A381C49075F6DC631AC3F7927583`
- KOLF2.1J native drug-program cache:
  `5675E6B25181CBD5400E15366B6DD2D748D8D79761B997916C3E5BB60DBAAF12`
- KOLF2.1J common basis source:
  `37C75930BFC6445E2BC1D9A999FCF0A71848219EB5ACC19EA19ACAD0ED3A1A2B`
- Recomputed KOLF2.1J common drug-program table:
  `5548C5FB15AD8C4B29C57E2822D5092F69FD6C235127B923A8B910182A5B8E2C`
