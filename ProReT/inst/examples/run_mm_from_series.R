library(ProReT)

files <- example_mm_files()
basis <- load_program_basis(
  files[["basis"]], gene_column = "entrez_id",
  gene_id_type = "ENTREZID", basis_type = "signed_template"
)
drugs <- readRDS(files[["drug_reference"]])

disease <- prepare_disease_from_series(
  files[["series_matrix"]],
  group_column = "title",
  groups = c(
    Disease = "^(New |Relapsed |Smoldering |rNew |MGUS )MM",
    Normal = "^Normal donor"
  ),
  target = "Disease", reference = "Normal",
  gene_id_type = "PROBEID",
  annotation_package = "hgu133a.db",
  normalization = "microarray_quantile"
)
disease_program <- project_disease_signature(
  disease, basis, drugs$core_genes, min_coverage = 0.80
)
result <- rank_drugs(
  disease_program, drugs,
  geometry = "orthogonal_subspace", ridge_fraction = 1e-3,
  min_cell_lines = 3, aggregator = "mean",
  cell_panel_policy = "coverage",
  min_cell_panel_coverage = 0.80,
  min_drug_panel_fraction = 0.60
)
head(result$ranking, 20)
