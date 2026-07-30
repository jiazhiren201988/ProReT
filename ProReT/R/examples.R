#' Paths to the bundled multiple-myeloma example
#'
#' @return Named paths to the GSE6477 Series Matrix, moderated disease
#'   signature, K562 basis and precomputed LINCS program reference.
#' @export
example_mm_files <- function() {
  keys <- c(
    series_matrix = "GSE6477_series_matrix.txt.gz",
    disease_signature = "MM_GSE6477_moderated_t.tsv",
    basis = "K562_program_50_common_basis.tsv.gz",
    drug_reference = "K562_LINCS_drug_program_reference.rds"
  )
  paths <- vapply(
    keys,
    function(x) system.file("extdata", x, package = "ProReT",
                            mustWork = TRUE),
    character(1L)
  )
  paths
}

#' Run the bundled real-data MM example
#'
#' This quick example uses the frozen moderated-t signature from GSE6477 and
#' the precomputed K562/LINCS program reference, so the 5.6 GB GCTX file is not
#' downloaded.
#'
#' @param output_dir Optional directory for the ranking and run manifest.
#' @return A \code{proret_result}.
#' @export
run_mm_example <- function(output_dir = NULL) {
  files <- example_mm_files()
  basis <- load_program_basis(
    files[["basis"]], gene_column = "entrez_id",
    gene_id_type = "ENTREZID", basis_type = "signed_template"
  )
  ds <- data.table::fread(files[["disease_signature"]])
  if (!all(c("entrez_id", "moderated_t") %in% names(ds))) {
    stop("Bundled MM disease signature has an invalid schema.",
         call. = FALSE)
  }
  signature <- stats::setNames(ds$moderated_t, as.character(ds$entrez_id))
  drugs <- readRDS(files[["drug_reference"]])
  if (!inherits(drugs, "proret_drug_programs")) {
    stop("Bundled drug reference has an invalid class.", call. = FALSE)
  }
  if (!identical(basis$audit$checksum, drugs$basis_checksum)) {
    stop("Bundled basis and drug reference checksums do not match.",
         call. = FALSE)
  }
  disease_program <- project_disease_signature(
    signature, basis, drugs$core_genes, min_coverage = 0.80
  )
  result <- rank_drugs(
    disease_program, drugs,
    geometry = "orthogonal_subspace", ridge_fraction = 1e-3,
    min_cell_lines = 3L, aggregator = "mean",
    cell_panel_policy = "coverage",
    min_cell_panel_coverage = 0.80,
    min_drug_panel_fraction = 0.60
  )
  if (!is.null(output_dir)) {
    output_dir <- dir_create(output_dir)
    write_tsv(result$ranking, file.path(output_dir, "MM_drug_ranking.tsv"))
    saveRDS(result, file.path(output_dir, "MM_proret_result.rds"))
    saveRDS(list(
      example = "GSE6477 MM versus normal donor",
      input_sha256 = vapply(files, file_sha256, character(1L)),
      basis_checksum = basis$audit$checksum,
      settings = result$settings
    ), file.path(output_dir, "MM_example_manifest.rds"))
  }
  result
}
