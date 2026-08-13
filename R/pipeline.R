#' Run the complete ProReT workflow
#'
#' @param series_matrix Local GEO Series Matrix.
#' @param basis_file User gene-by-program basis file.
#' @param output_dir Analysis output directory.
#' @param group_column,groups,target,reference Disease-design specification.
#' @param basis_gene_column,basis_gene_id_type,basis_type Basis specification.
#' @param disease_gene_id_type,annotation_package Disease row-ID mapping.
#' @param paired,block_column,covariates,normalization Disease-model settings.
#'   Raw RNA-seq counts are not accepted; count data require an upstream
#'   count-aware differential-expression model.
#' @param lincs_cache LINCS download/cache directory.
#' @param download_missing Download missing LINCS files.
#' @param gene_space,pert_type,time_hours,cell_lines,require_hq LINCS filters.
#' @param geometry,ridge_fraction,min_coverage,min_cell_lines,aggregator Model
#'   and ranking settings.
#' @return A \code{proret_result}.
#' @export
run_proret <- function(
    series_matrix, basis_file, output_dir,
    group_column, groups, target, reference,
    basis_gene_column = 1L, basis_gene_id_type = "SYMBOL",
    basis_type = "signed_template",
    disease_gene_id_type = "PROBEID", annotation_package = NULL,
    paired = FALSE, block_column = NULL, covariates = character(),
    normalization = "microarray_quantile",
    lincs_cache = default_lincs_cache(), download_missing = TRUE,
    gene_space = "BING", pert_type = "trt_cp", time_hours = 24,
    cell_lines = NULL, require_hq = FALSE,
    geometry = "orthogonal_subspace", ridge_fraction = 1e-3,
    min_coverage = 0.80, min_cell_lines = 3L, aggregator = "mean",
    cell_panel_policy = "coverage", min_cell_panel_coverage = 0.80,
    min_drug_panel_fraction = 0.60) {
  output_dir <- dir_create(output_dir)
  basis <- load_program_basis(
    basis_file, gene_column = basis_gene_column,
    gene_id_type = basis_gene_id_type, basis_type = basis_type
  )
  disease <- prepare_disease_from_series(
    series_matrix, group_column, groups, target, reference,
    gene_id_type = disease_gene_id_type,
    annotation_package = annotation_package, paired = paired,
    block_column = block_column, covariates = covariates,
    normalization = normalization,
    output_dir = file.path(output_dir, "disease")
  )
  paths <- tryCatch(
    check_lincs_files(lincs_cache),
    error = function(e) {
      if (!isTRUE(download_missing)) stop(e)
      download_lincs(lincs_cache)
    }
  )
  drug_programs <- project_lincs_gctx(
    paths[["gctx"]], paths[["sig_info"]], paths[["gene_info"]], basis,
    gene_space = gene_space, pert_type = pert_type,
    time_hours = time_hours, cell_lines = cell_lines, require_hq = require_hq
  )
  disease_program <- project_disease_signature(
    disease, basis, drug_programs$core_genes, min_coverage = min_coverage
  )
  result <- rank_drugs(
    disease_program, drug_programs, geometry, ridge_fraction,
    min_cell_lines, cell_lines, aggregator, cell_panel_policy,
    min_cell_panel_coverage, min_drug_panel_fraction
  )
  write_tsv(result$ranking, file.path(output_dir, "drug_ranking.tsv"))
  write_tsv(result$cell_scores, file.path(output_dir, "drug_cell_scores.tsv"))
  write_tsv(data.frame(
    program = names(disease_program$program_signature),
    score = disease_program$program_signature
  ), file.path(output_dir, "disease_program_signature.tsv"))
  manifest <- list(
    package_version = as.character(utils::packageVersion("ProReT")),
    run_time = format(Sys.time(), tz = "UTC"),
    series_matrix = normalizePath(series_matrix, winslash = "/"),
    series_sha256 = file_sha256(series_matrix),
    basis_file = normalizePath(basis_file, winslash = "/"),
    basis_sha256 = file_sha256(basis_file),
    basis_checksum = basis$audit$checksum,
    lincs_files = unname(paths),
    lincs_sha256 = vapply(paths, file_sha256, character(1L)),
    core_genes = length(drug_programs$core_genes),
    disease_coverage = disease_program$coverage,
    disease_normalization = disease$normalization,
    settings = result$settings,
    lincs_settings = drug_programs$settings
  )
  saveRDS(manifest, file.path(output_dir, "run_manifest.rds"))
  saveRDS(result, file.path(output_dir, "proret_result.rds"))
  result
}
