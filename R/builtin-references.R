#' List bundled program references
#'
#' @return A data frame describing the four bundled bases and their
#'   precomputed LINCS availability.
#' @export
available_builtin_references <- function() {
  data.frame(
    cell_line = c("K562", "K562", "KOLF", "KOLF"),
    variant = c("common", "complete", "common", "complete"),
    n_genes = c(5927L, 8139L, 5927L, 24781L),
    n_programs = c(50L, 50L, 30L, 30L),
    precomputed_lincs = c(TRUE, FALSE, TRUE, TRUE),
    reference_scope = c(
      "5927-gene common universe",
      "5927-gene BING intersection of the complete basis",
      "5927-gene common universe",
      "10129-gene BING intersection of the complete basis"
    ),
    stringsAsFactors = FALSE
  )
}

builtin_basis_filename <- function(cell_line, variant) {
  paste0(cell_line, "_program_",
         if (cell_line == "K562") "50" else "30", "_",
         variant, "_basis.tsv.gz")
}

#' Load a bundled K562 or KOLF program basis
#'
#' @param cell_line \code{"K562"} or \code{"KOLF"}.
#' @param variant \code{"common"} for the frozen 5,927-gene universe or
#'   \code{"complete"} for all mapped genes in the trained basis.
#' @param cache_dir Local reference cache.
#' @param download Download a missing reference from the versioned release.
#' @return A \code{proret_basis}.
#' @export
load_builtin_basis <- function(cell_line = c("K562", "KOLF"),
                               variant = c("common", "complete"),
                               cache_dir = default_reference_cache(),
                               download = TRUE) {
  cell_line <- match.arg(cell_line)
  variant <- match.arg(variant)
  path <- reference_file_path(
    builtin_basis_filename(cell_line, variant), cache_dir, download
  )
  load_program_basis(
    path, gene_column = "entrez_id", gene_id_type = "ENTREZID",
    basis_type = "signed_template"
  )
}

builtin_drug_filename <- function(cell_line, variant) {
  if (cell_line == "K562") {
    if (variant == "complete") return(NA_character_)
    return("K562_LINCS_drug_program_reference.rds")
  }
  paste0("KOLF_LINCS_", variant, "_drug_program_reference.rds")
}

#' Load bundled basis and matching precomputed LINCS projection
#'
#' The returned drug reference is checksum-bound to the returned basis. No
#' original GCTX download is required.
#'
#' @param cell_line \code{"K562"} or \code{"KOLF"}.
#' @param variant \code{"common"} or \code{"complete"}.
#' @param cache_dir Local reference cache.
#' @param download Download missing files from the versioned release.
#' @return A list with \code{basis} and \code{drug_reference}.
#' @export
load_builtin_reference <- function(cell_line = c("K562", "KOLF"),
                                   variant = c("common", "complete"),
                                   cache_dir = default_reference_cache(),
                                   download = TRUE) {
  cell_line <- match.arg(cell_line)
  variant <- match.arg(variant)
  drug_filename <- builtin_drug_filename(cell_line, variant)
  if (is.na(drug_filename)) {
    stop(
      "No checksum-matched precomputed LINCS projection is available for ",
      "K562 complete. Use K562 common, or rebuild a drug reference from the ",
      "complete basis with project_lincs_gctx().",
      call. = FALSE
    )
  }
  basis <- load_builtin_basis(cell_line, variant, cache_dir, download)
  path <- reference_file_path(
    drug_filename, cache_dir, download
  )
  drugs <- readRDS(path)
  if (!inherits(drugs, "proret_drug_programs")) {
    stop("Bundled LINCS reference has an invalid class.", call. = FALSE)
  }
  if (!all(drugs$core_genes %in% rownames(basis$weights))) {
    stop("Bundled LINCS core genes are absent from its basis.", call. = FALSE)
  }
  expected_programs <- colnames(basis$weights)
  if (!identical(drugs$program_names, expected_programs)) {
    stop("Bundled LINCS program coordinates do not match its basis.",
         call. = FALSE)
  }
  expected_projection <- basis_projection_checksum(basis, drugs$core_genes)
  if (is.null(drugs$projection_checksum) ||
      !identical(drugs$projection_checksum, expected_projection)) {
    stop("Bundled LINCS reference does not match the selected core basis.",
         call. = FALSE)
  }
  observed_payload <- drug_payload_checksum(
    drugs$signatures, drugs$core_genes, drugs$program_names,
    drugs$projection_checksum
  )
  if (is.null(drugs$payload_checksum) ||
      !identical(drugs$payload_checksum, observed_payload)) {
    stop("Bundled LINCS reference failed its payload integrity check.",
         call. = FALSE)
  }
  drugs$settings$selected_basis <- paste(cell_line, variant, sep = "_")
  list(basis = basis, drug_reference = drugs)
}
