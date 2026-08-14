#' Validate a gene-by-program basis
#'
#' @param basis Numeric matrix with genes in rows and programs in columns.
#' @param basis_type One of \code{"signed_template"}, \code{"nonnegative_spectra"}
#'   or \code{"pathway_membership"}.
#' @param min_programs Minimum number of programs.
#' @param min_genes Minimum number of unique genes.
#' @return A validated matrix with audit attributes.
#' @export
validate_program_basis <- function(
    basis,
    basis_type = c("signed_template", "nonnegative_spectra", "pathway_membership"),
    min_programs = 2L,
    min_genes = 20L) {
  basis_type <- match.arg(basis_type)
  basis <- as.matrix(basis)
  storage.mode(basis) <- "double"
  if (is.null(rownames(basis)) || any(!nzchar(rownames(basis)))) {
    stop("Basis rows must have non-empty gene identifiers.", call. = FALSE)
  }
  if (is.null(colnames(basis)) || any(!nzchar(colnames(basis)))) {
    stop("Basis columns must have non-empty program names.", call. = FALSE)
  }
  if (anyDuplicated(colnames(basis))) stop("Program names must be unique.", call. = FALSE)
  if (any(!is.finite(basis))) stop("Basis contains NA/NaN/Inf.", call. = FALSE)
  if (ncol(basis) < min_programs) stop("Too few programs.", call. = FALSE)
  if (basis_type != "signed_template" && any(basis < 0)) {
    stop("Negative weights are incompatible with basis_type='", basis_type, "'.",
         call. = FALSE)
  }
  if (basis_type == "signed_template" && !any(basis < 0)) {
    warning("The declared signed basis contains no negative weights.", call. = FALSE)
  }
  ids <- trimws(as.character(rownames(basis)))
  keep <- nzchar(ids)
  basis <- basis[keep, , drop = FALSE]
  ids <- ids[keep]
  if (anyDuplicated(ids)) {
    basis <- limma::avereps(basis, ID = ids)
  }
  if (nrow(basis) < min_genes) stop("Too few unique basis genes.", call. = FALSE)
  zero <- colSums(abs(basis)) == 0
  if (any(zero)) stop("All-zero programs: ", paste(colnames(basis)[zero], collapse = ", "),
                       call. = FALSE)
  qr_rank <- qr(basis)$rank
  gram <- crossprod(normalize_columns(basis))
  ev <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
  attr(basis, "proret_audit") <- list(
    basis_type = basis_type,
    n_genes = nrow(basis),
    n_programs = ncol(basis),
    signed_fraction = mean(basis < 0),
    matrix_rank = qr_rank,
    effective_rank = sum(ev)^2 / sum(ev^2),
    checksum = hash_numeric_matrix(basis, digits = 15L)
  )
  basis
}

map_ids_to_entrez <- function(ids, id_type, annotation_package = NULL) {
  id_type <- toupper(id_type)
  ids <- sub("\\..*$", "", trimws(as.character(ids)))
  if (id_type == "ENTREZID") return(ids)
  if (id_type == "PROBEID") {
    if (is.null(annotation_package)) {
      stop("annotation_package is required for PROBEID mapping.", call. = FALSE)
    }
    if (!requireNamespace(annotation_package, quietly = TRUE)) {
      stop("Install annotation package '", annotation_package, "'.", call. = FALSE)
    }
    db <- getExportedValue(annotation_package, annotation_package)
  } else {
    require_packages(c("AnnotationDbi", "org.Hs.eg.db"),
                     "human gene-ID mapping")
    db <- org.Hs.eg.db::org.Hs.eg.db
  }
  require_packages("AnnotationDbi", "gene-ID mapping")
  if (!id_type %in% AnnotationDbi::keytypes(db)) {
    stop("Unsupported key type ", id_type, " for the selected annotation database.",
         call. = FALSE)
  }
  unname(AnnotationDbi::mapIds(
    db, keys = ids, column = "ENTREZID", keytype = id_type,
    multiVals = "first"
  ))
}

#' Load a user-supplied program basis
#'
#' @param file Delimited text file containing genes by programs.
#' @param gene_column Gene identifier column name or position.
#' @param gene_id_type \code{ENTREZID}, \code{SYMBOL}, \code{ENSEMBL}, or
#'   \code{PROBEID}.
#' @param basis_type Declared basis semantics.
#' @param annotation_package Required for probe identifiers.
#' @param sep Field separator; tab is the default.
#' @return A \code{proret_basis} object.
#' @export
load_program_basis <- function(
    file,
    gene_column = 1L,
    gene_id_type = "SYMBOL",
    basis_type = "signed_template",
    annotation_package = NULL,
    sep = "\t") {
  file <- assert_file(file, "program basis")
  if (grepl("\\.gz$", file, ignore.case = TRUE)) {
    con <- gzfile(file, open = "rt")
    on.exit(close(con), add = TRUE)
    x <- utils::read.delim(
      con, sep = sep, check.names = FALSE, stringsAsFactors = FALSE
    )
  } else {
    x <- data.table::fread(
      file, sep = sep, data.table = FALSE, check.names = FALSE
    )
  }
  if (is.numeric(gene_column)) gene_column <- names(x)[gene_column]
  if (!gene_column %in% names(x)) stop("gene_column not found.", call. = FALSE)
  ids <- map_ids_to_entrez(x[[gene_column]], gene_id_type, annotation_package)
  mat <- as.matrix(x[, setdiff(names(x), gene_column), drop = FALSE])
  storage.mode(mat) <- "double"
  keep <- !is.na(ids) & nzchar(ids)
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- ids[keep]
  mat <- validate_program_basis(mat, basis_type = basis_type)
  out <- list(
    weights = mat,
    basis_type = basis_type,
    source_file = file,
    source_sha256 = file_sha256(file),
    audit = attr(mat, "proret_audit")
  )
  class(out) <- "proret_basis"
  out
}

#' @export
print.proret_basis <- function(x, ...) {
  cat("<proret_basis>\n",
      " genes: ", nrow(x$weights), "\n",
      " programs: ", ncol(x$weights), "\n",
      " type: ", x$basis_type, "\n",
      " checksum: ", x$audit$checksum, "\n", sep = "")
  invisible(x)
}

#' Build correlation-adjusted program geometry
#'
#' @param basis A \code{proret_basis} object or numeric matrix.
#' @param ridge_fraction Ridge as a fraction of mean Gram diagonal.
#' @param mode \code{"orthogonal_subspace"} or \code{"template"}.
#' @return Geometry transform and diagnostics.
#' @export
build_program_geometry <- function(
    basis,
    ridge_fraction = 1e-3,
    mode = c("orthogonal_subspace", "template")) {
  mode <- match.arg(mode)
  w <- if (inherits(basis, "proret_basis")) basis$weights else basis
  w <- normalize_columns(w)
  gram <- crossprod(w)
  p <- ncol(w)
  if (mode == "template") {
    return(list(mode = mode, transform = diag(p), gram = gram, ridge = 0,
                condition_number = kappa(gram), effective_rank = qr(gram)$rank))
  }
  lambda <- max(0, ridge_fraction) * mean(diag(gram))
  ee <- eigen(gram + diag(lambda, p), symmetric = TRUE)
  floor_ev <- max(ee$values) * 1e-10
  transform <- ee$vectors %*%
    diag(1 / sqrt(pmax(ee$values, floor_ev)), p) %*% t(ee$vectors)
  list(
    mode = mode,
    transform = transform,
    gram = gram,
    ridge = lambda,
    condition_number = max(ee$values) / min(ee$values),
    effective_rank = sum(ee$values)^2 / sum(ee$values^2)
  )
}
