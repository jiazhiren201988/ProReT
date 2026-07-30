#' LINCS reference manifest
#'
#' @return A data frame containing the frozen LINCS Level 5 files used by
#'   ProReT.
#' @export
lincs_manifest <- function() {
  data.frame(
    key = c("gctx", "sig_info", "gene_info"),
    filename = c(
      "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx",
      "GSE70138_Broad_LINCS_sig_info.txt",
      "GSE92742_Broad_LINCS_gene_info.txt"
    ),
    url = c(
      paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/",
             "GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"),
      paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/",
             "GSE70138/suppl/GSE70138_Broad_LINCS_sig_info.txt"),
      paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/",
             "GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt")
    ),
    expected_bytes = c(5824120558, 18260895, 621887),
    stringsAsFactors = FALSE
  )
}

default_lincs_cache <- function() {
  if (getRversion() >= "4.0.0") {
    tools::R_user_dir("ProReT", "cache")
  } else {
    file.path(path.expand("~"), ".cache", "ProReT")
  }
}

#' Check a local LINCS reference
#'
#' @param cache_dir Directory containing the three manifest files.
#' @param strict_size Require the frozen file sizes to match.
#' @return Named normalized paths.
#' @export
check_lincs_files <- function(cache_dir = default_lincs_cache(),
                              strict_size = TRUE) {
  manifest <- lincs_manifest()
  paths <- file.path(cache_dir, manifest$filename)
  missing <- !file.exists(paths)
  if (any(missing)) {
    stop("Missing LINCS file(s): ",
         paste(manifest$filename[missing], collapse = ", "),
         ". Run download_lincs().", call. = FALSE)
  }
  if (isTRUE(strict_size)) {
    observed <- as.numeric(file.info(paths)$size)
    bad <- observed != manifest$expected_bytes
    if (any(bad)) {
      stop("LINCS size check failed for: ",
           paste(manifest$filename[bad], collapse = ", "),
           ". Remove incomplete files and download again.", call. = FALSE)
    }
  }
  stats::setNames(
    normalizePath(paths, winslash = "/", mustWork = TRUE),
    manifest$key
  )
}

#' Download and cache the frozen LINCS reference
#'
#' The approximately 5.6 GB GCTX matrix is stored in a user cache and is never
#' bundled inside the installed R package.
#'
#' @param cache_dir Cache directory.
#' @param overwrite Replace existing files.
#' @param verify Verify expected file sizes after download.
#' @return Named paths returned by \code{check_lincs_files()}.
#' @export
download_lincs <- function(cache_dir = default_lincs_cache(),
                           overwrite = FALSE,
                           verify = TRUE) {
  cache_dir <- dir_create(cache_dir)
  manifest <- lincs_manifest()
  for (i in seq_len(nrow(manifest))) {
    dest <- file.path(cache_dir, manifest$filename[[i]])
    if (file.exists(dest) && !isTRUE(overwrite)) next
    partial <- paste0(dest, ".part")
    if (file.exists(partial)) unlink(partial)
    message("Downloading ", manifest$filename[[i]], " ...")
    if (requireNamespace("curl", quietly = TRUE)) {
      curl::curl_download(manifest$url[[i]], partial, quiet = FALSE,
                          mode = "wb")
    } else {
      utils::download.file(manifest$url[[i]], partial, mode = "wb",
                           quiet = FALSE)
    }
    if (!file.rename(partial, dest)) {
      stop("Could not finalize download: ", dest, call. = FALSE)
    }
  }
  check_lincs_files(cache_dir, strict_size = verify)
}

standardize_lincs_metadata <- function(metadata) {
  md <- data.table::as.data.table(metadata)
  columns <- list(
    sig_id = c("sig_id", "signature_id"),
    pert_id = c("pert_id", "compound_id"),
    pert_iname = c("pert_iname", "compound_name"),
    cell_iname = c("cell_iname", "cell_id"),
    pert_type = c("pert_type"),
    pert_itime = c("pert_itime", "pert_time"),
    pert_idose = c("pert_idose", "pert_dose"),
    is_hiq = c("is_hiq", "is_hq", "high_quality")
  )
  out <- list()
  for (nm in names(columns)) {
    col <- pick_col(md, columns[[nm]], required = nm %in%
                      c("sig_id", "pert_id", "pert_iname", "cell_iname"),
                    label = paste0("LINCS ", nm))
    out[[nm]] <- if (is.null(col)) rep(NA_character_, nrow(md)) else md[[col]]
  }
  data.table::as.data.table(out)
}

parse_lincs_hours <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(x))))
}

filter_lincs_metadata <- function(metadata, pert_type = "trt_cp",
                                  time_hours = 24,
                                  cell_lines = NULL,
                                  require_hq = FALSE) {
  md <- standardize_lincs_metadata(metadata)
  wanted_pert_type <- pert_type
  wanted_hours <- time_hours
  if (all(is.na(md$pert_type))) {
    if (!is.null(pert_type)) {
      stop("LINCS pert_type is required but absent.", call. = FALSE)
    }
  } else if (!is.null(wanted_pert_type)) {
    md <- md[as.character(pert_type) %in% wanted_pert_type]
  }
  if (all(is.na(md$pert_itime))) {
    if (!is.null(time_hours)) {
      stop("LINCS pert_itime is required but absent.", call. = FALSE)
    }
  } else if (!is.null(wanted_hours)) {
    md <- md[parse_lincs_hours(pert_itime) %in% as.numeric(wanted_hours)]
  }
  if (!is.null(cell_lines)) md <- md[cell_iname %in% cell_lines]
  if (isTRUE(require_hq)) {
    if (all(is.na(md$is_hiq))) {
      stop("HQ filtering was requested, but no quality field exists.",
           call. = FALSE)
    }
    md <- md[as_flag(is_hiq)]
  }
  if (!nrow(md)) stop("No LINCS signatures remain after filtering.",
                       call. = FALSE)
  md
}

read_lincs_gene_space <- function(gene_info_file, gene_space = "BING") {
  gi <- data.table::fread(assert_file(gene_info_file, "LINCS gene info"))
  id_col <- pick_col(gi, c("pr_gene_id", "gene_id"), label = "LINCS gene ID")
  ids <- as.character(gi[[id_col]])
  if (toupper(gene_space) == "ALL") return(ids)
  if (toupper(gene_space) != "BING") {
    stop("gene_space must be 'BING' or 'ALL'.", call. = FALSE)
  }
  flag <- pick_col(gi, c("pr_is_bing", "is_bing"),
                   label = "LINCS BING flag")
  ids[as_flag(gi[[flag]])]
}

aggregate_program_instances <- function(scores, metadata) {
  programs <- rownames(scores)
  dt <- data.table::as.data.table(t(scores))
  data.table::setnames(dt, programs)
  dt <- cbind(metadata, dt)
  dt[, drug := canonicalize_drug_name(pert_iname)]
  dt[is.na(drug) | !nzchar(drug),
     drug := paste0("pertid_", gsub("[^a-z0-9]", "",
                                    tolower(as.character(pert_id))))]
  dt[, c(list(
    pert_iname = as.character(pert_iname[[1L]]),
    n_instances = .N,
    n_doses = data.table::uniqueN(pert_idose[!is.na(pert_idose)])
  ), lapply(.SD, stats::median, na.rm = TRUE)),
  by = .(pert_id, drug, cell_iname), .SDcols = programs]
}

make_drug_program_object <- function(signatures, core_genes, basis,
                                     settings, cache_key = NULL) {
  out <- list(
    signatures = signatures,
    core_genes = core_genes,
    program_names = colnames(basis$weights),
    basis_checksum = basis$audit$checksum,
    settings = settings,
    cache_key = cache_key
  )
  class(out) <- "proret_drug_programs"
  out
}

#' Project an in-memory LINCS-like matrix
#'
#' This function is useful for unit tests and smaller custom perturbation
#' references. The expression matrix must be genes by signatures.
#'
#' @param expression Gene-by-signature numeric matrix.
#' @param metadata Signature metadata.
#' @param basis A \code{proret_basis}.
#' @param gene_space_ids Optional permitted gene universe.
#' @param pert_type,time_hours,cell_lines,require_hq LINCS filters.
#' @return A \code{proret_drug_programs} object.
#' @export
project_lincs_matrix <- function(
    expression, metadata, basis, gene_space_ids = NULL,
    pert_type = "trt_cp", time_hours = 24, cell_lines = NULL,
    require_hq = FALSE) {
  if (!inherits(basis, "proret_basis")) {
    stop("basis must be a proret_basis object.", call. = FALSE)
  }
  x <- as.matrix(expression)
  storage.mode(x) <- "double"
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("expression needs gene row names and signature column names.",
         call. = FALSE)
  }
  if (any(!is.finite(x))) stop("LINCS expression contains NA/NaN/Inf.",
                               call. = FALSE)
  md <- filter_lincs_metadata(metadata, pert_type, time_hours, cell_lines,
                              require_hq)
  md <- md[sig_id %in% colnames(x)]
  if (!nrow(md)) stop("No metadata signatures occur in expression.",
                       call. = FALSE)
  x <- x[, md$sig_id, drop = FALSE]
  allowed <- gene_space_ids %||% rownames(x)
  core <- Reduce(intersect, list(rownames(basis$weights), rownames(x),
                                 as.character(allowed)))
  if (length(core) < 20L) stop("Fewer than 20 shared basis/LINCS genes.",
                               call. = FALSE)
  w <- normalize_columns(basis$weights[core, , drop = FALSE])
  scores <- crossprod(w, x[core, , drop = FALSE])
  signatures <- aggregate_program_instances(scores, md)
  settings <- list(pert_type = pert_type, time_hours = time_hours,
                   cell_lines = cell_lines, require_hq = require_hq,
                   gene_space = if (is.null(gene_space_ids)) "input" else "custom")
  make_drug_program_object(signatures, core, basis, settings)
}

#' Project LINCS Level 5 GCTX signatures in chunks
#'
#' @param gctx_file,sig_info_file,gene_info_file Frozen LINCS files.
#' @param basis A \code{proret_basis}.
#' @param gene_space \code{"BING"} (default) or \code{"ALL"}.
#' @param pert_type,time_hours,cell_lines,require_hq LINCS filters.
#' @param chunk_size Number of signatures parsed per GCTX block.
#' @param cache_dir Projection-cache directory.
#' @param overwrite Recompute an existing projection cache.
#' @return A \code{proret_drug_programs} object.
#' @export
project_lincs_gctx <- function(
    gctx_file, sig_info_file, gene_info_file, basis,
    gene_space = "BING", pert_type = "trt_cp", time_hours = 24,
    cell_lines = NULL, require_hq = FALSE, chunk_size = 3000L,
    cache_dir = file.path(default_lincs_cache(), "projections"),
    overwrite = FALSE) {
  if (!requireNamespace("cmapR", quietly = TRUE)) {
    stop("Package 'cmapR' is required to read LINCS GCTX files. ",
         "Install it with BiocManager::install('cmapR').", call. = FALSE)
  }
  if (!inherits(basis, "proret_basis")) {
    stop("basis must be a proret_basis object.", call. = FALSE)
  }
  gctx_file <- assert_file(gctx_file, "LINCS GCTX")
  sig_info_file <- assert_file(sig_info_file, "LINCS signature info")
  gene_info_file <- assert_file(gene_info_file, "LINCS gene info")
  md <- data.table::fread(sig_info_file)
  md <- filter_lincs_metadata(md, pert_type, time_hours, cell_lines, require_hq)
  allowed <- read_lincs_gene_space(gene_info_file, gene_space)
  rid <- as.character(cmapR::read_gctx_ids(gctx_file, dim = "row"))
  cid <- as.character(cmapR::read_gctx_ids(gctx_file, dim = "col"))
  md <- md[sig_id %in% cid]
  if (!nrow(md)) stop("No filtered signatures are present in the GCTX.",
                       call. = FALSE)
  core <- Reduce(intersect, list(rownames(basis$weights), rid, allowed))
  if (length(core) < 20L) stop("Fewer than 20 shared basis/LINCS genes.",
                               call. = FALSE)
  settings <- list(
    basis_checksum = basis$audit$checksum,
    gctx_sha256 = file_sha256(gctx_file),
    sig_info_sha256 = file_sha256(sig_info_file),
    gene_info_sha256 = file_sha256(gene_info_file),
    gene_space = toupper(gene_space), pert_type = pert_type,
    time_hours = time_hours, cell_lines = sort(cell_lines),
    require_hq = require_hq, core_genes = core
  )
  cache_key <- hash_object(settings)
  cache_dir <- dir_create(cache_dir)
  cache_file <- file.path(cache_dir, paste0("drug_programs_", cache_key, ".rds"))
  if (file.exists(cache_file) && !isTRUE(overwrite)) return(readRDS(cache_file))

  w <- normalize_columns(basis$weights[core, , drop = FALSE])
  selected <- intersect(cid, md$sig_id)
  blocks <- split(selected, ceiling(seq_along(selected) / as.integer(chunk_size)))
  projected <- vector("list", length(blocks))
  for (i in seq_along(blocks)) {
    message("Projecting LINCS block ", i, "/", length(blocks))
    parsed <- cmapR::parse_gctx(gctx_file, rid = core, cid = blocks[[i]],
                               matrix_only = TRUE)
    g <- parsed@mat
    g <- g[core, blocks[[i]], drop = FALSE]
    projected[[i]] <- crossprod(w, g)
  }
  scores <- do.call(cbind, projected)
  md <- md[match(colnames(scores), sig_id)]
  signatures <- aggregate_program_instances(scores, md)
  out <- make_drug_program_object(signatures, core, basis, settings, cache_key)
  saveRDS(out, cache_file)
  out
}

#' @export
print.proret_drug_programs <- function(x, ...) {
  cat("<proret_drug_programs>\n",
      " drug-cell profiles: ", nrow(x$signatures), "\n",
      " core genes: ", length(x$core_genes), "\n",
      " programs: ", length(x$program_names), "\n", sep = "")
  invisible(x)
}
