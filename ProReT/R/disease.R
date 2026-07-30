#' Read a local GEO Series Matrix file
#'
#' @param series_matrix Path to a GEO Series Matrix text file, optionally gzip
#'   compressed.
#' @return A list containing an expression matrix, phenotype table, platform
#'   annotation name and source checksum.
#' @export
read_geo_series <- function(series_matrix) {
  require_packages(c("GEOquery", "Biobase"), "reading a GEO Series Matrix")
  series_matrix <- assert_file(series_matrix, "GEO Series Matrix")
  obj <- GEOquery::getGEO(
    filename = series_matrix, GSEMatrix = TRUE, getGPL = FALSE
  )
  eset <- if (methods::is(obj, "ExpressionSet")) obj else obj[[1L]]
  if (!methods::is(eset, "ExpressionSet")) {
    stop("The file did not produce a GEO ExpressionSet.", call. = FALSE)
  }
  expr <- Biobase::exprs(eset)
  pd <- Biobase::pData(eset)
  if (!ncol(expr) || !nrow(expr)) stop("Series Matrix has no expression data.",
                                       call. = FALSE)
  if (!identical(colnames(expr), rownames(pd))) {
    common <- intersect(colnames(expr), rownames(pd))
    if (length(common) < 3L) stop("Cannot align expression and phenotype samples.",
                                  call. = FALSE)
    expr <- expr[, common, drop = FALSE]
    pd <- pd[common, , drop = FALSE]
  }
  list(
    expression = expr,
    phenotype = pd,
    annotation = Biobase::annotation(eset),
    source_file = series_matrix,
    source_sha256 = file_sha256(series_matrix)
  )
}

assign_groups <- function(values, patterns) {
  if (is.null(names(patterns)) || any(!nzchar(names(patterns)))) {
    stop("groups must be a named character vector of regular expressions.",
         call. = FALSE)
  }
  values <- as.character(values)
  match_matrix <- vapply(
    patterns,
    function(p) grepl(p, values, ignore.case = TRUE, perl = TRUE),
    logical(length(values))
  )
  if (is.null(dim(match_matrix))) match_matrix <- matrix(match_matrix, ncol = 1L)
  nmatch <- rowSums(match_matrix)
  if (any(nmatch > 1L)) {
    bad <- head(values[nmatch > 1L], 5L)
    stop("Group patterns overlap for sample labels: ",
         paste(bad, collapse = " | "), call. = FALSE)
  }
  out <- rep(NA_character_, length(values))
  for (j in seq_along(patterns)) out[match_matrix[, j]] <- names(patterns)[j]
  out
}

normalize_expression <- function(expr, method, force_log2 = NA) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  if (any(!is.finite(expr))) {
    stop("Expression matrix contains NA/NaN/Inf.", call. = FALSE)
  }
  method <- match.arg(method, c("microarray_quantile", "log2_only", "none"))
  if (method == "none") return(expr)
  qx <- stats::quantile(expr, c(0, .25, .5, .75, .99, 1), na.rm = TRUE)
  auto_log <- qx[[5L]] > 100 || (qx[[6L]] - qx[[1L]]) > 50
  do_log <- if (is.na(force_log2)) auto_log else isTRUE(force_log2)
  if (do_log) expr <- log2(pmax(expr, 0) + 1)
  if (method == "microarray_quantile") {
    expr <- limma::normalizeBetweenArrays(expr, method = "quantile")
  }
  expr
}

#' Fit a moderated disease signature
#'
#' @param expression Numeric gene-by-sample matrix.
#' @param sample_metadata Data frame with one row per sample.
#' @param sample_column Column containing sample identifiers.
#' @param group_column Column containing assigned group labels.
#' @param target Disease group.
#' @param reference Control group.
#' @param gene_id_type Row-identifier type.
#' @param annotation_package Probe annotation package for PROBEID input.
#' @param paired Whether to fit a paired fixed-effect design.
#' @param block_column Pairing/block column.
#' @param covariates Explicit nuisance covariate columns.
#' @param normalization Expression normalization rule.
#' @param force_log2 Override automatic log2 detection.
#' @return A \code{proret_disease} object.
#' @export
fit_disease_signature <- function(
    expression,
    sample_metadata,
    sample_column = "sample",
    group_column = "group",
    target,
    reference,
    gene_id_type = "ENTREZID",
    annotation_package = NULL,
    paired = FALSE,
    block_column = NULL,
    covariates = character(),
    normalization = "microarray_quantile",
    force_log2 = NA) {
  require_packages("limma", "moderated disease-signature estimation")
  expr <- as.matrix(expression)
  storage.mode(expr) <- "double"
  md <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)
  required <- c(sample_column, group_column, covariates)
  if (isTRUE(paired)) required <- c(required, block_column)
  missing <- setdiff(required, names(md))
  if (length(missing)) stop("Missing metadata columns: ",
                            paste(missing, collapse = ", "), call. = FALSE)
  samples <- intersect(colnames(expr), as.character(md[[sample_column]]))
  if (length(samples) < 4L) stop("Fewer than four aligned samples.", call. = FALSE)
  md <- md[match(samples, md[[sample_column]]), , drop = FALSE]
  expr <- expr[, samples, drop = FALSE]
  groups <- as.character(md[[group_column]])
  keep <- groups %in% c(target, reference)
  md <- md[keep, , drop = FALSE]
  expr <- expr[, keep, drop = FALSE]
  groups <- groups[keep]
  if (sum(groups == target) < 2L || sum(groups == reference) < 2L) {
    stop("Each contrast group requires at least two samples.", call. = FALSE)
  }

  if (isTRUE(paired)) {
    if (is.null(block_column)) stop("block_column is required for paired=TRUE.",
                                    call. = FALSE)
    block <- as.character(md[[block_column]])
    tab <- table(block, groups)
    complete <- rownames(tab)[tab[, target] > 0 & tab[, reference] > 0]
    keep <- !is.na(block) & block %in% complete
    md <- md[keep, , drop = FALSE]
    expr <- expr[, keep, drop = FALSE]
    groups <- groups[keep]
    block <- factor(block[keep])
    if (nlevels(block) < 3L) stop("Fewer than three complete paired blocks.",
                                  call. = FALSE)
  } else {
    block <- NULL
  }

  expr <- normalize_expression(expr, normalization, force_log2)
  ids <- map_ids_to_entrez(rownames(expr), gene_id_type, annotation_package)
  keep_gene <- !is.na(ids) & nzchar(ids) &
    apply(expr, 1L, function(z) all(is.finite(z)))
  if (sum(keep_gene) < 20L) stop("Fewer than 20 genes mapped to Entrez.",
                                 call. = FALSE)
  expr <- limma::avereps(expr[keep_gene, , drop = FALSE], ID = ids[keep_gene])

  safe_reference <- make.names(reference)
  safe_target <- make.names(target)
  if (safe_reference == safe_target) {
    stop("target and reference collapse to the same syntactic name.",
         call. = FALSE)
  }
  group_factor <- factor(groups, levels = c(reference, target))
  design <- stats::model.matrix(~ 0 + group_factor)
  colnames(design) <- c(safe_reference, safe_target)
  if (!is.null(block)) {
    block_mm <- stats::model.matrix(~ block)[, -1L, drop = FALSE]
    design <- cbind(design, block_mm)
  }
  if (length(covariates)) {
    cov_df <- md[, covariates, drop = FALSE]
    if (anyNA(cov_df)) stop("Covariates contain NA.", call. = FALSE)
    cov_mm <- stats::model.matrix(~ ., cov_df)
    cov_mm <- cov_mm[, colnames(cov_mm) != "(Intercept)", drop = FALSE]
    if (ncol(cov_mm)) design <- cbind(design, cov_mm)
  }
  colnames(design) <- make.unique(gsub("[^[:alnum:]_]", "_", colnames(design)))
  if (qr(design)$rank < ncol(design)) {
    stop("Disease design matrix is rank deficient.", call. = FALSE)
  }
  contrast <- limma::makeContrasts(
    contrasts = paste0(safe_target, "-", safe_reference), levels = design
  )
  fit <- limma::lmFit(expr, design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = 1L, number = Inf, sort.by = "P",
                        adjust.method = "BH")
  tt <- data.table::as.data.table(tt, keep.rownames = "entrez_id")
  tt[, score := t]
  out <- list(
    signature = stats::setNames(as.numeric(fit$t[, 1L]), rownames(fit$t)),
    deg = tt,
    expression = expr,
    sample_metadata = md,
    design = design,
    contrast = paste0(target, "_vs_", reference),
    target = target,
    reference = reference
  )
  class(out) <- "proret_disease"
  out
}

#' Prepare a disease signature from a GEO Series Matrix
#'
#' @param series_matrix Local Series Matrix file.
#' @param group_column GEO phenotype column used for grouping.
#' @param groups Named regular-expression vector, for example
#'   \code{c(Control="normal", Disease="tumou?r")}.
#' @param target Disease group name.
#' @param reference Control group name.
#' @param gene_id_type Expression-row identifier type.
#' @param annotation_package Probe annotation package for microarrays.
#' @param paired Whether samples are paired.
#' @param block_column GEO phenotype block column.
#' @param covariates Explicit GEO phenotype covariates.
#' @param normalization Normalization strategy.
#' @param output_dir Optional output directory.
#' @return A \code{proret_disease} object.
#' @export
prepare_disease_from_series <- function(
    series_matrix,
    group_column,
    groups,
    target,
    reference,
    gene_id_type = "PROBEID",
    annotation_package = NULL,
    paired = FALSE,
    block_column = NULL,
    covariates = character(),
    normalization = "microarray_quantile",
    force_log2 = NA,
    output_dir = NULL) {
  geo <- read_geo_series(series_matrix)
  pd <- geo$phenotype
  if (!group_column %in% names(pd)) {
    stop("group_column '", group_column, "' not found. Available columns: ",
         paste(names(pd), collapse = ", "), call. = FALSE)
  }
  assigned <- assign_groups(pd[[group_column]], groups)
  md <- pd
  md$sample <- rownames(pd)
  md$group <- assigned
  n_unmatched <- sum(is.na(md$group))
  if (n_unmatched) {
    message("Excluding ", n_unmatched, " sample(s) unmatched by group patterns.")
  }
  out <- fit_disease_signature(
    geo$expression, md, sample_column = "sample", group_column = "group",
    target = target, reference = reference, gene_id_type = gene_id_type,
    annotation_package = annotation_package, paired = paired,
    block_column = block_column, covariates = covariates,
    normalization = normalization, force_log2 = force_log2
  )
  out$series_matrix <- geo$source_file
  out$series_sha256 <- geo$source_sha256
  out$group_patterns <- groups
  if (!is.null(output_dir)) {
    output_dir <- dir_create(output_dir)
    write_tsv(out$deg, file.path(output_dir, "disease_DEG.tsv"))
    write_tsv(out$sample_metadata, file.path(output_dir, "sample_metadata.tsv"))
    write_tsv(data.frame(sample_index = seq_len(nrow(out$design)), out$design),
              file.path(output_dir, "design_matrix.tsv"))
    saveRDS(out, file.path(output_dir, "disease_signature.rds"))
  }
  out
}

#' @export
print.proret_disease <- function(x, ...) {
  cat("<proret_disease>\n",
      " contrast: ", x$contrast, "\n",
      " genes: ", length(x$signature), "\n",
      " samples: ", nrow(x$sample_metadata), "\n", sep = "")
  invisible(x)
}
