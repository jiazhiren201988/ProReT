#' Project a disease signature into program space
#'
#' Observed moderated t statistics are standardized across observed core genes;
#' unmeasured core genes are assigned zero only after standardization.
#'
#' @param disease A \code{proret_disease} object or named numeric signature.
#' @param basis A \code{proret_basis}.
#' @param core_genes Frozen genes shared by the basis and drug reference.
#' @param min_coverage Minimum fraction of core genes observed in the disease.
#' @param min_observed Minimum observed genes.
#' @return A \code{proret_disease_program} object.
#' @export
project_disease_signature <- function(disease, basis, core_genes,
                                      min_coverage = 0.80,
                                      min_observed = 20L) {
  if (!inherits(basis, "proret_basis")) {
    stop("basis must be a proret_basis object.", call. = FALSE)
  }
  sig <- if (inherits(disease, "proret_disease")) disease$signature else disease
  if (is.null(names(sig)) || any(!nzchar(names(sig)))) {
    stop("Disease signature must be a named numeric vector.", call. = FALSE)
  }
  if (any(!is.finite(sig))) stop("Disease signature contains NA/NaN/Inf.",
                                 call. = FALSE)
  core <- intersect(as.character(core_genes), rownames(basis$weights))
  if (length(core) != length(unique(core_genes))) {
    stop("core_genes are not all represented in the selected basis.",
         call. = FALSE)
  }
  observed <- intersect(core, names(sig))
  coverage <- length(observed) / length(core)
  if (length(observed) < min_observed) {
    stop("Too few observed core genes: ", length(observed), call. = FALSE)
  }
  if (coverage < min_coverage) {
    stop(sprintf("Disease core-gene coverage %.1f%% is below %.1f%%.",
                 100 * coverage, 100 * min_coverage), call. = FALSE)
  }
  z <- stats::setNames(rep(0, length(core)), core)
  z[observed] <- zscore_vector(sig[observed])
  w <- normalize_columns(basis$weights[core, , drop = FALSE])
  ps <- as.numeric(crossprod(w, z))
  names(ps) <- colnames(w)
  out <- list(
    program_signature = ps, standardized_gene_signature = z,
    core_genes = core, observed_genes = observed, coverage = coverage,
    W_core = w, basis_checksum = basis$audit$checksum,
    program_names = colnames(w)
  )
  class(out) <- "proret_disease_program"
  out
}

negative_cosine <- function(a, b) {
  den <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (!is.finite(den) || den <= 0) return(NA_real_)
  -sum(a * b) / den
}

#' Score disease-drug reversal in program space
#'
#' @param disease_projection A \code{proret_disease_program}.
#' @param drug_programs A \code{proret_drug_programs}.
#' @param geometry \code{"orthogonal_subspace"} (default) or
#'   \code{"template"}.
#' @param ridge_fraction Ridge fraction used for Gram whitening.
#' @return Drug-cell reversal scores.
#' @export
score_reversal <- function(
    disease_projection, drug_programs,
    geometry = c("orthogonal_subspace", "template"),
    ridge_fraction = 1e-3) {
  geometry <- match.arg(geometry)
  if (!inherits(disease_projection, "proret_disease_program")) {
    stop("disease_projection has the wrong class.", call. = FALSE)
  }
  if (!inherits(drug_programs, "proret_drug_programs")) {
    stop("drug_programs has the wrong class.", call. = FALSE)
  }
  if (!identical(disease_projection$basis_checksum,
                 drug_programs$basis_checksum)) {
    stop("Disease and drug projections use different basis checksums.",
         call. = FALSE)
  }
  if (!identical(disease_projection$core_genes, drug_programs$core_genes)) {
    stop("Disease and drug projections use different core genes.",
         call. = FALSE)
  }
  programs <- disease_projection$program_names
  if (!identical(programs, drug_programs$program_names)) {
    stop("Disease and drug program coordinates do not match.", call. = FALSE)
  }
  geo <- build_program_geometry(disease_projection$W_core,
                                ridge_fraction, geometry)
  d <- as.numeric(disease_projection$program_signature %*% geo$transform)
  dt <- data.table::copy(drug_programs$signatures)
  missing <- setdiff(programs, names(dt))
  if (length(missing)) stop("Drug program columns are missing.", call. = FALSE)
  drug_coords <- as.matrix(dt[, programs, with = FALSE]) %*% geo$transform
  reversal <- apply(drug_coords, 1L, negative_cosine, b = d)
  dt[, reversal_score := reversal]
  attr(dt, "geometry") <- geo
  dt
}

#' Aggregate drug-cell reversal scores
#'
#' @param cell_scores Output from \code{score_reversal()}.
#' @param min_cell_lines Minimum eligible cell lines per drug.
#' @param cell_panel Optional fixed cell panel. If absent, all observed cell
#'   lines are eligible.
#' @param aggregator \code{"mean"} or \code{"median"}.
#' @return Compound-level ranking table.
#' @export
aggregate_drug_scores <- function(cell_scores, min_cell_lines = 3L,
                                  cell_panel = NULL,
                                  aggregator = c("mean", "median"),
                                  cell_panel_policy = c("coverage", "all",
                                                        "configured"),
                                  min_cell_panel_coverage = 0.80,
                                  min_drug_panel_fraction = 0.60) {
  aggregator <- match.arg(aggregator)
  cell_panel_policy <- match.arg(cell_panel_policy)
  dt <- data.table::as.data.table(cell_scores)
  required <- c("drug", "cell_iname", "reversal_score")
  missing <- setdiff(required, names(dt))
  if (length(missing)) stop("Missing score columns: ",
                            paste(missing, collapse = ", "), call. = FALSE)
  dt[is.na(drug) | !nzchar(drug),
     drug := paste0("pertid_", gsub("[^a-z0-9]", "",
                                    tolower(as.character(pert_id))))]
  display <- dt[, .(
    pert_id = paste(sort(unique(as.character(pert_id))), collapse = ","),
    pert_iname = names(sort(table(pert_iname), decreasing = TRUE))[[1L]]
  ), by = drug]
  count_col <- intersect(c("n_signatures", "n_instances"), names(dt))
  if (length(count_col)) {
    count_col <- count_col[[1L]]
    collapsed <- dt[, .(
      reversal_score = stats::median(reversal_score, na.rm = TRUE),
      total_signatures = sum(get(count_col), na.rm = TRUE)
    ), by = .(drug, cell_iname)]
  } else {
    collapsed <- dt[, .(
      reversal_score = stats::median(reversal_score, na.rm = TRUE),
      total_signatures = .N
    ), by = .(drug, cell_iname)]
  }
  dt <- merge(collapsed, display, by = "drug", all.x = TRUE, sort = FALSE)
  dc_unique <- unique(dt[, .(drug, cell_iname)])
  n_drugs <- data.table::uniqueN(dc_unique$drug)
  cell_coverage <- dc_unique[, .(n_drugs = .N), by = cell_iname]
  # Use a differently named scalar to avoid data.table name masking.
  coverage_denominator <- n_drugs
  cell_coverage[, coverage_fraction := n_drugs / coverage_denominator]
  if (!is.null(cell_panel)) cell_panel_policy <- "configured"
  panel <- switch(
    cell_panel_policy,
    coverage = cell_coverage[
      coverage_fraction >= min_cell_panel_coverage, cell_iname
    ],
    all = cell_coverage$cell_iname,
    configured = as.character(cell_panel %||% character())
  )
  panel <- intersect(unique(panel), unique(dt$cell_iname))
  if (length(panel) < min_cell_lines) {
    stop("The selected common cell panel has fewer than min_cell_lines.",
         call. = FALSE)
  }
  cell_coverage[, selected_for_panel := cell_iname %in% panel]
  dt <- dt[cell_iname %in% panel]
  fun <- if (aggregator == "mean") mean else stats::median
  out <- dt[is.finite(reversal_score), .(
    reversal_score = fun(reversal_score),
    n_cell_lines = data.table::uniqueN(cell_iname),
    panel_size = length(panel),
    panel_fraction = data.table::uniqueN(cell_iname) / length(panel),
    total_signatures = sum(total_signatures, na.rm = TRUE),
    cell_score_sd = if (.N > 1L) stats::sd(reversal_score) else NA_real_,
    pert_id = as.character(pert_id[[1L]]),
    pert_iname = as.character(pert_iname[[1L]])
  ), by = drug]
  required_fraction <- if (cell_panel_policy == "all") 0 else
    min_drug_panel_fraction
  out <- out[
    n_cell_lines >= as.integer(min_cell_lines) &
      panel_fraction >= required_fraction
  ]
  data.table::setorder(out, -reversal_score, drug)
  out[, rank := seq_len(.N)]
  attr(out, "cell_panel") <- sort(panel)
  attr(out, "cell_coverage") <- cell_coverage[order(-coverage_fraction)]
  out[]
}

#' Rank drugs from disease and drug program projections
#'
#' @param disease_projection A projected disease.
#' @param drug_programs Projected LINCS reference.
#' @param geometry,ridge_fraction,min_cell_lines,cell_panel,aggregator Scoring
#'   and aggregation options.
#' @return A \code{proret_result}.
#' @export
rank_drugs <- function(
    disease_projection, drug_programs,
    geometry = "orthogonal_subspace", ridge_fraction = 1e-3,
    min_cell_lines = 3L, cell_panel = NULL, aggregator = "mean",
    cell_panel_policy = "coverage", min_cell_panel_coverage = 0.80,
    min_drug_panel_fraction = 0.60) {
  cell_scores <- score_reversal(disease_projection, drug_programs,
                                geometry, ridge_fraction)
  ranking <- aggregate_drug_scores(
    cell_scores, min_cell_lines, cell_panel, aggregator,
    cell_panel_policy, min_cell_panel_coverage, min_drug_panel_fraction
  )
  out <- list(
    ranking = ranking, cell_scores = cell_scores,
    disease_projection = disease_projection,
    settings = list(geometry = geometry, ridge_fraction = ridge_fraction,
                    min_cell_lines = min_cell_lines,
                    cell_panel = attr(ranking, "cell_panel"),
                    aggregator = aggregator,
                    cell_panel_policy = cell_panel_policy,
                    min_cell_panel_coverage = min_cell_panel_coverage,
                    min_drug_panel_fraction = min_drug_panel_fraction)
  )
  class(out) <- "proret_result"
  out
}

#' @export
print.proret_result <- function(x, ...) {
  cat("<proret_result>\n",
      " ranked drugs: ", nrow(x$ranking), "\n",
      " core-gene coverage: ",
      sprintf("%.1f%%", 100 * x$disease_projection$coverage), "\n", sep = "")
  invisible(x)
}
