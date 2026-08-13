args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(paste(
    "Usage: Rscript data-raw/build_builtin_program_references.R",
    "K562_FRAMEWORK KOLF_FRAMEWORK"
  ))
}
k562 <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
kolf <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
extdata <- file.path(root, "inst", "extdata")
dir.create(extdata, recursive = TRUE, showWarnings = FALSE)

for (file in list.files(file.path(root, "R"), full.names = TRUE,
                        pattern = "[.]R$")) {
  source(file)
}

write_basis <- function(weights, gene_ids, program_names, destination) {
  stopifnot(nrow(weights) == length(gene_ids),
            ncol(weights) == length(program_names))
  con <- gzfile(destination, open = "wt", compression = 9)
  on.exit(close(con))
  out <- data.frame(
    entrez_id = gene_ids, weights,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  names(out) <- c("entrez_id", program_names)
  utils::write.table(out, con, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE)
}

compress_basis <- function(source, destination) {
  x <- data.table::fread(source)
  con <- gzfile(destination, open = "wt", compression = 9)
  on.exit(close(con))
  utils::write.table(x, con, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE)
}

k562_model_file <- file.path(
  k562, "results_publication_v2_20260723", "cache", "cnmf_program_model.rds"
)
kolf_model_file <- file.path(
  kolf, "results_kolf_publication", "cache", "cnmf_program_model.rds"
)
k562_model <- readRDS(k562_model_file)
kolf_model <- readRDS(kolf_model_file)
write_basis(
  k562_model$weights, k562_model$gene_ids, k562_model$program_names,
  file.path(extdata, "K562_program_50_complete_basis.tsv.gz")
)
write_basis(
  kolf_model$weights, kolf_model$gene_ids, kolf_model$program_names,
  file.path(extdata, "KOLF_program_30_complete_basis.tsv.gz")
)

kolf_common_source <- file.path(
  kolf, "results_kolf_publication", "same_moa_retrieval_v1",
  "frozen_inputs", "KOLF_program_30_common_basis.tsv"
)
compress_basis(
  kolf_common_source,
  file.path(extdata, "KOLF_program_30_common_basis.tsv.gz")
)

wrap_legacy <- function(legacy, basis, core, source_file, label) {
  dt <- data.table::as.data.table(legacy)
  dt[, drug := canonicalize_drug_name(pert_iname)]
  dt[is.na(drug) | !nzchar(drug),
     drug := paste0("pertid_", gsub("[^a-z0-9]", "",
                                    tolower(as.character(pert_id))))]
  programs <- colnames(basis$weights)
  stopifnot(all(programs %in% names(dt)))
  out <- list(
    signatures = dt,
    core_genes = as.character(core),
    program_names = programs,
    basis_checksum = basis$audit$checksum,
    settings = list(
      source = label, gene_space = "BING", pert_type = "trt_cp",
      time_hours = 24, condition_aggregation = "median",
      source_sha256 = file_sha256(source_file)
    ),
    cache_key = hash_object(list(
      basis$audit$checksum, file_sha256(source_file), core
    ))
  )
  out$projection_checksum <- basis_projection_checksum(basis, core)
  out$payload_checksum <- drug_payload_checksum(
    out$signatures, out$core_genes, out$program_names,
    out$projection_checksum
  )
  class(out) <- "proret_drug_programs"
  out
}

kolf_complete_basis <- load_program_basis(
  file.path(extdata, "KOLF_program_30_complete_basis.tsv.gz"),
  "entrez_id", "ENTREZID", "signed_template"
)
kolf_native_file <- file.path(
  kolf, "results_kolf_publication", "cache", "drug_program_signatures.rds"
)
kolf_native_raw <- readRDS(kolf_native_file)
kolf_native <- wrap_legacy(
  kolf_native_raw, kolf_complete_basis,
  attr(kolf_native_raw, "core_genes"), kolf_native_file,
  "KOLF complete basis x GSE70138 LINCS Level 5"
)
saveRDS(
  kolf_native,
  file.path(extdata, "KOLF_LINCS_complete_drug_program_reference.rds"),
  compress = "xz"
)

kolf_common_basis <- load_program_basis(
  file.path(extdata, "KOLF_program_30_common_basis.tsv.gz"),
  "entrez_id", "ENTREZID", "signed_template"
)
kolf_common_file <- file.path(
  root, "data-raw", "KOLF_common_drug_programs.tsv.gz"
)
con <- gzfile(kolf_common_file, open = "rt")
kolf_common_raw <- data.table::as.data.table(utils::read.delim(
  con, check.names = FALSE, stringsAsFactors = FALSE
))
close(con)
kolf_common <- wrap_legacy(
  kolf_common_raw, kolf_common_basis,
  rownames(kolf_common_basis$weights), kolf_common_file,
  "KOLF common 5927-gene basis x GSE70138 LINCS Level 5"
)
saveRDS(
  kolf_common,
  file.path(extdata, "KOLF_LINCS_common_drug_program_reference.rds"),
  compress = "xz"
)
