# Add content-derived integrity fields to the frozen drug-program references.

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run this script from the ProReT repository root.")
}

source_files <- list.files(file.path(root, "R"), pattern = "[.]R$",
                           full.names = TRUE)
for (file in source_files) source(file)

reference_dir <- commandArgs(trailingOnly = TRUE)
if (!length(reference_dir)) {
  reference_dir <- file.path(root, "inst", "extdata")
} else {
  reference_dir <- normalizePath(reference_dir[[1L]], winslash = "/",
                                 mustWork = TRUE)
}
specification <- list(
  list(
    drug = "K562_LINCS_drug_program_reference.rds",
    basis = "K562_program_50_common_basis.tsv.gz"
  ),
  list(
    drug = "KOLF_LINCS_common_drug_program_reference.rds",
    basis = "KOLF_program_30_common_basis.tsv.gz"
  ),
  list(
    drug = "KOLF_LINCS_complete_drug_program_reference.rds",
    basis = "KOLF_program_30_complete_basis.tsv.gz"
  )
)

for (item in specification) {
  if (!all(file.exists(file.path(reference_dir,
                                 c(item$basis, item$drug))))) {
    next
  }
  basis <- load_program_basis(
    file.path(reference_dir, item$basis), gene_column = "entrez_id",
    gene_id_type = "ENTREZID", basis_type = "signed_template"
  )
  path <- file.path(reference_dir, item$drug)
  object <- readRDS(path)
  object$basis_checksum <- basis$audit$checksum
  object$projection_checksum <- basis_projection_checksum(
    basis, object$core_genes
  )
  object$payload_checksum <- drug_payload_checksum(
    object$signatures, object$core_genes, object$program_names,
    object$projection_checksum
  )
  saveRDS(object, path, version = 3)
  verified <- readRDS(path)
  verified_checksum <- drug_payload_checksum(
    verified$signatures, verified$core_genes, verified$program_names,
    verified$projection_checksum
  )
  if (!identical(verified$payload_checksum, verified_checksum)) {
    stop("Payload checksum did not survive serialization: ", path)
  }
  message("Updated integrity metadata: ", path)
}
