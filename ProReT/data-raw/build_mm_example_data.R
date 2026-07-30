args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript data-raw/build_mm_example_data.R FRAMEWORK_DIR")
}

framework <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
package_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
extdata <- file.path(package_root, "inst", "extdata")
dir.create(extdata, recursive = TRUE, showWarnings = FALSE)

for (file in list.files(file.path(package_root, "R"), full.names = TRUE,
                        pattern = "[.]R$")) {
  source(file)
}

basis_source <- file.path(
  framework, "results_publication_v2_20260723", "same_moa_retrieval_v1",
  "frozen_inputs", "K562_program_50_common_basis.tsv"
)
drug_source <- file.path(
  framework, "results_publication_v2_20260723", "cache",
  "drug_program_signatures.rds"
)
disease_source <- file.path(
  framework, "results_publication_v2_20260723", "GSE6477-MM",
  "disease_gene_signatures_moderated_t.tsv"
)
series_source <- file.path(
  framework, "data", "raw", "GSE6477_series_matrix.txt.gz"
)
stopifnot(file.exists(basis_source), file.exists(drug_source),
          file.exists(disease_source), file.exists(series_source))

basis_gz <- file.path(extdata, "K562_program_50_common_basis.tsv.gz")
input <- file(basis_source, open = "rt")
output <- gzfile(basis_gz, open = "wt", compression = 9)
writeLines(readLines(input, warn = FALSE), output)
close(input)
close(output)

disease <- data.table::fread(disease_source)
data.table::setnames(disease, c("entrez_id", "moderated_t"))
data.table::fwrite(
  disease, file.path(extdata, "MM_GSE6477_moderated_t.tsv"), sep = "\t"
)
file.copy(
  series_source,
  file.path(extdata, "GSE6477_series_matrix.txt.gz"),
  overwrite = TRUE
)

basis <- load_program_basis(
  basis_gz, gene_column = "entrez_id",
  gene_id_type = "ENTREZID", basis_type = "signed_template"
)
legacy_raw <- readRDS(drug_source)
core <- as.character(attr(legacy_raw, "core_genes"))
legacy <- data.table::as.data.table(legacy_raw)
legacy[, drug := canonicalize_drug_name(pert_iname)]
legacy[is.na(drug) | !nzchar(drug),
       drug := paste0("pertid_", gsub("[^a-z0-9]", "",
                                      tolower(as.character(pert_id))))]
programs <- colnames(basis$weights)
stopifnot(all(programs %in% names(legacy)))
stopifnot(identical(core, rownames(basis$weights)))

drug_reference <- list(
  signatures = legacy,
  core_genes = core,
  program_names = programs,
  basis_checksum = basis$audit$checksum,
  settings = list(
    source = "GSE70138 LINCS Level 5, precomputed 36-chunk projection",
    gene_space = "BING",
    pert_type = "trt_cp",
    time_hours = 24,
    condition_aggregation = "median",
    basis_kind = "signed_gene_score",
    source_sha256 = file_sha256(drug_source)
  ),
  cache_key = hash_object(list(
    basis$audit$checksum, file_sha256(drug_source), core
  ))
)
class(drug_reference) <- "proret_drug_programs"
saveRDS(
  drug_reference,
  file.path(extdata, "K562_LINCS_drug_program_reference.rds"),
  compress = "xz"
)
