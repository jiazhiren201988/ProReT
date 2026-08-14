`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

assert_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(name, " must be one non-empty character string.", call. = FALSE)
  }
  invisible(x)
}

assert_file <- function(path, label = "file") {
  assert_scalar_character(path, label)
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

require_packages <- function(packages, purpose = "this operation") {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing)) {
    stop("Install required package(s) for ", purpose, ": ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

dir_create <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create directory: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

pick_col <- function(x, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit)) return(hit[[1L]])
  if (required) {
    stop("Cannot find ", label, ". Tried: ", paste(candidates, collapse = ", "),
         ". Available: ", paste(names(x), collapse = ", "), call. = FALSE)
  }
  NULL
}

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

zscore_vector <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (sum(ok) < 3L) stop("At least three finite values are required.", call. = FALSE)
  s <- stats::sd(x[ok])
  if (!is.finite(s) || s == 0) stop("Signature has zero finite variance.", call. = FALSE)
  out <- rep(NA_real_, length(x))
  out[ok] <- (x[ok] - mean(x[ok])) / s
  out
}

normalize_columns <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  norms <- sqrt(colSums(x^2))
  if (any(!is.finite(norms) | norms <= 0)) {
    stop("Every program must have a finite, non-zero L2 norm.", call. = FALSE)
  }
  sweep(x, 2L, norms, "/")
}

hash_object <- function(x) {
  # Serialization version 2 uses the portable XDR representation and avoids
  # checksum drift between supported R releases and operating systems.
  digest::digest(
    x, algo = "sha256", serialize = TRUE, serializeVersion = 2L
  )
}

hash_numeric_matrix <- function(x, digits = 12L) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  canonical_values <- formatC(
    as.vector(x), digits = digits, format = "e", decimal.mark = "."
  )
  hash_object(list(
    as.integer(dim(x)), enc2utf8(rownames(x)), enc2utf8(colnames(x)),
    canonical_values
  ))
}

file_sha256 <- function(path) {
  digest::digest(file = assert_file(path), algo = "sha256", serialize = FALSE)
}

file_sha512 <- function(path) {
  digest::digest(file = assert_file(path), algo = "sha512", serialize = FALSE)
}

basis_projection_checksum <- function(basis, core_genes) {
  if (!inherits(basis, "proret_basis")) {
    stop("basis must be a proret_basis object.", call. = FALSE)
  }
  core_genes <- as.character(core_genes)
  if (anyDuplicated(core_genes) || !all(core_genes %in% rownames(basis$weights))) {
    stop("core_genes must be unique and present in the basis.", call. = FALSE)
  }
  # Bind the reference to the source weights and coordinate order. Hashing
  # normalized weights would incorporate platform-specific BLAS rounding even
  # though those differences are numerically immaterial to the projection.
  w <- basis$weights[core_genes, , drop = FALSE]
  hash_numeric_matrix(w, digits = 12L)
}

drug_payload_checksum <- function(signatures, core_genes, program_names,
                                  projection_checksum) {
  if (!is.list(signatures)) {
    stop("signatures must be a data-frame-like list.", call. = FALSE)
  }
  signature_names <- names(signatures)
  signature_nrow <- if (length(signatures)) length(signatures[[1L]]) else 0L
  if (length(signatures) &&
      any(vapply(signatures, length, integer(1L)) != signature_nrow)) {
    stop("All signature columns must have equal length.", call. = FALSE)
  }
  canonical_columns <- lapply(seq_along(signatures), function(i) {
    column <- signatures[[i]]
    if (is.factor(column)) column <- as.character(column)
    if (inherits(column, c("POSIXct", "POSIXlt", "Date"))) {
      column <- as.character(column)
    }
    if (is.character(column)) column <- enc2utf8(column)
    attributes(column) <- NULL
    column
  })
  hash_object(list(
    as.integer(signature_nrow), enc2utf8(signature_names),
    vapply(canonical_columns, hash_object, character(1L)),
    as.character(core_genes), as.character(program_names),
    as.character(projection_checksum)
  ))
}

canonicalize_drug_name <- function(x) {
  aliases <- c(
    glibenclamide = "glyburide", salbutamol = "albuterol",
    adrenaline = "epinephrine", noradrenaline = "norepinephrine",
    noradrenalin = "norepinephrine",
    paracetamol = "acetaminophen", acetylsalicylicacid = "aspirin"
  )
  salts <- paste(c(
    "mesylate", "hydrochloride", "hcl", "sulfate", "sulphate", "citrate",
    "tartrate", "fumarate", "succinate", "acetate", "phosphate", "nitrate",
    "maleate", "besylate", "tosylate", "esylate", "lactate", "oxalate",
    "malonate", "adipate", "stearate", "gluconate", "hydrobromide",
    "dihydrochloride", "diacetate", "dimesylate", "ditosylate", "camsylate",
    "napadisylate", "sodium", "potassium", "calcium", "magnesium", "zinc",
    "aluminum", "tromethamine", "hydrate", "monohydrate", "dihydrate",
    "trihydrate", "sesquihydrate", "ethanolamide", "benzathine", "procaine",
    "depot", "delayedrelease", "extendedrelease"
  ), collapse = "|")
  y <- tolower(as.character(x))
  y <- gsub(paste0("\\b(", salts, ")\\b"), " ", y, perl = TRUE)
  y <- gsub("[^a-z0-9]", "", y)
  replace <- match(y, names(aliases))
  y[!is.na(replace)] <- unname(aliases[replace[!is.na(replace)]])
  y
}

write_tsv <- function(x, path) {
  dir_create(dirname(path))
  data.table::fwrite(data.table::as.data.table(x), path, sep = "\t", na = "NA")
  invisible(path)
}
