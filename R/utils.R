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

hash_text_fields <- function(x) {
  x <- as.character(x)
  missing <- is.na(x)
  x[missing] <- ""
  x <- enc2utf8(x)
  con <- rawConnection(raw(0L), open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.integer(length(x)), con, size = 4L, endian = "big")
  for (i in seq_along(x)) {
    writeBin(as.raw(missing[[i]]), con)
    bytes <- charToRaw(x[[i]])
    writeBin(as.integer(length(bytes)), con, size = 4L, endian = "big")
    if (length(bytes)) writeBin(bytes, con)
  }
  digest::digest(rawConnectionValue(con), algo = "sha256", serialize = FALSE)
}

hash_atomic_column <- function(x) {
  if (is.factor(x) || inherits(x, c("POSIXct", "POSIXlt", "Date"))) {
    x <- as.character(x)
  }
  if (is.character(x)) {
    return(hash_text_fields(c("character", x)))
  }
  if (is.logical(x)) {
    missing <- is.na(x)
    x[missing] <- FALSE
    payload <- c(charToRaw("logical"), as.raw(missing), as.raw(x))
    return(digest::digest(payload, algo = "sha256", serialize = FALSE))
  }
  if (is.integer(x)) {
    missing <- is.na(x)
    x[missing] <- 0L
    con <- rawConnection(raw(0L), open = "wb")
    on.exit(close(con), add = TRUE)
    writeBin(charToRaw("integer"), con)
    writeBin(as.raw(missing), con)
    writeBin(x, con, size = 4L, endian = "big")
    return(digest::digest(rawConnectionValue(con), algo = "sha256",
                          serialize = FALSE))
  }
  if (is.numeric(x)) {
    missing <- is.na(x)
    x[missing] <- 0
    con <- rawConnection(raw(0L), open = "wb")
    on.exit(close(con), add = TRUE)
    writeBin(charToRaw("double"), con)
    writeBin(as.raw(missing), con)
    writeBin(as.double(x), con, size = 8L, endian = "big")
    return(digest::digest(rawConnectionValue(con), algo = "sha256",
                          serialize = FALSE))
  }
  hash_text_fields(c(typeof(x), as.character(x)))
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
  if (!is.null(basis$source_sha256) && nzchar(basis$source_sha256)) {
    return(hash_text_fields(c(
      "projection_basis_file_v1", basis$source_sha256,
      as.character(core_genes), enc2utf8(colnames(basis$weights)),
      basis$basis_type
    )))
  }
  # In-memory bases do not have a source file. Their canonicalized weights are
  # used as a fallback for references created and consumed in the same session.
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
  column_hashes <- vapply(signatures, hash_atomic_column, character(1L))
  hash_text_fields(c(
    "drug_payload_v2", as.character(signature_nrow),
    enc2utf8(signature_names), column_hashes,
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
