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

hash_object <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

file_sha256 <- function(path) {
  digest::digest(file = assert_file(path), algo = "sha256", serialize = FALSE)
}

canonicalize_drug_name <- function(x) {
  aliases <- c(
    glibenclamide = "glyburide", salbutamol = "albuterol",
    adrenaline = "epinephrine", noradrenaline = "norepinephrine",
    paracetamol = "acetaminophen", acetylsalicylicacid = "aspirin"
  )
  salts <- paste(c(
    "mesylate", "hydrochloride", "hcl", "sulfate", "sulphate", "citrate",
    "tartrate", "fumarate", "succinate", "acetate", "phosphate", "maleate",
    "besylate", "tosylate", "sodium", "potassium", "calcium", "hydrate",
    "monohydrate", "dihydrate"
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
