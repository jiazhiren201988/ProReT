#' Manifest for frozen ProReT reference files
#'
#' @return A data frame containing versioned file names, SHA-256 checksums and
#'   release URLs for the frozen program and drug-program references.
#' @export
proret_reference_manifest <- function() {
  base_url <- paste0(
    "https://github.com/jiazhiren201988/ProReT/releases/download/",
    "v0.1.3/"
  )
  x <- data.frame(
    filename = c(
      "K562_program_50_common_basis.tsv.gz",
      "K562_program_50_complete_basis.tsv.gz",
      "K562_LINCS_drug_program_reference.rds",
      "KOLF_program_30_common_basis.tsv.gz",
      "KOLF_program_30_complete_basis.tsv.gz",
      "KOLF_LINCS_common_drug_program_reference.rds",
      "KOLF_LINCS_complete_drug_program_reference.rds"
    ),
    bytes = c(2563164, 3523234, 6042040, 1545134, 6627566, 3646616,
              3683353),
    sha256 = c(
      "a12c9e04212618d5b319cf184d207f39c1065673d5fd301487a88e33302b3d1c",
      "63ccfa3123b390a3098ca8747bb21460950cc64488dfae76c47834348ce29a6c",
      "23253cb12002c41321ba8d90e2791285111f3583e1a32bd76063515d9f1030b9",
      "1dadbbe2154d0475cd410f0909fd554702e93118908b16175a396896b466c921",
      "b6dd7fa60e9547cad5030d812650d4c338a0bbae71ef9ba08a6ed688f46e0aca",
      "ebbb3fc75e5fd29b9077f29e861944a8b070b9a6a48141cd09aa435ae8faed2b",
      "ad9e26f76f300f4f15c37e48f9802aa632b9ccaff3b26351ffcd08d76765d10b"
    ),
    stringsAsFactors = FALSE
  )
  x$url <- paste0(base_url, x$filename)
  x
}

default_reference_cache <- function() {
  configured <- Sys.getenv("PRORET_REFERENCE_CACHE", unset = "")
  if (nzchar(configured)) return(configured)
  if (getRversion() >= "4.0.0") {
    file.path(tools::R_user_dir("ProReT", "cache"), "references", "v0.1.3")
  } else {
    file.path(path.expand("~"), ".cache", "ProReT", "references", "v0.1.3")
  }
}

reference_file_path <- function(filename,
                                cache_dir = default_reference_cache(),
                                download = TRUE) {
  bundled <- system.file("extdata", filename, package = "ProReT")
  if (nzchar(bundled)) return(bundled)
  cached <- file.path(cache_dir, filename)
  if (!file.exists(cached) && isTRUE(download)) {
    download_references(files = filename, cache_dir = cache_dir)
  }
  assert_file(cached, paste0("ProReT reference '", filename, "'"))
}

#' Download frozen ProReT program references
#'
#' Files are obtained from the versioned GitHub release, written atomically
#' and accepted only after file-size and SHA-256 verification.
#'
#' @param files Optional file names from \code{proret_reference_manifest()}.
#'   The default downloads all frozen references.
#' @param cache_dir Versioned local reference cache.
#' @param overwrite Replace verified local files.
#' @return Named normalized paths.
#' @export
download_references <- function(
    files = proret_reference_manifest()$filename,
    cache_dir = default_reference_cache(), overwrite = FALSE) {
  manifest <- proret_reference_manifest()
  unknown <- setdiff(files, manifest$filename)
  if (length(unknown)) {
    stop("Unknown reference file(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  manifest <- manifest[match(files, manifest$filename), , drop = FALSE]
  cache_dir <- dir_create(cache_dir)
  paths <- character(nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    dest <- file.path(cache_dir, manifest$filename[[i]])
    valid <- file.exists(dest) &&
      as.numeric(file.info(dest)$size) == manifest$bytes[[i]] &&
      identical(tolower(file_sha256(dest)), manifest$sha256[[i]])
    if (valid && !isTRUE(overwrite)) {
      paths[[i]] <- dest
      next
    }
    partial <- paste0(dest, ".part")
    if (file.exists(partial)) unlink(partial)
    message("Downloading ", manifest$filename[[i]], " ...")
    if (requireNamespace("curl", quietly = TRUE)) {
      curl::curl_download(manifest$url[[i]], partial, mode = "wb",
                          quiet = FALSE)
    } else {
      utils::download.file(manifest$url[[i]], partial, mode = "wb",
                           quiet = FALSE)
    }
    observed_size <- as.numeric(file.info(partial)$size)
    observed_hash <- file_sha256(partial)
    if (observed_size != manifest$bytes[[i]] ||
        !identical(tolower(observed_hash), manifest$sha256[[i]])) {
      unlink(partial)
      stop("Reference integrity check failed: ", manifest$filename[[i]],
           call. = FALSE)
    }
    if (!file.rename(partial, dest)) {
      unlink(partial)
      stop("Could not finalize reference file: ", dest, call. = FALSE)
    }
    paths[[i]] <- dest
  }
  stats::setNames(
    normalizePath(paths, winslash = "/", mustWork = TRUE),
    manifest$filename
  )
}
