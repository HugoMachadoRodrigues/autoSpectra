# R/ossl_data.R — Download and manage OSSL v1.2 spectral data
#
# Data source: Open Soil Spectral Library (OSSL)
# https://docs.soilspectroscopy.org
#
# Files (Google Cloud Storage, public):
#   VisNIR L0: ossl_visnir_L0_v1.2.csv.gz   (350-2500 nm, reflectance)
#   MIR   L0: ossl_mir_L0_v1.2.csv.gz        (600-4000 cm-1, absorbance)
#   Soillab L1: ossl_soillab_L1_v1.2.csv.gz  (harmonized soil properties)
#
# Joining key: id.layer_uuid_txt

.OSSL_BASE <- "https://storage.googleapis.com/soilspec4gg-public/"

.ossl_files <- list(
  visnir  = "ossl_visnir_L0_v1.2.csv.gz",
  mir     = "ossl_mir_L0_v1.2.csv.gz",
  soillab = "ossl_soillab_L1_v1.2.csv.gz",
  # NeoSpectra handheld library (Sanderman et al. 2023), same bucket
  neospectra_nir     = "neospectra_nir_v1.2.csv.gz",
  neospectra_soillab = "neospectra_soillab_v1.2.csv.gz",
  neospectra_soilsite = "neospectra_soilsite_v1.2.csv.gz"
)

#' Default local cache directory for OSSL data
#' @export
ossl_cache_dir <- function() {
  d <- file.path(tools::R_user_dir("autoSpectra", "data"), "ossl_v1.2")
  dir_create(d)
  d
}

#' Check that a gzipped OSSL cache file is complete
#'
#' An interrupted download leaves a truncated \code{.csv.gz} behind. Because the
#' old cache check was a bare \code{file.exists()}, that stub would then be
#' treated as a valid cache forever, and \code{fread()} would either fail or --
#' worse -- silently read a partial database.
#'
#' Note that R's own \code{gzfile()} is NOT sufficient here: it reads a truncated
#' member up to the cut and then reports a clean EOF, with no error or warning.
#' We therefore verify the CRC32/ISIZE trailer with \code{gzip -t}, which is a
#' real gzip implementation. Where \code{gzip} is unavailable (some Windows
#' setups) we fall back to a size sanity check, and rely on the atomic
#' \code{.part} -> rename in \code{ossl_download_file()} for correctness.
#'
#' @param path Path to a \code{.csv.gz} file
#' @return \code{TRUE} if the file exists and its gzip stream is complete
#' @keywords internal
.ossl_gz_ok <- function(path) {
  if (!file.exists(path) || file.size(path) < 1024) return(FALSE)

  gz <- Sys.which("gzip")
  if (nzchar(gz)) {
    status <- suppressWarnings(
      system2(gz, args = c("-t", shQuote(path)),
              stdout = FALSE, stderr = FALSE)
    )
    return(identical(as.integer(status), 0L))
  }

  # No gzip binary: we cannot verify the trailer. Accept the file, since
  # ossl_download_file() only ever publishes fully-transferred downloads.
  TRUE
}

#' Download one OSSL component file to the cache directory
#'
#' Downloads from the official OSSL distribution bucket documented at
#' \url{https://docs.soilspectroscopy.org/db-access.html}. autoSpectra never
#' bundles or re-hosts OSSL data — each user fetches it from the official
#' source into their own local cache.
#'
#' The file is written to a temporary path and only moved into place once the
#' HTTP status and the gzip stream both check out, so an interrupted download
#' can never masquerade as a valid cache entry.
#'
#' @param component One of "visnir", "mir", "soillab", or the NeoSpectra
#'   library files "neospectra_nir", "neospectra_soillab", "neospectra_soilsite"
#' @param cache_dir Local directory for caching (default: ossl_cache_dir())
#' @param force Re-download even if a valid file already exists
#' @return Invisible path to the downloaded file
#' @export
ossl_download_file <- function(component = c("visnir", "mir", "soillab",
                                             "neospectra_nir", "neospectra_soillab",
                                             "neospectra_soilsite"),
                               cache_dir = ossl_cache_dir(),
                               force = FALSE) {
  component <- match.arg(component)
  fname <- .ossl_files[[component]]
  dest  <- file.path(cache_dir, fname)

  if (!force && file.exists(dest)) {
    if (.ossl_gz_ok(dest)) {
      message("OSSL [", component, "] already cached: ", dest)
      return(invisible(dest))
    }
    message("OSSL [", component, "] cache is truncated - re-downloading: ", dest)
  }

  url <- paste0(.OSSL_BASE, fname)
  message("Downloading OSSL [", component, "] from:\n  ", url)

  tmp <- paste0(dest, ".part")
  on.exit(unlink(tmp), add = TRUE)

  resp <- httr::GET(url, httr::write_disk(tmp, overwrite = TRUE),
                    httr::progress())

  if (httr::http_error(resp))
    stop("OSSL download failed for '", fname, "': HTTP ",
         httr::status_code(resp), "\n  ", url)

  if (!.ossl_gz_ok(tmp))
    stop("OSSL download for '", fname, "' is incomplete or corrupt ",
         "(gzip stream does not terminate cleanly). Re-run to retry.")

  if (!file.rename(tmp, dest))
    stop("Could not move the downloaded OSSL file into the cache: ", dest)

  message("Saved to: ", dest)
  invisible(dest)
}

#' Download all three OSSL components (VisNIR, MIR, soillab)
#'
#' @param cache_dir Local directory for caching
#' @param force Re-download even if files already exist
#' @param components Which components to download (default: all three)
#' @return Invisible named list of file paths
#' @export
ossl_download <- function(cache_dir = ossl_cache_dir(),
                          force = FALSE,
                          components = c("visnir", "mir", "soillab")) {
  paths <- lapply(components, ossl_download_file,
                  cache_dir = cache_dir, force = force)
  names(paths) <- components
  invisible(paths)
}

# ---- Internal column helpers ---------------------------------------------

# VisNIR spectral columns: scan_visnir.<wl>_ref   (underscore after "scan")
# Metadata columns use scan.visnir.<key> and are intentionally excluded.
.visnir_cols <- function(df) grep("^scan_visnir\\.", names(df), value = TRUE)

# MIR spectral columns: scan_mir.<wn>_abs
.mir_cols <- function(df) grep("^scan_mir\\.", names(df), value = TRUE)

# Extract numeric position (nm or cm-1) from an OSSL spectral column name.
# Handles both formats that have appeared in OSSL releases:
#   scan_visnir.<wl>_ref   (v1.2 actual)
#   scan.visnir.<wl>_ref   (older drafts / docs)
.wl_from_col <- function(cols) {
  as.numeric(sub(".*\\.([0-9]+(?:\\.[0-9]+)?)_.*$", "\\1", cols))
}

.drop_empty_spectra <- function(M) {
  keep <- rowSums(is.finite(M)) > 0L
  out  <- M[keep, , drop = FALSE]
  attr(out, "n_dropped_empty") <- sum(!keep)
  out
}

.ossl_spectral_availability_report <- function(spec_raw,
                                               sensor_type = c("visnir", "mir"),
                                               valid_ids = NULL) {
  sensor_type <- match.arg(sensor_type)
  dataset_col <- "dataset.code_ascii_txt"
  id_col <- "id.layer_uuid_txt"

  if (!dataset_col %in% names(spec_raw) || !id_col %in% names(spec_raw)) {
    return(NULL)
  }

  dataset_labels <- as.character(spec_raw[[dataset_col]])
  dataset_labels[is.na(dataset_labels) | !nzchar(dataset_labels)] <- "Unknown dataset"

  if (is.null(valid_ids)) {
    valid_ids <- character(0)
  } else {
    valid_ids <- unique(as.character(valid_ids))
  }

  valid_mask <- spec_raw[[id_col]] %in% valid_ids
  counts_raw <- sort(table(dataset_labels), decreasing = TRUE)
  counts_valid <- tapply(as.integer(valid_mask), dataset_labels, sum, simplify = TRUE)
  counts_valid <- counts_valid[names(counts_raw)]
  counts_valid[is.na(counts_valid)] <- 0L

  counts <- data.frame(
    dataset = names(counts_raw),
    n_raw = as.integer(counts_raw),
    n_valid = as.integer(counts_valid),
    stringsAsFactors = FALSE
  )
  counts$n_empty <- counts$n_raw - counts$n_valid
  counts$status <- ifelse(
    counts$n_valid == 0L,
    "no_valid_spectra",
    ifelse(counts$n_empty == 0L, "all_valid", "partial")
  )
  counts$possible_upstream_inconsistency <- sensor_type == "visnir" &
    counts$n_valid == 0L &
    counts$dataset %in% c("AFSIS1.SSL", "AFSIS2.SSL")

  list(
    sensor_type = sensor_type,
    n_rows_raw = nrow(spec_raw),
    n_rows_valid = length(valid_ids),
    n_rows_empty = nrow(spec_raw) - length(valid_ids),
    n_datasets_raw = nrow(counts),
    n_datasets_valid = sum(counts$n_valid > 0L),
    n_datasets_empty = sum(counts$n_valid == 0L),
    counts = counts
  )
}

.ossl_metadata_cols <- function(sensor_type, spec_raw) {
  wanted <- c(
    "dataset.code_ascii_txt",
    "id.scan_local_c"
  )

  if (identical(sensor_type, "visnir")) {
    wanted <- c(
      wanted,
      "scan.visnir.model.name_utf8_txt",
      "scan.visnir.model.code_any_txt",
      "scan.visnir.method.preparation_any_txt"
    )
  } else {
    wanted <- c(
      wanted,
      "scan.mir.model.name_utf8_txt",
      "scan.mir.model.code_any_txt",
      "scan.mir.method.preparation_any_txt"
    )
  }

  intersect(wanted, names(spec_raw))
}

# ---- Load functions -----------------------------------------------------

#' Load a cached OSSL CSV (gzipped) as a data.table
#'
#' @param component "visnir", "mir", or "soillab"
#' @param cache_dir Cache directory
#' @return data.table
#' @export
ossl_load_raw <- function(component = c("visnir", "mir", "soillab"),
                          cache_dir = ossl_cache_dir()) {
  component <- match.arg(component)
  path <- file.path(cache_dir, .ossl_files[[component]])
  if (!file.exists(path)) {
    stop("File not found: ", path, "\nRun ossl_download('", component, "') first.")
  }
  # fread() reads .gz files through R.utils; without it the user gets a
  # data.table error that does not say what to install.
  if (!requireNamespace("R.utils", quietly = TRUE))
    stop("Reading the OSSL .csv.gz files needs the 'R.utils' package: ",
         "install.packages('R.utils')", call. = FALSE)
  message("Loading OSSL [", component, "] ...")
  data.table::fread(path, data.table = FALSE)
}

#' Extract OSSL VisNIR spectra as a matrix aligned to the OSSL standard grid
#'
#' Returns a numeric matrix with rows = samples, columns = wavelengths in nm
#' (350 to 2500, step 2). The row names are the layer UUIDs.
#'
#' @param visnir_df Raw VisNIR data frame from ossl_load_raw("visnir")
#' @return Numeric matrix (n_samples × 1076 wavelengths)
#' @export
ossl_visnir_matrix <- function(visnir_df) {
  cols <- .visnir_cols(visnir_df)
  if (length(cols) == 0) stop("No VisNIR columns found (expected 'scan.visnir.*')")
  wl   <- .wl_from_col(cols)
  # Sort by wavelength
  ord  <- order(wl)
  M    <- as.matrix(visnir_df[, cols[ord]])
  rownames(M) <- visnir_df[["id.layer_uuid_txt"]]
  colnames(M) <- as.character(wl[ord])
  # Replace any zero/negative reflectance with a small positive value
  M[is.finite(M) & M <= 0] <- 1e-6
  .drop_empty_spectra(M)
}

#' Extract OSSL MIR spectra as a matrix aligned to the OSSL standard grid
#'
#' Returns a numeric matrix with rows = samples, columns = wavenumbers in cm-1
#' (600 to 4000, step 2). Values are already in absorbance.
#'
#' @param mir_df Raw MIR data frame from ossl_load_raw("mir")
#' @return Numeric matrix (n_samples × 1701 wavenumbers)
#' @export
ossl_mir_matrix <- function(mir_df) {
  cols <- .mir_cols(mir_df)
  if (length(cols) == 0) stop("No MIR columns found (expected 'scan.mir.*')")
  wn   <- .wl_from_col(cols)
  ord  <- order(wn)
  M    <- as.matrix(mir_df[, cols[ord]])
  rownames(M) <- mir_df[["id.layer_uuid_txt"]]
  colnames(M) <- as.character(wn[ord])
  .drop_empty_spectra(M)
}

#' Extract OSSL Level-1 soil lab data for a set of properties
#'
#' Each property is read from one fixed OSSL column (\code{.ossl_l1_columns}),
#' so its units are known and match \code{ossl_l1_labels}.
#'
#' @param soillab_df Raw soillab data frame from ossl_load_raw("soillab")
#' @param properties Character vector of simplified property names
#'   (from ossl_l1_properties). Default: all.
#' @return Data frame with columns: id.layer_uuid_txt + one column per property
#' @export
ossl_soillab <- function(soillab_df,
                         properties = ossl_l1_properties) {
  id_col <- "id.layer_uuid_txt"
  out <- data.frame(id.layer_uuid_txt = soillab_df[[id_col]],
                    stringsAsFactors = FALSE)
  for (prop in properties) {
    col <- .ossl_l1_columns[prop]
    if (is.na(col) || !col %in% names(soillab_df)) {
      warning("OSSL soillab has no column for '", prop, "' (expected ",
              if (is.na(col)) "a mapping in .ossl_l1_columns" else col, ").",
              call. = FALSE)
      out[[prop]] <- NA_real_
      next
    }
    out[[prop]] <- suppressWarnings(as.numeric(soillab_df[[col]]))
  }
  out
}

#' Join OSSL spectra with soil lab data for model training
#'
#' @param spectra_mat Matrix from ossl_visnir_matrix() or ossl_mir_matrix()
#' @param soillab_df Data frame from ossl_soillab()
#' @param properties Properties to include; default all available
#' @return Data frame with Soil_ID + spectral columns + property columns
#' @export
ossl_join <- function(spectra_mat, soillab_df,
                      properties = ossl_l1_properties) {
  common_ids <- intersect(rownames(spectra_mat),
                          soillab_df[["id.layer_uuid_txt"]])
  if (length(common_ids) == 0)
    stop("No matching layer IDs between spectra and soillab.")

  spec_df <- as.data.frame(spectra_mat[common_ids, , drop = FALSE],
                           check.names = FALSE)
  spec_df[["Soil_ID"]] <- common_ids

  lab_sub <- soillab_df[soillab_df[["id.layer_uuid_txt"]] %in% common_ids, ]
  lab_sub <- lab_sub[match(common_ids, lab_sub[["id.layer_uuid_txt"]]), ]

  props_present <- intersect(properties, names(lab_sub))
  cbind(spec_df[, c("Soil_ID", colnames(spectra_mat[common_ids, , drop = FALSE]))],
        lab_sub[, props_present, drop = FALSE])
}

#' Convenience: load and join all OSSL data for a given sensor type
#'
#' Downloads data if not yet cached, then returns a joined data frame ready
#' for model training.
#'
#' @param sensor_type "visnir", "mir" or "neospectra" (the OSSL NeoSpectra
#'   library, see \code{neospectra_prepare()})
#' @param cache_dir Cache directory
#' @param properties Soil properties to include
#' @param download_if_missing Automatically download if files are absent
#' @return Data frame with Soil_ID + spectral + soil property columns
#' @export
ossl_prepare <- function(sensor_type = c("visnir", "mir", "neospectra"),
                         cache_dir = ossl_cache_dir(),
                         properties = ossl_l1_properties,
                         download_if_missing = TRUE) {
  sensor_type <- match.arg(sensor_type)
  if (sensor_type == "neospectra")
    return(neospectra_prepare(cache_dir, properties, download_if_missing))

  spec_file <- file.path(cache_dir, .ossl_files[[sensor_type]])
  lab_file  <- file.path(cache_dir, .ossl_files[["soillab"]])

  if (download_if_missing) {
    if (!file.exists(spec_file)) ossl_download_file(sensor_type, cache_dir)
    if (!file.exists(lab_file))  ossl_download_file("soillab", cache_dir)
  }

  spec_raw <- ossl_load_raw(sensor_type, cache_dir)
  lab_raw  <- ossl_load_raw("soillab",   cache_dir)

  lab_df <- ossl_soillab(lab_raw, properties = properties)

  M <- if (sensor_type == "visnir") {
    ossl_visnir_matrix(spec_raw)
  } else {
    ossl_mir_matrix(spec_raw)
  }

  availability_report <- .ossl_spectral_availability_report(
    spec_raw = spec_raw,
    sensor_type = sensor_type,
    valid_ids = rownames(M)
  )

  n_dropped_empty <- attr(M, "n_dropped_empty")
  if (isTRUE(n_dropped_empty > 0L)) {
    message("  Dropped ", format(n_dropped_empty, big.mark = ","),
            " rows with no finite ", toupper(sensor_type), " spectral values")
  }

  joined <- ossl_join(M, lab_df, properties = properties)

  id_col <- "id.layer_uuid_txt"
  meta_cols <- .ossl_metadata_cols(sensor_type, spec_raw)
  if (length(meta_cols) > 0L && id_col %in% names(spec_raw)) {
    ds_map      <- spec_raw[, c(id_col, meta_cols), drop = FALSE]
    row_match   <- match(joined[["Soil_ID"]], ds_map[[id_col]])
    for (cn in meta_cols)
      joined[[cn]] <- ds_map[[cn]][row_match]
  }

  attr(joined, "ossl_availability_report") <- availability_report
  attr(joined, "ossl_n_dropped_empty") <- if (length(n_dropped_empty) == 0L ||
                                               is.null(n_dropped_empty)) {
    0L
  } else {
    as.integer(n_dropped_empty)
  }

  joined
}

#' Load the OSSL NeoSpectra library for model training
#'
#' Reads the NeoSpectra handheld NIR library of the OSSL (1350-2550 nm, several
#' units per soil; Sanderman et al. 2023) from the official OSSL storage. One
#' row per scan. \code{sample_group} names the soil, so every scan of a soil
#' stays on one side of any split.
#'
#' By default the 90 African soils (Ghana, Nigeria, Kenya) are left out: they
#' are the external test set of Safanelli et al. (2025, Table 6), not part of
#' their cross-validation.
#'
#' The NeoSpectra lab file holds 22 properties. The soils are KSSL archive
#' samples that are also in the OSSL lab table, under an id that the public
#' files do not link. Properties missing from the NeoSpectra file are taken
#' from the OSSL lab table for soils that match it on six lab values (organic
#' carbon, nitrogen, pH, clay, sand, CEC) within reporting precision and match
#' exactly one OSSL row: 1974 of the 2096 soils, the sample count of the paper's
#' NeoSpectra models. A deliberately shifted match finds none, so the matches
#' are not coincidences.
#'
#' @param cache_dir Cache directory
#' @param properties Soil properties to include
#' @param download_if_missing Automatically download if files are absent
#' @param include_africa Keep the African soils (default \code{FALSE})
#' @return Data frame with \code{Soil_ID} (scan id), spectral columns named by
#'   wavelength, \code{sample_group}, \code{country},
#'   \code{dataset.code_ascii_txt} and one column per property
#' @export
neospectra_prepare <- function(cache_dir = ossl_cache_dir(),
                               properties = ossl_l1_properties,
                               download_if_missing = TRUE,
                               include_africa = FALSE) {
  comps <- c("neospectra_nir", "neospectra_soillab", "neospectra_soilsite")
  for (cp in comps) {
    if (!file.exists(file.path(cache_dir, .ossl_files[[cp]]))) {
      if (!download_if_missing) stop("Missing ", .ossl_files[[cp]], " in ", cache_dir)
      ossl_download_file(cp, cache_dir)
    }
  }
  rd <- function(cp) data.table::fread(file.path(cache_dir, .ossl_files[[cp]]),
                                       showProgress = FALSE, data.table = FALSE)
  nir  <- rd("neospectra_nir"); lab <- rd("neospectra_soillab"); site <- rd("neospectra_soilsite")
  nir  <- nir[!duplicated(nir$id.scan_local_c), , drop = FALSE]
  cols <- grep("^scan_nir[.][0-9]+_ref$", names(nir), value = TRUE)
  wl   <- as.numeric(sub("^scan_nir[.]([0-9]+)_ref$", "\\1", cols))
  M    <- as.matrix(nir[, cols, drop = FALSE]); colnames(M) <- wl
  sid  <- nir$id.sample_local_c
  country <- site$location.country_iso.3166_txt[match(sid, site$id.sample_local_c)]
  keep <- include_africa | !(country %in% c("GHA", "NGA", "KEN"))
  out <- data.frame(Soil_ID = as.character(nir$id.scan_local_c), M, check.names = FALSE,
                    stringsAsFactors = FALSE)
  out$sample_group <- paste0("NEO:", sid)
  out$country <- country
  out$dataset.code_ascii_txt <- "NEOSPECTRA"
  li <- match(sid, lab$id.sample_local_c)
  missing <- character(0)
  for (prop in properties) {
    col <- .ossl_l1_columns[prop]
    if (!is.na(col) && col %in% names(lab)) {
      out[[prop]] <- suppressWarnings(as.numeric(lab[[col]][li]))
    } else {
      out[[prop]] <- NA_real_
      if (!is.na(col)) missing <- c(missing, prop)
    }
  }
  if (length(missing)) {
    ol <- ossl_load_raw("soillab", cache_dir)
    m  <- .neospectra_ossl_match(lab, ol)            # NeoSpectra lab row -> OSSL lab row
    hit <- m[li]
    for (prop in missing) {
      col <- .ossl_l1_columns[[prop]]
      if (col %in% names(ol))
        out[[prop]] <- suppressWarnings(as.numeric(ol[[col]][hit]))
    }
    message("  NeoSpectra: ", length(missing), " properties taken from the OSSL lab table (",
            sum(!is.na(m)), " of ", nrow(lab), " soils matched)")
  }
  out <- out[keep, , drop = FALSE]
  message("  NeoSpectra: ", nrow(out), " scans of ", length(unique(out$sample_group)),
          " soils", if (!include_africa) " (African soils left out)")
  out
}

# Match NeoSpectra lab rows to OSSL lab rows on six lab values. Returns, for
# each NeoSpectra lab row, the OSSL lab row index, or NA unless exactly one
# OSSL row agrees on all six within the NeoSpectra reporting precision.
.neospectra_ossl_match <- function(neo_lab, ossl_lab) {
  keys <- c(oc = 0.051, n.tot = 0.0051, ph.h2o = 0.0051, clay.tot = 0.051,
            sand.tot = 0.051, cec = 0.051)
  cols <- unname(.ossl_l1_columns[names(keys)])
  N <- as.matrix(neo_lab[, cols]); O <- as.matrix(ossl_lab[, cols])
  storage.mode(N) <- "double"; storage.mode(O) <- "double"
  ok_o <- which(stats::complete.cases(O))
  O <- O[ok_o, , drop = FALSE]
  ord <- order(O[, 1]); O <- O[ord, , drop = FALSE]; ok_o <- ok_o[ord]
  vapply(seq_len(nrow(N)), function(i) {
    x <- N[i, ]
    if (anyNA(x)) return(NA_integer_)
    lo <- findInterval(x[1] - keys[1] - 1e-9, O[, 1]) + 1L
    hi <- findInterval(x[1] + keys[1], O[, 1])
    if (hi < lo) return(NA_integer_)
    cand <- lo:hi
    for (j in 2:length(keys))
      cand <- cand[abs(O[cand, j] - x[j]) <= keys[j]]
    if (length(cand) == 1L) ok_o[cand] else NA_integer_
  }, integer(1))
}
