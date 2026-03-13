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
  soillab = "ossl_soillab_L1_v1.2.csv.gz"
)

#' Default local cache directory for OSSL data
#' @export
ossl_cache_dir <- function() {
  d <- file.path(tools::R_user_dir("autoSpectra", "data"), "ossl_v1.2")
  dir_create(d)
  d
}

#' Download one OSSL component file to the cache directory
#'
#' @param component One of "visnir", "mir", or "soillab"
#' @param cache_dir Local directory for caching (default: ossl_cache_dir())
#' @param force Re-download even if file already exists
#' @return Invisible path to the downloaded file
#' @export
ossl_download_file <- function(component = c("visnir", "mir", "soillab"),
                               cache_dir = ossl_cache_dir(),
                               force = FALSE) {
  component <- match.arg(component)
  fname <- .ossl_files[[component]]
  dest  <- file.path(cache_dir, fname)
  if (file.exists(dest) && !force) {
    message("OSSL [", component, "] already cached: ", dest)
    return(invisible(dest))
  }
  url <- paste0(.OSSL_BASE, fname)
  message("Downloading OSSL [", component, "] from:\n  ", url)
  httr::GET(url, httr::write_disk(dest, overwrite = TRUE),
            httr::progress())
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

# VisNIR columns: scan.visnir.XXX_ref
.visnir_cols <- function(df) grep("^scan\\.visnir\\.", names(df), value = TRUE)

# MIR columns: scan.mir.XXX_abs
.mir_cols <- function(df) grep("^scan\\.mir\\.", names(df), value = TRUE)

# Extract numeric position (nm or cm-1) from an OSSL spectral column name
.wl_from_col <- function(cols) {
  as.numeric(gsub("^scan\\.(visnir|mir)\\.(\\d+(?:\\.\\d+)?)_.*$", "\\2", cols))
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
  M
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
  M
}

#' Extract OSSL Level-1 soil lab data for a set of properties
#'
#' Selects the best representative column for each requested property by
#' picking the column with the fewest NAs among candidates matching the
#' simplified property name prefix.
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
    # Find candidate columns: match on property prefix (e.g. "oc." or "oc_")
    pat  <- paste0("^", gsub("\\.", "\\\\.", prop), "[._]")
    cands <- grep(pat, names(soillab_df), value = TRUE)
    if (length(cands) == 0) {
      out[[prop]] <- NA_real_
      next
    }
    # Pick candidate with most non-NA values
    na_counts <- sapply(cands, function(cn) sum(is.na(soillab_df[[cn]])))
    best <- cands[which.min(na_counts)]
    out[[prop]] <- suppressWarnings(as.numeric(soillab_df[[best]]))
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
#' @param sensor_type "visnir" or "mir"
#' @param cache_dir Cache directory
#' @param properties Soil properties to include
#' @param download_if_missing Automatically download if files are absent
#' @return Data frame with Soil_ID + spectral + soil property columns
#' @export
ossl_prepare <- function(sensor_type = c("visnir", "mir"),
                         cache_dir = ossl_cache_dir(),
                         properties = ossl_l1_properties,
                         download_if_missing = TRUE) {
  sensor_type <- match.arg(sensor_type)

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

  ossl_join(M, lab_df, properties = properties)
}
