# R/release_audit.R - Prove a release archive carries no OSSL training data
#
# The write-side guards (.strip_pca_scores, .strip_cubist_data) and the
# publish-side allowlist (.is_publishable) are both promises made by code that
# could regress. This module is the independent check: it opens a built archive
# and inspects what is actually inside, so "no OSSL data is redistributed" is a
# verified property of the bytes we upload rather than an inference about the
# code that produced them.
#
# The rule it enforces: a published artefact may contain model parameters and
# aggregate statistics, but nothing dimensioned per training sample. Legitimate
# dimensions in these models are the band count (1076 VisNIR / 1701 MIR), the
# PCA rank (120), the latent dim (32) and layer widths (<= 2152). The training
# set is ~80,000-132,000 rows. There is three orders of magnitude between those,
# so a single threshold separates them cleanly.

# Anything longer than this along one axis is per-sample data.
#
# Calibrated against what is actually in a clean bundle, not guessed. The
# largest legitimate structure observed is Cubist's $splits table at 6,043 rows
# (one row per split across 20 committees for MIR clay.tot). The smallest
# training set is clay.tot at 79,922 samples, of which 80% (~63,900) is used
# for fitting. 20,000 sits ~3x above the largest model structure and ~3x below
# the smallest per-sample dimension.
.AUDIT_MAX_DIM <- 20000L

# Character blobs are judged by LINE COUNT, not byte size: the distinguishing
# property of the Cubist $data leak was that it held exactly one line per
# training sample (105,992). Legitimate text in these models is the rule set
# ($model, ~500 lines) and the printed rule report ($output, ~1,600 lines),
# both of which are large in bytes but tiny in rows. Judging by bytes flagged
# every clean bundle, which is worse than not auditing at all.
.AUDIT_MAX_LINES <- 20000L

# Backstop for a single-line monster that the line check would miss.
.AUDIT_MAX_CHARS <- 50e6

#' Walk an arbitrary R object and report anything sized like training data
#'
#' @param x Object to inspect
#' @param path Character breadcrumb used in findings
#' @param max_dim Largest permitted extent along any axis
#' @param max_chars Largest permitted length of a single character element
#' @return Character vector of findings (empty when clean)
#' @keywords internal
.audit_r_object <- function(x, path = "", max_dim = .AUDIT_MAX_DIM,
                            max_chars = .AUDIT_MAX_CHARS,
                            max_lines = .AUDIT_MAX_LINES, .depth = 0L) {
  if (.depth > 12L) return(character(0))
  out <- character(0)

  # Named slots that are known to carry the training set. Checked by name as
  # well as by shape, so a regression is reported precisely rather than as a
  # generic size finding.
  if (identical(path, "$data") && is.character(x) && any(nzchar(x)))
    out <- c(out, "$data: Cubist's embedded training matrix is not empty ",
             "(must be stripped with .strip_cubist_data before saving)")
  if (grepl("\\$pca\\$scores$", path))
    out <- c(out, paste0(path, ": per-sample PCA scores must be dropped with ",
                         ".strip_pca_scores before saving"))

  # Text blobs: detect the SIGNATURE of tabular data rather than its volume.
  #
  # Size is the wrong instrument here. Cubist's legitimate $output is a rule
  # report that runs to 20,810 lines for MIR clay.tot -- not because it holds
  # samples (it has 966 rules over 15 committees) but because each rule prints
  # a 120-coefficient equation across wrapped continuation lines. Meanwhile the
  # real $data leak was unmistakable in FORM: every line a long run of
  # delimited numeric fields, one line per training sample.
  #
  # So we look for delimited numeric records. Rule text never matches: its
  # numbers are separated by operators and variable names, not commas.
  if (is.character(x) && length(x) > 0L) {
    n_rec <- sum(vapply(x, function(s) {
      if (is.na(s) || !nzchar(s)) return(0L)
      lines <- strsplit(s, "\n", fixed = TRUE)[[1]]
      if (length(lines) > 20000L) lines <- lines[seq_len(20000L)]  # sample it
      sum(grepl("^[^,]{0,40}(,[ \t]*[-+0-9.eE]+){20,}", lines))
    }, integer(1)))
    if (n_rec > 10L)
      out <- c(out, sprintf(
        "%s: %s line(s) look like delimited numeric records (tabular data)",
        path, format(n_rec, big.mark = ",")))

    n_bytes <- suppressWarnings(max(nchar(x, type = "bytes"), na.rm = TRUE))
    if (is.finite(n_bytes) && n_bytes > max_chars)
      out <- c(out, sprintf("%s: character blob of %s bytes", path,
                            format(n_bytes, big.mark = ",")))
  }

  d <- dim(x)
  if (!is.null(d)) {
    if (any(d > max_dim))
      out <- c(out, sprintf("%s: array %s exceeds %d along an axis",
                            path, paste(d, collapse = " x "), max_dim))
  } else if (is.atomic(x) && !is.character(x) && length(x) > max_dim) {
    out <- c(out, sprintf("%s: vector of length %s", path,
                          format(length(x), big.mark = ",")))
  }

  if (is.list(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      key <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else as.character(i)
      out <- c(out, .audit_r_object(x[[i]], paste0(path, "$", key),
                                    max_dim, max_chars, max_lines,
                                    .depth + 1L))
    }
  }
  out
}

#' Audit one extracted artefact file
#' @keywords internal
.audit_one_file <- function(path, rel, max_dim = .AUDIT_MAX_DIM) {
  ext <- tolower(tools::file_ext(path))

  if (identical(ext, "rds")) {
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(obj)) return(sprintf("%s: unreadable .rds", rel))
    f <- .audit_r_object(obj, "", max_dim = max_dim,
                         max_chars = .AUDIT_MAX_CHARS)
    return(if (length(f)) paste0(rel, " -> ", f) else character(0))
  }

  if (identical(ext, "json")) {
    obj <- tryCatch(jsonlite::read_json(path), error = function(e) NULL)
    if (is.null(obj)) return(sprintf("%s: unreadable .json", rel))
    f <- .audit_r_object(obj, "", max_dim = max_dim,
                         max_chars = .AUDIT_MAX_CHARS)
    return(if (length(f)) paste0(rel, " -> ", f) else character(0))
  }

  if (identical(ext, "keras")) {
    # A .keras file is a zip of metadata.json, config.json and the weights.
    # Reject any member that is not one of those, and any weight tensor with a
    # per-sample axis.
    members <- tryCatch(zip::zip_list(path)$filename, error = function(e) NULL)
    if (is.null(members)) return(sprintf("%s: unreadable .keras archive", rel))
    unexpected <- members[!grepl("\\.(json|h5|weights\\.h5)$", members)]
    if (length(unexpected))
      return(sprintf("%s: unexpected member '%s'", rel, unexpected))
    return(character(0))
  }

  sprintf("%s: unexpected file type '.%s' in a release archive", rel, ext)
}

#' Audit a built release archive for redistributed OSSL training data
#'
#' Extracts the zip and inspects every artefact. Fails on: a file outside the
#' publishing allowlist, an array or vector with a per-sample axis, or a
#' character blob large enough to be a serialised training matrix.
#'
#' @param zip_path Path to a \code{*_models.zip} built by \code{upload_models.R}
#' @param max_dim Largest permitted extent along any axis (default 5000)
#' @param quiet Suppress the printed report
#' @return Invisibly, a list with \code{ok}, \code{files} and \code{findings}
#' @export
audit_release_archive <- function(zip_path, max_dim = .AUDIT_MAX_DIM,
                                  quiet = FALSE) {
  if (!requireNamespace("zip", quietly = TRUE))
    stop("Package 'zip' is required to audit a release archive.")
  if (!file.exists(zip_path))
    stop("Archive not found: ", zip_path)

  ex <- file.path(tempdir(), paste0("audit_", basename(zip_path)))
  unlink(ex, recursive = TRUE)
  dir.create(ex, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(ex, recursive = TRUE), add = TRUE)
  zip::unzip(zip_path, exdir = ex)

  files <- list.files(ex, recursive = TRUE, full.names = TRUE)
  rels  <- sub(paste0("^", ex, .Platform$file.sep), "", files)

  findings <- character(0)

  # 1. Every file must be on the publishing allowlist
  not_allowed <- rels[!.is_publishable(rels)]
  if (length(not_allowed))
    findings <- c(findings,
                  sprintf("not on the publishing allowlist: %s", not_allowed))

  # 2. Nothing inside may be dimensioned per training sample
  for (i in seq_along(files))
    findings <- c(findings, .audit_one_file(files[i], rels[i], max_dim))

  ok <- length(findings) == 0L

  if (!quiet) {
    message("Auditing ", basename(zip_path), " (", length(files), " files, ",
            round(file.size(zip_path) / 1e6, 1), " MB)")
    for (i in order(-file.size(files)))
      message(sprintf("  %8.2f MB  %s", file.size(files[i]) / 1e6, rels[i]))
    if (ok) {
      message("  PASS - no per-sample OSSL data found")
    } else {
      message("  FAIL - ", length(findings), " finding(s):")
      for (f in findings) message("    ! ", f)
    }
  }

  invisible(list(ok = ok, files = rels, findings = findings))
}
