#!/usr/bin/env Rscript
# Split an older models/<family>/cv10_summary.json into cv10/<prop>.json, the
# per-property layout train_ossl.R --cv10 resumes from.
for (fs in Sys.glob(file.path(if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "models", "*", "cv10_summary.json"))) {
  cv <- jsonlite::read_json(fs)$cv10
  d <- file.path(dirname(fs), "cv10"); dir.create(d, showWarnings = FALSE)
  for (p in names(cv)) {
    f <- file.path(d, paste0(p, ".json"))
    if (!file.exists(f) && !is.null(cv[[p]]))
      jsonlite::write_json(cv[[p]], f, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  }
  message(fs, ": ", length(cv), " properties")
}
