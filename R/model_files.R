.model_file_extensions <- c(".keras", ".h5")

.model_file_pattern <- "\\.(keras|h5)$"

.default_model_path <- function(path_base) {
  paste0(path_base, ".keras")
}

.resolve_model_path <- function(path_base,
                                must_exist = TRUE,
                                extensions = .model_file_extensions) {
  candidates <- paste0(path_base, extensions)

  if (!must_exist)
    return(candidates[[1L]])

  hits <- candidates[file.exists(candidates)]
  if (length(hits) > 0L)
    return(hits[[1L]])

  stop(
    "Model file not found for base '", path_base, "'. Checked: ",
    paste(basename(candidates), collapse = ", ")
  )
}

.list_model_files <- function(path, recursive = FALSE, full.names = FALSE) {
  if (!dir.exists(path))
    return(character(0))

  list.files(
    path,
    pattern = .model_file_pattern,
    recursive = recursive,
    full.names = full.names
  )
}
