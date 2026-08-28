# Replicate-aware OSSL subsampling helpers used by train_ossl.R.

.normalize_ossl_sensor <- function(x,
                                   missing_label = "Unknown / unspecified sensor") {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- missing_label
  x
}

.allocate_proportional_counts <- function(sizes, n) {
  sizes <- as.integer(sizes)
  out   <- integer(length(sizes))

  if (length(sizes) == 0L || n <= 0L || sum(sizes) <= 0L)
    return(out)

  raw <- n * sizes / sum(sizes)
  out <- floor(raw)

  left <- n - sum(out)
  if (left > 0L) {
    frac <- raw - out
    ord  <- order(frac, sizes, decreasing = TRUE)
    for (i in ord) {
      if (left == 0L) break
      if (out[i] < sizes[i]) {
        out[i] <- out[i] + 1L
        left   <- left - 1L
      }
    }
  }

  if (sum(out) < n) {
    ord <- order(sizes - out, decreasing = TRUE)
    for (i in ord) {
      if (sum(out) == n) break
      if (out[i] < sizes[i]) out[i] <- out[i] + 1L
    }
  }

  out
}

.allocate_balanced_counts <- function(capacities, n, base = NULL) {
  capacities <- as.integer(capacities)
  out <- if (is.null(base)) integer(length(capacities)) else as.integer(base)

  if (length(capacities) == 0L || n <= 0L)
    return(out)

  out <- pmin(out, capacities)
  target_n <- min(as.integer(n), sum(capacities))
  left <- target_n - sum(out)

  if (left <= 0L)
    return(out)

  while (left > 0L) {
    active <- which(out < capacities)
    if (length(active) == 0L) break

    min_count <- min(out[active])
    tied <- active[out[active] == min_count]
    tied <- tied[order(capacities[tied] - out[tied], decreasing = TRUE)]

    for (i in tied) {
      if (left == 0L) break
      if (out[i] < capacities[i]) {
        out[i] <- out[i] + 1L
        left   <- left - 1L
      }
    }
  }

  out
}

.extract_group_value <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 1L) return(NA_real_)
  as.numeric(x[1L])
}

.assign_quantile_bins <- function(x,
                                  n_bins = 4L,
                                  prefix = "Q",
                                  missing_label = "missing") {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(missing_label, length(x))
  ok  <- is.finite(x)

  if (!any(ok))
    return(out)

  n_bins <- max(1L, as.integer(n_bins))
  probs  <- seq(0, 1, length.out = n_bins + 1L)
  brks   <- unique(stats::quantile(
    x[ok],
    probs  = probs,
    na.rm  = TRUE,
    names  = FALSE,
    type   = 7
  ))

  if (length(brks) <= 1L) {
    out[ok] <- paste0(prefix, "1")
    return(out)
  }

  out[ok] <- paste0(
    prefix,
    cut(
      x[ok],
      breaks         = brks,
      include.lowest = TRUE,
      labels         = FALSE
    )
  )
  out
}

select_ossl_training_subset <- function(df,
                                        sample_n,
                                        sample_id_col = "Soil_ID",
                                        sensor_col    = NULL,
                                        dataset_col   = NULL,
                                        target_col    = NULL,
                                        target_bins   = 4L,
                                        balance_datasets = !is.null(dataset_col),
                                        balance_target   = !is.null(target_col),
                                        seed          = 42L) {
  if (is.null(sample_n))
    return(list(
      data               = df,
      sample_ids         = unique(as.character(df[[sample_id_col]])),
      sample_meta        = NULL,
      n_selected_samples = length(unique(df[[sample_id_col]])),
      n_selected_rows    = nrow(df),
      sensor_counts      = integer(0),
      dataset_counts     = integer(0),
      target_counts      = integer(0),
      kept_replicates    = anyDuplicated(df[[sample_id_col]]) > 0L,
      all_sensors_kept   = NA,
      all_datasets_kept  = NA,
      used_target_balance = FALSE,
      target_col         = NULL
    ))

  if (!sample_id_col %in% names(df))
    stop("Sample ID column not found: ", sample_id_col)

  sample_n <- as.integer(sample_n)
  if (!is.finite(sample_n) || sample_n < 1L)
    stop("sample_n must be a positive integer.")

  sample_ids <- as.character(df[[sample_id_col]])
  unique_ids <- unique(sample_ids)

  if (sample_n >= length(unique_ids)) {
    sensor_counts <- integer(0)
    dataset_counts <- integer(0)
    target_counts <- integer(0)

    if (!is.null(sensor_col) && sensor_col %in% names(df)) {
      sensor_vals <- .normalize_ossl_sensor(df[[sensor_col]])
      sensor_counts <- sort(table(sensor_vals[match(unique_ids, sample_ids)]),
                            decreasing = TRUE)
    }
    if (!is.null(dataset_col) && dataset_col %in% names(df)) {
      dataset_vals <- as.character(df[[dataset_col]])
      dataset_counts <- sort(table(dataset_vals[match(unique_ids, sample_ids)]),
                             decreasing = TRUE)
    }
    if (!is.null(target_col) && target_col %in% names(df)) {
      target_vals <- suppressWarnings(as.numeric(df[[target_col]]))
      target_by_id <- tapply(target_vals, sample_ids, .extract_group_value, simplify = TRUE)
      target_bins_out <- .assign_quantile_bins(
        target_by_id[unique_ids],
        n_bins = target_bins
      )
      target_counts <- sort(table(target_bins_out[target_bins_out != "missing"]),
                            decreasing = FALSE)
    }

    return(list(
      data               = df,
      sample_ids         = unique_ids,
      sample_meta        = NULL,
      n_selected_samples = length(unique_ids),
      n_selected_rows    = nrow(df),
      sensor_counts      = sensor_counts,
      dataset_counts     = dataset_counts,
      target_counts      = target_counts,
      kept_replicates    = anyDuplicated(sample_ids) > 0L,
      all_sensors_kept   = TRUE,
      all_datasets_kept  = TRUE,
      used_target_balance = FALSE,
      target_col         = if (!is.null(target_col) && target_col %in% names(df)) target_col else NULL
    ))
  }

  group_meta <- data.frame(
    sample_id = unique_ids,
    stringsAsFactors = FALSE
  )

  if (!is.null(sensor_col) && sensor_col %in% names(df)) {
    sensor_vals <- .normalize_ossl_sensor(df[[sensor_col]])
    sensor_sets <- tapply(sensor_vals, sample_ids, function(x) unique(x), simplify = FALSE)
    group_meta$primary_sensor <- vapply(sensor_sets[unique_ids], `[`, character(1), 1L)
    group_meta$sensor_signature <- vapply(
      sensor_sets[unique_ids],
      function(x) paste(sort(unique(x)), collapse = " + "),
      character(1)
    )
  }

  if (!is.null(dataset_col) && dataset_col %in% names(df)) {
    dataset_vals <- as.character(df[[dataset_col]])
    dataset_sets <- tapply(dataset_vals, sample_ids, function(x) unique(x), simplify = FALSE)
    group_meta$primary_dataset <- vapply(dataset_sets[unique_ids], `[`, character(1), 1L)
  }

  if (!is.null(target_col) && target_col %in% names(df)) {
    target_vals <- suppressWarnings(as.numeric(df[[target_col]]))
    target_sets <- tapply(target_vals, sample_ids, .extract_group_value, simplify = TRUE)
    group_meta$target_value <- as.numeric(target_sets[unique_ids])
  }

  if (!"primary_sensor" %in% names(group_meta))
    group_meta$primary_sensor <- rep("all", nrow(group_meta))
  if (!"primary_dataset" %in% names(group_meta))
    group_meta$primary_dataset <- rep("all", nrow(group_meta))
  if (!"target_value" %in% names(group_meta))
    group_meta$target_value <- rep(NA_real_, nrow(group_meta))

  use_target_balance <- isTRUE(balance_target) && any(is.finite(group_meta$target_value))
  if (use_target_balance) {
    target_bins_out <- rep("missing", nrow(group_meta))
    if (isTRUE(balance_datasets) && !"all" %in% unique(group_meta$primary_dataset)) {
      idx_by_dataset <- split(seq_len(nrow(group_meta)), group_meta$primary_dataset)
      for (idx in idx_by_dataset) {
        target_bins_out[idx] <- .assign_quantile_bins(
          group_meta$target_value[idx],
          n_bins = target_bins
        )
      }
    } else {
      target_bins_out <- .assign_quantile_bins(
        group_meta$target_value,
        n_bins = target_bins
      )
    }
    group_meta$target_bin <- target_bins_out
  } else {
    group_meta$target_bin <- rep("all", nrow(group_meta))
  }

  group_meta$stratum <- paste(group_meta$primary_sensor,
                              group_meta$primary_dataset,
                              group_meta$target_bin,
                              sep = " || ")

  set.seed(as.integer(seed))

  selected_ids <- character(0)
  selection_meta <- group_meta
  if (isTRUE(use_target_balance)) {
    selection_meta <- group_meta[
      is.finite(group_meta$target_value),
      ,
      drop = FALSE
    ]
  }
  if (nrow(selection_meta) == 0L)
    selection_meta <- group_meta

  sensor_levels <- unique(selection_meta$primary_sensor)
  if (!all(sensor_levels == "all") && sample_n >= length(sensor_levels)) {
    mandatory_ids <- unlist(lapply(
      split(selection_meta$sample_id, selection_meta$primary_sensor),
      function(ids) sample(ids, 1L)
    ), use.names = FALSE)
    selected_ids <- unique(mandatory_ids)
  }

  remaining <- selection_meta[!selection_meta$sample_id %in% selected_ids, , drop = FALSE]
  remaining_n <- sample_n - length(selected_ids)

  if (remaining_n > 0L) {
    if (isTRUE(balance_datasets) && !"all" %in% unique(selection_meta$primary_dataset)) {
      dataset_sizes <- table(selection_meta$primary_dataset)
      dataset_levels <- names(dataset_sizes)
      base_counts <- table(factor(
        selection_meta$primary_dataset[selection_meta$sample_id %in% selected_ids],
        levels = dataset_levels
      ))
      target_counts <- .allocate_balanced_counts(
        capacities = as.integer(dataset_sizes),
        n          = sample_n,
        base       = as.integer(base_counts)
      )
      names(target_counts) <- dataset_levels
      remaining_counts <- pmax(target_counts - as.integer(base_counts), 0L)
      names(remaining_counts) <- dataset_levels

      sampled_ids <- unlist(lapply(dataset_levels, function(ds) {
        n_take <- remaining_counts[[ds]]
        if (n_take < 1L) return(character(0))

        ds_all  <- selection_meta[selection_meta$primary_dataset == ds, , drop = FALSE]
        ds_rows <- remaining[remaining$primary_dataset == ds, , drop = FALSE]

        if (isTRUE(use_target_balance)) {
          bin_levels <- sort(unique(ds_all$target_bin))
          bin_capacities <- table(factor(ds_all$target_bin, levels = bin_levels))
          bin_base <- table(factor(
            ds_all$target_bin[ds_all$sample_id %in% selected_ids],
            levels = bin_levels
          ))
          bin_targets <- .allocate_balanced_counts(
            capacities = as.integer(bin_capacities),
            n          = target_counts[[ds]],
            base       = as.integer(bin_base)
          )
          names(bin_targets) <- bin_levels
          bin_remaining <- pmax(bin_targets - as.integer(bin_base), 0L)
          names(bin_remaining) <- bin_levels

          unlist(lapply(bin_levels, function(bin_id) {
            n_take_bin <- bin_remaining[[bin_id]]
            if (n_take_bin < 1L) return(character(0))

            bin_rows <- ds_rows[ds_rows$target_bin == bin_id, , drop = FALSE]
            if (nrow(bin_rows) < 1L) return(character(0))

            sensor_groups <- split(bin_rows$sample_id, bin_rows$primary_sensor)
            alloc_sensor <- .allocate_proportional_counts(
              lengths(sensor_groups),
              n_take_bin
            )
            names(alloc_sensor) <- names(sensor_groups)

            unlist(Map(
              function(ids, n_take_sensor) {
                if (n_take_sensor < 1L) return(character(0))
                sample(ids, n_take_sensor)
              },
              sensor_groups,
              alloc_sensor
            ), use.names = FALSE)
          }), use.names = FALSE)
        } else {
          strata  <- split(ds_rows$sample_id, ds_rows$stratum)
          alloc   <- .allocate_proportional_counts(lengths(strata), n_take)
          names(alloc) <- names(strata)

          unlist(Map(
            function(ids, n_take_stratum) {
              if (n_take_stratum < 1L) return(character(0))
              sample(ids, n_take_stratum)
            },
            strata,
            alloc
          ), use.names = FALSE)
        }
      }), use.names = FALSE)
    } else {
      strata <- split(remaining$sample_id, remaining$stratum)
      alloc  <- .allocate_proportional_counts(lengths(strata), remaining_n)
      names(alloc) <- names(strata)

      sampled_ids <- unlist(Map(
        function(ids, n_take) {
          if (n_take < 1L) return(character(0))
          sample(ids, n_take)
        },
        strata,
        alloc
      ), use.names = FALSE)
    }

    selected_ids <- unique(c(selected_ids, sampled_ids))
  }

  if (length(selected_ids) < sample_n) {
    pool <- setdiff(selection_meta$sample_id, selected_ids)
    if (length(pool) < (sample_n - length(selected_ids)))
      pool <- setdiff(group_meta$sample_id, selected_ids)
    selected_ids <- c(selected_ids, sample(pool, sample_n - length(selected_ids)))
  }

  selected_ids <- selected_ids[seq_len(sample_n)]
  keep_rows    <- sample_ids %in% selected_ids
  out_df       <- df[keep_rows, , drop = FALSE]
  out_meta     <- group_meta[match(selected_ids, group_meta$sample_id), , drop = FALSE]
  ref_meta     <- selection_meta

  sensor_counts <- sort(table(out_meta$primary_sensor), decreasing = TRUE)
  dataset_counts <- sort(table(out_meta$primary_dataset), decreasing = TRUE)
  target_counts <- sort(
    table(out_meta$target_bin[out_meta$target_bin != "all" & out_meta$target_bin != "missing"]),
    decreasing = FALSE
  )
  all_sensors_kept <- setequal(unique(ref_meta$primary_sensor),
                               unique(out_meta$primary_sensor))
  all_datasets_kept <- setequal(unique(ref_meta$primary_dataset),
                                unique(out_meta$primary_dataset))

  # Preserve original row order and keep all spectra for the selected sample IDs.
  out_df <- out_df[order(match(as.character(out_df[[sample_id_col]]), selected_ids)), , drop = FALSE]

  list(
    data               = out_df,
    sample_ids         = selected_ids,
    sample_meta        = out_meta,
    n_selected_samples = length(selected_ids),
    n_selected_rows    = nrow(out_df),
    sensor_counts      = sensor_counts,
    dataset_counts     = dataset_counts,
    target_counts      = target_counts,
    kept_replicates    = anyDuplicated(as.character(out_df[[sample_id_col]])) > 0L,
    all_sensors_kept   = all_sensors_kept,
    all_datasets_kept  = all_datasets_kept,
    used_target_balance = use_target_balance,
    target_col         = if (isTRUE(use_target_balance)) target_col else NULL
  )
}
