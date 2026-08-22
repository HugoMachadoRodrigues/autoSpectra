# R/train_cv.R — Cross-validation helpers for autoSpectra
#
# These functions are in the package (not in train_ossl.R) so that callr
# worker subprocesses can find them after library(autoSpectra) / load_all().

# ---- LODO k-fold indices -----------------------------------------------------

#' Create Leave-One-Dataset-Out (LODO) cross-validation fold indices
#'
#' Assigns each OSSL contributing dataset to one fold.  If there are fewer
#' unique datasets than \code{k}, \code{k} is automatically reduced so that
#' every fold contains at least one dataset (true LODO — no random fallback).
#'
#' @param dataset_col Character vector of dataset labels (one per sample).
#' @param k           Maximum number of folds (default 10).
#' @param seed        RNG seed for random dataset ordering (default 42).
#' @return List of length \code{min(k, n_datasets)}; each element is an integer
#'   vector of test-set row indices.
#' @export
stratified_kfold <- function(dataset_col, k = 10L, seed = 42L) {
  set.seed(seed)
  datasets <- unique(dataset_col)
  k        <- min(k, length(datasets))   # true LODO: one fold per dataset
  folds    <- vector("list", k)
  ds_order <- sample(datasets)
  ds_fold  <- setNames(((seq_along(ds_order) - 1L) %% k) + 1L, ds_order)
  for (fold in seq_len(k))
    folds[[fold]] <- which(ds_fold[dataset_col] == fold)
  folds
}

.row_stratified_kfold <- function(strata, k = 10L, seed = 42L) {
  set.seed(seed)
  n <- length(strata)
  if (n < 1L)
    return(list(integer(0)))

  k <- max(1L, min(suppressWarnings(as.integer(k)[1L]), n))
  folds <- vector("list", k)
  strata <- as.character(strata)
  strata[is.na(strata) | !nzchar(strata)] <- "__MISSING__"
  strata_levels <- sample(unique(strata))

  for (lev in strata_levels) {
    idx <- sample(which(strata == lev))
    if (length(idx) < 1L)
      next
    fold_ids <- ((seq_along(idx) - 1L) %% k) + 1L
    for (fold in seq_len(k)) {
      folds[[fold]] <- c(folds[[fold]], idx[fold_ids == fold])
    }
  }

  lapply(folds, function(x) sort(unique(as.integer(x))))
}

.grouped_stratified_kfold <- function(group_col, y = NULL, k = 10L, seed = 42L) {
  group_col <- as.character(group_col)
  missing_group <- is.na(group_col) | !nzchar(group_col)
  if (any(missing_group)) {
    group_col[missing_group] <- paste0("__ROW__", which(missing_group))
  }
  group_levels <- unique(group_col)
  if (length(group_levels) < 1L)
    return(list(integer(0)))

  group_target <- vapply(group_levels, function(g) {
    vals <- suppressWarnings(as.numeric(y[group_col == g]))
    vals <- vals[is.finite(vals)]
    if (length(vals) < 1L) NA_real_ else stats::median(vals)
  }, numeric(1))
  group_bins <- make_target_bins(group_target, n_bins = 4L)
  group_bins[!is.finite(group_bins)] <- 0L
  group_folds <- .row_stratified_kfold(
    strata = paste0("bin_", group_bins),
    k = min(k, length(group_levels)),
    seed = seed
  )
  lapply(group_folds, function(idx_group) {
    keep_groups <- group_levels[idx_group]
    which(group_col %in% keep_groups)
  })
}

.grouped_split_idx <- function(group_col, y, seed = 42L,
                               p_train = 0.6, p_cal = 0.2, p_test = 0.2) {
  stopifnot(abs(p_train + p_cal + p_test - 1) < 1e-9)
  group_col <- as.character(group_col)
  missing_group <- is.na(group_col) | !nzchar(group_col)
  if (any(missing_group)) {
    group_col[missing_group] <- paste0("__ROW__", which(missing_group))
  }
  group_levels <- unique(group_col)
  if (length(group_levels) < 2L) {
    idx <- split_idx(length(group_col), seed = seed, p_train = p_train, p_cal = p_cal, p_test = p_test)
    return(idx)
  }

  group_target <- vapply(group_levels, function(g) {
    vals <- suppressWarnings(as.numeric(y[group_col == g]))
    vals <- vals[is.finite(vals)]
    if (length(vals) < 1L) NA_real_ else stats::median(vals)
  }, numeric(1))
  group_bins <- make_target_bins(group_target, n_bins = 4L)
  group_bins[!is.finite(group_bins)] <- 0L

  set.seed(seed)
  train_groups <- integer(0)
  calib_groups <- integer(0)
  test_groups  <- integer(0)
  probs <- c(train = p_train, calib = p_cal, test = p_test)

  for (bin in unique(group_bins)) {
    idx <- which(group_bins == bin)
    idx <- sample(idx)
    m <- length(idx)
    raw_counts <- probs * m
    counts <- floor(raw_counts)
    remainder <- m - sum(counts)
    if (remainder > 0L) {
      ord <- order(raw_counts - counts, decreasing = TRUE)
      for (j in seq_len(remainder)) {
        counts[ord[((j - 1L) %% length(ord)) + 1L]] <- counts[ord[((j - 1L) %% length(ord)) + 1L]] + 1L
      }
    }
    if (m >= 3L) {
      zero_idx <- which(counts == 0L & probs > 0)
      for (z in zero_idx) {
        donor <- which.max(counts)
        if (counts[donor] > 1L) {
          counts[donor] <- counts[donor] - 1L
          counts[z] <- counts[z] + 1L
        }
      }
    }
    i1 <- counts[[1L]]
    i2 <- counts[[2L]]
    train_groups <- c(train_groups, idx[seq_len(i1)])
    if (i2 > 0L)
      calib_groups <- c(calib_groups, idx[seq.int(i1 + 1L, i1 + i2)])
    if ((m - i1 - i2) > 0L)
      test_groups <- c(test_groups, idx[seq.int(i1 + i2 + 1L, m)])
  }

  if (length(calib_groups) < 1L && length(train_groups) > 1L) {
    calib_groups <- train_groups[[length(train_groups)]]
    train_groups <- train_groups[-length(train_groups)]
  }
  if (length(test_groups) < 1L && length(train_groups) > 1L) {
    test_groups <- train_groups[[length(train_groups)]]
    train_groups <- train_groups[-length(train_groups)]
  }

  list(
    train = which(group_col %in% group_levels[train_groups]),
    calib = which(group_col %in% group_levels[calib_groups]),
    test  = which(group_col %in% group_levels[test_groups])
  )
}

.resolve_cv_folds <- function(dataset_col, y = NULL, k = 10L, seed = 42L,
                              strategy = c("grouped_sample_kfold", "lodo_dataset"),
                              group_col = NULL) {
  strategy <- match.arg(strategy)
  dataset_levels <- sort(unique(as.character(dataset_col)))
  dataset_levels <- dataset_levels[!is.na(dataset_levels) & nzchar(dataset_levels)]

  if (identical(strategy, "grouped_sample_kfold")) {
    if (is.null(group_col))
      group_col <- seq_along(dataset_col)
    return(list(
      folds = .grouped_stratified_kfold(group_col = group_col, y = y, k = k, seed = seed),
      strategy = "grouped_sample_kfold",
      strategy_detail = "Grouped stratified K-fold by Soil_ID",
      dataset_levels = dataset_levels
    ))
  }

  if (length(dataset_levels) >= 2L) {
    return(list(
      folds = stratified_kfold(dataset_col, k = k, seed = seed),
      strategy = "lodo_dataset",
      strategy_detail = "LODO — stratified by OSSL contributing dataset",
      dataset_levels = dataset_levels
    ))
  }

  if (length(dataset_levels) < 1L) {
    return(list(
      folds = list(integer(0)),
      strategy = "unavailable",
      strategy_detail = "No dataset labels available for cross-validation",
      dataset_levels = dataset_levels
    ))
  }

  target_bins <- make_target_bins(y, n_bins = 4L)
  target_bins[!is.finite(target_bins)] <- 0L
  folds <- .row_stratified_kfold(
    strata = paste0("bin_", target_bins),
    k = k,
    seed = seed
  )
  list(
    folds = folds,
    strategy = "within_dataset_kfold",
    strategy_detail = paste0("Stratified K-fold within dataset (", dataset_levels[[1L]], ")"),
    dataset_levels = dataset_levels
  )
}

resolve_lodo_fold_count <- function(dataset_col,
                                    requested_k = 10L,
                                    quick = FALSE,
                                    quick_full_lodo_threshold = 5L) {
  datasets <- unique(as.character(dataset_col))
  n_ds <- length(datasets)

  if (n_ds < 1L)
    return(1L)

  requested_k <- suppressWarnings(as.integer(requested_k))
  if (!is.finite(requested_k) || requested_k < 1L)
    requested_k <- n_ds

  if (isTRUE(quick) && n_ds <= quick_full_lodo_threshold)
    return(n_ds)

  min(requested_k, n_ds)
}

.collapse_unique_labels <- function(x) {
  x <- unique(as.character(x))
  x <- x[nzchar(x) & !is.na(x)]
  if (length(x) < 1L) NA_character_ else paste(sort(x), collapse = " + ")
}

.summarise_fold_metric_rows <- function(fold_rows) {
  metric_names <- c("RMSE", "Bias", "R2", "RPIQ", "CCC")
  out <- list(
    n_fold_results    = nrow(fold_rows),
    n_unique_datasets = length(unique(fold_rows$held_out_dataset))
  )
  for (metric in metric_names) {
    vals <- fold_rows[[metric]]
    out[[paste0(metric, "_mean")]] <- mean(vals, na.rm = TRUE)
    out[[paste0(metric, "_sd")]] <- if (nrow(fold_rows) > 1L) {
      stats::sd(vals, na.rm = TRUE)
    } else {
      NA_real_
    }
  }
  out
}

.build_fold_group_summary <- function(fold_summary, strategy = "lodo_dataset") {
  rows <- list()

  add_group <- function(grouping, group, idx) {
    idx <- which(idx %in% TRUE)
    if (length(idx) < 1L)
      return(invisible(NULL))
    sm <- .summarise_fold_metric_rows(fold_summary[idx, , drop = FALSE])
    rows[[length(rows) + 1L]] <<- data.frame(
      grouping          = grouping,
      group             = group,
      n_fold_results    = sm$n_fold_results,
      n_unique_datasets = sm$n_unique_datasets,
      RMSE_mean         = sm$RMSE_mean,
      RMSE_sd           = sm$RMSE_sd,
      Bias_mean         = sm$Bias_mean,
      Bias_sd           = sm$Bias_sd,
      R2_mean           = sm$R2_mean,
      R2_sd             = sm$R2_sd,
      RPIQ_mean         = sm$RPIQ_mean,
      RPIQ_sd           = sm$RPIQ_sd,
      CCC_mean          = sm$CCC_mean,
      CCC_sd            = sm$CCC_sd,
      stringsAsFactors  = FALSE
    )
    invisible(NULL)
  }

  if (identical(strategy, "lodo_dataset")) {
    held_out_dataset <- as.character(fold_summary$held_out_dataset)
    with_icraf <- grepl("ICRAF", held_out_dataset, fixed = TRUE)
    add_group("held_out_dataset", "with_icraf", with_icraf)
    add_group("held_out_dataset", "without_icraf", !with_icraf)
  }

  if ("sensor_unseen" %in% names(fold_summary)) {
    unseen <- !is.na(fold_summary$sensor_unseen) & fold_summary$sensor_unseen
    seen   <- !is.na(fold_summary$sensor_unseen) & !fold_summary$sensor_unseen
    add_group("sensor_status", "sensor_unseen", unseen)
    add_group("sensor_status", "sensor_seen", seen)
  }

  if (length(rows) < 1L) {
    return(data.frame(
      grouping = character(0),
      group = character(0),
      n_fold_results = integer(0),
      n_unique_datasets = integer(0),
      RMSE_mean = numeric(0),
      RMSE_sd = numeric(0),
      Bias_mean = numeric(0),
      Bias_sd = numeric(0),
      R2_mean = numeric(0),
      R2_sd = numeric(0),
      RPIQ_mean = numeric(0),
      RPIQ_sd = numeric(0),
      CCC_mean = numeric(0),
      CCC_sd = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.build_dataset_summary <- function(fold_summary, strategy = "lodo_dataset") {
  if (!identical(strategy, "lodo_dataset")) {
    return(data.frame(
      held_out_dataset = character(0),
      held_out_sensor = character(0),
      sensor_unseen = logical(0),
      n_fold_results = integer(0),
      n_eval_rows = integer(0),
      n_holdout_rows = integer(0),
      n_unique_datasets = integer(0),
      RMSE_mean = numeric(0),
      RMSE_sd = numeric(0),
      Bias_mean = numeric(0),
      Bias_sd = numeric(0),
      R2_mean = numeric(0),
      R2_sd = numeric(0),
      RPIQ_mean = numeric(0),
      RPIQ_sd = numeric(0),
      CCC_mean = numeric(0),
      CCC_sd = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  dataset_rows <- split(fold_summary, as.character(fold_summary$held_out_dataset))
  rows <- lapply(names(dataset_rows), function(dataset_name) {
    ds_rows <- dataset_rows[[dataset_name]]
    sm <- .summarise_fold_metric_rows(ds_rows)
    unseen_vals <- ds_rows$sensor_unseen
    unseen_vals <- unseen_vals[!is.na(unseen_vals)]
    sensor_unseen <- if (length(unseen_vals) < 1L) {
      NA
    } else if (all(unseen_vals)) {
      TRUE
    } else if (all(!unseen_vals)) {
      FALSE
    } else {
      NA
    }
    data.frame(
      held_out_dataset = dataset_name,
      held_out_sensor  = if ("held_out_sensor" %in% names(ds_rows))
        .collapse_unique_labels(ds_rows$held_out_sensor) else NA_character_,
      sensor_unseen    = sensor_unseen,
      n_fold_results   = sm$n_fold_results,
      n_eval_rows      = sum(ds_rows$n_test, na.rm = TRUE),
      n_holdout_rows   = if ("n_holdout" %in% names(ds_rows))
        sum(ds_rows$n_holdout, na.rm = TRUE) else sum(ds_rows$n_test, na.rm = TRUE),
      n_unique_datasets = sm$n_unique_datasets,
      RMSE_mean        = sm$RMSE_mean,
      RMSE_sd          = sm$RMSE_sd,
      Bias_mean        = sm$Bias_mean,
      Bias_sd          = sm$Bias_sd,
      R2_mean          = sm$R2_mean,
      R2_sd            = sm$R2_sd,
      RPIQ_mean        = sm$RPIQ_mean,
      RPIQ_sd          = sm$RPIQ_sd,
      CCC_mean         = sm$CCC_mean,
      CCC_sd           = sm$CCC_sd,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$held_out_dataset), , drop = FALSE]
}

.summarise_cv_folds <- function(fold_summary) {
  out <- list(
    RMSE_mean = mean(fold_summary$RMSE, na.rm = TRUE),
    RMSE_sd   = stats::sd(fold_summary$RMSE, na.rm = TRUE),
    Bias_mean = mean(fold_summary$Bias, na.rm = TRUE),
    Bias_sd   = stats::sd(fold_summary$Bias, na.rm = TRUE),
    R2_mean   = mean(fold_summary$R2,   na.rm = TRUE),
    R2_sd     = stats::sd(fold_summary$R2,   na.rm = TRUE),
    RPIQ_mean = mean(fold_summary$RPIQ, na.rm = TRUE),
    RPIQ_sd   = stats::sd(fold_summary$RPIQ, na.rm = TRUE),
    CCC_mean  = mean(fold_summary$CCC,  na.rm = TRUE),
    CCC_sd    = stats::sd(fold_summary$CCC,  na.rm = TRUE),
    n_folds   = nrow(fold_summary),
    fold_RMSE = fold_summary$RMSE,
    fold_Bias = fold_summary$Bias,
    fold_R2   = fold_summary$R2,
    fold_RPIQ = fold_summary$RPIQ,
    fold_CCC  = fold_summary$CCC,
    fold_dataset = fold_summary$held_out_dataset,
    fold_n_test  = fold_summary$n_test
  )
  if ("n_holdout" %in% names(fold_summary))
    out$fold_n_holdout <- fold_summary$n_holdout
  if ("held_out_sensor" %in% names(fold_summary))
    out$fold_sensor <- fold_summary$held_out_sensor
  if ("train_sensors" %in% names(fold_summary))
    out$fold_train_sensors <- fold_summary$train_sensors
  if ("sensor_unseen" %in% names(fold_summary))
    out$fold_sensor_unseen <- fold_summary$sensor_unseen
  out
}

.aggregate_cv_runs <- function(cv_runs) {
  if (length(cv_runs) < 1L)
    stop("cv_runs must contain at least one repeat result.")

  repeat_summary <- do.call(
    rbind,
    lapply(cv_runs, function(x) x$repeat_summary)
  )
  fold_summary <- do.call(
    rbind,
    lapply(cv_runs, function(x) x$fold_summary)
  )
  rownames(repeat_summary) <- NULL
  rownames(fold_summary) <- NULL

  first_summary <- cv_runs[[1L]]$summary
  first_strategy <- cv_runs[[1L]]$strategy %||% "lodo_dataset"
  first_strategy_detail <- cv_runs[[1L]]$strategy_detail %||%
    "LODO — stratified by OSSL contributing dataset"
  metric_names <- c("RMSE", "Bias", "R2", "RPIQ", "CCC")

  out <- list(
    n_folds        = first_summary$n_folds,
    n_repeats      = nrow(repeat_summary),
    strategy       = first_strategy,
    strategy_detail = first_strategy_detail,
    repeat_summary = repeat_summary,
    fold_summary   = fold_summary,
    dataset_summary = .build_dataset_summary(fold_summary, strategy = first_strategy),
    fold_group_summary = .build_fold_group_summary(fold_summary, strategy = first_strategy),
    fold_dataset   = fold_summary$held_out_dataset,
    fold_n_test    = fold_summary$n_test,
    fold_n_holdout = if ("n_holdout" %in% names(fold_summary)) fold_summary$n_holdout else fold_summary$n_test,
    fold_RMSE      = fold_summary$RMSE,
    fold_Bias      = fold_summary$Bias,
    fold_R2        = fold_summary$R2,
    fold_RPIQ      = fold_summary$RPIQ,
    fold_CCC       = fold_summary$CCC
  )

  for (metric in metric_names) {
    repeat_col <- paste0(metric, "_mean")
    repeat_vals <- repeat_summary[[repeat_col]]
    fold_vals <- fold_summary[[metric]]
    out[[repeat_col]] <- mean(repeat_vals, na.rm = TRUE)
    out[[paste0(metric, "_sd")]] <- if (nrow(repeat_summary) > 1L) {
      stats::sd(repeat_vals, na.rm = TRUE)
    } else {
      first_summary[[paste0(metric, "_sd")]]
    }
    out[[paste0(metric, "_repeat_sd")]] <- if (nrow(repeat_summary) > 1L) {
      stats::sd(repeat_vals, na.rm = TRUE)
    } else {
      NA_real_
    }
    out[[paste0(metric, "_fold_sd")]] <- stats::sd(fold_vals, na.rm = TRUE)
    out[[paste0("repeat_", metric)]] <- repeat_vals
  }

  out
}

.run_single_lodo_cv <- function(X, y, dataset_col,
                                sample_id = NULL,
                                sensor_labels = NULL,
                                cv_eval_mask = NULL,
                                k,
                                epochs, batch, latent,
                                pat_es, pat_lr,
                                kl_beta, kl_warmup,
                                mask = NULL,
                                cv_strategy = "grouped_sample_kfold",
                                fold_seed = 42L,
                                train_seed = NULL,
                                repeat_id = 1L) {
  fold_info <- .resolve_cv_folds(
    dataset_col,
    y = y,
    k = k,
    seed = fold_seed,
    strategy = cv_strategy,
    group_col = sample_id
  )
  folds <- fold_info$folds
  k     <- length(folds)
  masked <- !is.null(mask)
  if (!is.null(sensor_labels)) {
    sensor_labels <- as.character(sensor_labels)
    if (length(sensor_labels) != nrow(X)) {
      stop("sensor_labels must have length equal to nrow(X).")
    }
    sensor_labels <- .normalize_ossl_sensor(sensor_labels)
  }
  if (!is.null(cv_eval_mask)) {
    cv_eval_mask <- as.logical(cv_eval_mask)
    if (length(cv_eval_mask) != nrow(X)) {
      stop("cv_eval_mask must have length equal to nrow(X).")
    }
    cv_eval_mask[is.na(cv_eval_mask)] <- FALSE
  }

  fold_summary <- data.frame(
    repeat_id        = rep.int(as.integer(repeat_id), k),
    fold             = seq_len(k),
    fold_seed        = rep.int(as.integer(fold_seed), k),
    train_seed       = rep.int(
      if (is.null(train_seed)) NA_integer_ else as.integer(train_seed),
      k
    ),
    held_out_dataset = rep(NA_character_, k),
    held_out_sensor  = rep(NA_character_, k),
    train_sensors    = rep(NA_character_, k),
    sensor_unseen    = rep(NA, k),
    n_holdout        = integer(k),
    n_test           = integer(k),
    RMSE             = rep(NA_real_, k),
    Bias             = rep(NA_real_, k),
    R2               = rep(NA_real_, k),
    RPIQ             = rep(NA_real_, k),
    CCC              = rep(NA_real_, k),
    stringsAsFactors = FALSE
  )

  mean_y <- mean(y, na.rm = TRUE)
  sd_y   <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z <- (y - mean_y) / sd_y

  for (fold in seq_len(k)) {
    te_idx <- folds[[fold]]
    te_eval_idx <- if (is.null(cv_eval_mask)) {
      te_idx
    } else {
      te_idx[cv_eval_mask[te_idx]]
    }
    fold_summary$held_out_dataset[fold] <- if (identical(fold_info$strategy, "grouped_sample_kfold")) {
      paste0("Fold ", fold)
    } else {
      paste(
        sort(unique(as.character(dataset_col[te_idx]))),
        collapse = " + "
      )
    }
    fold_summary$n_holdout[fold] <- length(te_idx)
    fold_summary$n_test[fold] <- length(te_eval_idx)

    if (length(te_idx) == 0L)
      next
    if (length(te_eval_idx) == 0L)
      next

    tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
    fold_label <- fold_summary$held_out_dataset[fold]
    message(sprintf(
      "      CV progress: repeat %d fold %d/%d [%s] started (train=%d, eval=%d)",
      repeat_id, fold, k, fold_label, length(tr_idx), length(te_eval_idx)
    ))
    if (!is.null(sensor_labels)) {
      test_sensors <- unique(sensor_labels[te_idx])
      train_sensors <- unique(sensor_labels[tr_idx])
      fold_summary$held_out_sensor[fold] <- .collapse_unique_labels(test_sensors)
      fold_summary$train_sensors[fold] <- .collapse_unique_labels(train_sensors)
      fold_summary$sensor_unseen[fold] <- length(setdiff(test_sensors, train_sensors)) > 0L
    }
    X_tr   <- X[tr_idx, , drop = FALSE]; y_tr <- y_z[tr_idx]
    X_te   <- X[te_eval_idx, , drop = FALSE]; y_te <- y[te_eval_idx]

    fold_train_seed <- if (is.null(train_seed)) {
      NULL
    } else {
      as.integer(train_seed + fold - 1L)
    }
    .autospectra_set_random_seed(fold_train_seed)

    if (masked) {
      mask_tr <- mask[tr_idx, , drop = FALSE]
      mask_te <- mask[te_eval_idx, , drop = FALSE]
      X_in_tr <- cbind(X_tr, mask_tr)
      X_in_te <- cbind(X_te, mask_te)
      mdl_fold <- build_maskedVAE(ncol(X_tr), latent_dim = latent)
    } else {
      X_in_tr  <- X_tr
      X_in_te  <- X_te
      mdl_fold <- build_soilVAE(ncol(X_tr), latent_dim = latent)
    }

    cbs <- list(
      keras3::callback_early_stopping(
        monitor = "val_loss", patience = pat_es,
        restore_best_weights = TRUE),
      keras3::callback_reduce_lr_on_plateau(
        monitor = "val_loss", patience = pat_lr,
        factor = 0.5, min_lr = 1e-5),
      .vae_env$KLAnnealingCallback(
        target_beta   = kl_beta,
        warmup_epochs = as.integer(kl_warmup))
    )
    suppressMessages(
      .keras_fit_quiet(
        mdl_fold,
        x = as.matrix(X_in_tr),
        y_recon = as.matrix(X_tr),
        y_pred = y_tr,
        epochs = epochs, batch_size = batch,
        validation_split = 0.15,
        callbacks = cbs, verbose = 0
      )
    )
    yhat_z <- .extract_prediction(mdl_fold, X_in_te)
    yhat   <- yhat_z * sd_y + mean_y
    mets   <- metrics_from_y(y_te, yhat)
    fold_summary$RMSE[fold] <- mets$RMSE
    fold_summary$Bias[fold] <- mets$Bias
    fold_summary$R2[fold]   <- mets$R2
    fold_summary$RPIQ[fold] <- mets$RPIQ
    fold_summary$CCC[fold]  <- mets$CCC
    message(sprintf(
      "      CV progress: repeat %d fold %d/%d [%s] done (R2=%.3f, RMSE=%.3f)",
      repeat_id, fold, k, fold_label, mets$R2, mets$RMSE
    ))
    keras3::clear_session()
    rm(mdl_fold); gc(verbose = FALSE)
  }

  repeat_metrics <- .summarise_cv_folds(fold_summary)
  repeat_row <- data.frame(
    repeat_id   = as.integer(repeat_id),
    fold_seed   = as.integer(fold_seed),
    train_seed  = if (is.null(train_seed)) NA_integer_ else as.integer(train_seed),
    n_folds     = repeat_metrics$n_folds,
    RMSE_mean   = repeat_metrics$RMSE_mean,
    Bias_mean   = repeat_metrics$Bias_mean,
    R2_mean     = repeat_metrics$R2_mean,
    RPIQ_mean   = repeat_metrics$RPIQ_mean,
    CCC_mean    = repeat_metrics$CCC_mean,
    stringsAsFactors = FALSE
  )

  list(
    summary        = repeat_metrics,
    fold_summary   = fold_summary,
    repeat_summary = repeat_row,
    strategy       = fold_info$strategy,
    strategy_detail = fold_info$strategy_detail
  )
}

# ---- Single-property CV training ---------------------------------------------

#' Train and cross-validate a soilVAE (or maskedVAE) for one soil property
#'
#' Runs grouped or LODO cross-validation depending on \code{cv_strategy}, then
#' optionally trains a final model on all data (80/10/10 grouped split) and saves
#' weights, scaler, conformal quantiles, and latent statistics to disk.
#'
#' @param X           Numeric matrix (rows = samples, cols = preprocessed bands).
#'   For masked VAE this is the SNV-normalised spectrum with zeros for absent bands.
#' @param y           Numeric vector of raw soil-property values (length = nrow(X)).
#' @param dataset_col Character vector of OSSL dataset labels (length = nrow(X)).
#' @param family_id   Family identifier string used for the output sub-directory.
#' @param prop        Soil property name (used for file names).
#' @param k           Maximum CV folds (default 10; reduced automatically if fewer groups).
#' @param epochs      Max training epochs (default 80).
#' @param batch       Mini-batch size (default 32).
#' @param latent      Latent space dimension (default 16).
#' @param pat_es      Early-stopping patience in epochs (default 10).
#' @param pat_lr      LR-reduction patience in epochs (default 5).
#' @param kl_beta     Target KL weight after annealing (default 4e-5).
#' @param kl_warmup   Epochs to ramp beta from 0 to \code{kl_beta} (default 30).
#' @param cv_repeats  Number of repeated CV runs to average (default 1).
#' @param seed        Optional integer seed used for fold assignment and model
#'   initialisation. \code{NULL} keeps the legacy stochastic behaviour.
#' @param min_n       Minimum samples required to attempt training (default 50).
#' @param out_dir     Root output directory (default "models").
#' @param cv_only     If \code{TRUE}, skip final-model training and return only
#'   CV metrics (used in hyperparameter grid search).
#' @param mask        Binary mask matrix, same shape as \code{X} (1 = valid band,
#'   0 = absent).  \code{NULL} (default) uses the standard soilVAE; a matrix
#'   activates the Masked VAE (\code{build_maskedVAE}).
#' @param sample_id   Optional grouping vector (typically \code{Soil_ID}) used
#'   to keep replicate spectra together across train/calibration/test splits.
#' @param sensor_labels Optional character vector of sensor labels (length = nrow(X)).
#' @param cv_eval_mask Optional logical vector (length = nrow(X)) marking which
#'   rows are eligible for CV evaluation. Training still excludes all rows from
#'   the held-out dataset fold; this mask only controls which held-out rows are
#'   scored. Useful for masked VAE training with augmentation but evaluation on
#'   the original full-spectrum copy only.
#' @param cv_eval_mode Character label saved with the CV output to describe the
#'   evaluation regime (for example \code{"full_spectrum_original"} or
#'   \code{"stress_test"}).
#' @param cv_strategy Cross-validation regime. Default
#'   \code{"grouped_sample_kfold"} uses grouped stratified K-fold by
#'   \code{sample_id}; \code{"lodo_dataset"} keeps the harsher
#'   leave-one-dataset-out benchmark.
#' @return Named list with elements \code{cv}, \code{test}, and \code{n}, or
#'   \code{NULL} if skipped.
#' @export
train_property_cv <- function(X, y, dataset_col, family_id, prop,
                              k         = 10L,
                              epochs    = 80L,
                               batch     = 32L,
                               latent    = 16L,
                               pat_es    = 10L,
                               pat_lr    = 5L,
                               kl_beta   = 1e-3,
                               kl_warmup = 30L,
                               cv_repeats = 1L,
                               seed      = NULL,
                               min_n     = 50L,
                              out_dir   = "models",
                              cv_only   = FALSE,
                              mask      = NULL,
                              sample_id = NULL,
                              sensor_labels = NULL,
                              cv_eval_mask = NULL,
                              cv_eval_mode = NULL,
                              cv_strategy = "grouped_sample_kfold") {

  n <- nrow(X)
  if (n < min_n) {
    message("    SKIP: need >= ", min_n, ", have ", n)
    return(NULL)
  }
  if (is.null(sample_id))
    sample_id <- seq_len(n)

  cv_repeats <- max(1L, suppressWarnings(as.integer(cv_repeats)[1L]))
  base_seed  <- .autospectra_seed_value(seed)
  masked     <- !is.null(mask)
  repeat_fold_seeds <- if (is.null(base_seed)) {
    42L + seq_len(cv_repeats) - 1L
  } else {
    base_seed + ((seq_len(cv_repeats) - 1L) * 1000L)
  }
  repeat_train_seeds <- if (is.null(base_seed)) {
    rep(NA_integer_, cv_repeats)
  } else {
    repeat_fold_seeds
  }

  cv_runs <- lapply(seq_len(cv_repeats), function(rep_id) {
    .run_single_lodo_cv(
      X = X,
      y = y,
      dataset_col = dataset_col,
      sample_id = sample_id,
      sensor_labels = sensor_labels,
      cv_eval_mask = cv_eval_mask,
      k = k,
      epochs = epochs,
      batch = batch,
      latent = latent,
      pat_es = pat_es,
      pat_lr = pat_lr,
      kl_beta = kl_beta,
      kl_warmup = kl_warmup,
      mask = mask,
      cv_strategy = cv_strategy,
      fold_seed = repeat_fold_seeds[[rep_id]],
      train_seed = if (is.na(repeat_train_seeds[[rep_id]])) NULL else repeat_train_seeds[[rep_id]],
      repeat_id = rep_id
    )
  })
  cv <- .aggregate_cv_runs(cv_runs)
  cv$eval_mode <- cv_eval_mode %||% if (!is.null(cv_eval_mask)) {
    "full_spectrum_original"
  } else if (masked) {
    "stress_test"
  } else {
    "standard"
  }

  if (cv_only) return(list(cv = cv, n = n))

  # ---- Final model: train on ALL data, 80/10/10 split ----------------------
  split_seed <- if (is.null(base_seed)) 42L else as.integer(base_seed + 900000L)
  sp <- .grouped_split_idx(
    group_col = sample_id,
    y = y,
    seed = split_seed,
    p_train = 0.8,
    p_cal = 0.1,
    p_test = 0.1
  )
  message(sprintf(
    "      Final fit: train=%d, calib=%d, test=%d",
    length(sp$train), length(sp$calib), length(sp$test)
  ))
  mean_y <- mean(y, na.rm = TRUE)
  sd_y   <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z <- (y - mean_y) / sd_y
  X_tr <- X[sp$train, , drop = FALSE]; y_tr <- y_z[sp$train]
  X_ca <- X[sp$calib, , drop = FALSE]; y_ca <- y[sp$calib]
  X_te <- X[sp$test,  , drop = FALSE]; y_te <- y[sp$test]

  .autospectra_set_random_seed(if (is.null(base_seed)) NULL else base_seed + 100000L)
  if (masked) {
    X_in_tr_f <- cbind(X_tr, mask[sp$train, , drop = FALSE])
    X_in_ca_f <- cbind(X_ca, mask[sp$calib, , drop = FALSE])
    X_in_te_f <- cbind(X_te, mask[sp$test,  , drop = FALSE])
    mdl_final <- build_maskedVAE(ncol(X_tr), latent_dim = latent)
  } else {
    X_in_tr_f <- X_tr
    X_in_ca_f <- X_ca
    X_in_te_f <- X_te
    mdl_final <- build_soilVAE(ncol(X_tr), latent_dim = latent)
  }

  cbs_final <- list(
    keras3::callback_early_stopping(
      monitor = "val_loss", patience = pat_es,
      restore_best_weights = TRUE),
    keras3::callback_reduce_lr_on_plateau(
      monitor = "val_loss", patience = pat_lr,
      factor = 0.5, min_lr = 1e-5),
    .vae_env$KLAnnealingCallback(
      target_beta   = kl_beta,
      warmup_epochs = as.integer(kl_warmup))
  )
  suppressMessages(
    .keras_fit_quiet(
      mdl_final,
      x = as.matrix(X_in_tr_f),
      y_recon = as.matrix(X_tr),
      y_pred = y_tr,
      epochs = epochs, batch_size = batch,
      validation_split = 0.15,
      callbacks = cbs_final, verbose = 0
    )
  )

  # Save model + scaler
  model_dir_prop <- file.path(out_dir, family_id, "models")
  dir_create(model_dir_prop)
  base <- file.path(model_dir_prop, prop)
  keras3::save_model(mdl_final, .default_model_path(base), overwrite = TRUE)
  save_scaler(mean_y, sd_y, path_base = base)

  # Conformal calibration (absolute residuals on calib set)
  yhat_ca <- .extract_prediction(mdl_final, X_in_ca_f) * sd_y + mean_y
  res_ca  <- abs(y_ca - yhat_ca)
  q90 <- as.numeric(stats::quantile(res_ca, 0.90, na.rm = TRUE))
  q95 <- as.numeric(stats::quantile(res_ca, 0.95, na.rm = TRUE))

  # Test-set metrics (held-out 10%)
  yhat_te <- .extract_prediction(mdl_final, X_in_te_f) * sd_y + mean_y
  mets_te <- metrics_from_y(y_te, yhat_te)

  # Latent statistics for Mahalanobis applicability domain
  encoder <- keras3::keras_model(inputs  = mdl_final$input,
                                 outputs = mdl_final$get_layer("z_mean")$output)
  Z_tr  <- .keras_predict_quiet(encoder, X_in_tr_f, verbose = 0)
  mu_z  <- colMeans(Z_tr, na.rm = TRUE)
  Sig_z <- stats::cov(Z_tr, use = "pairwise.complete.obs")
  Sig_z <- as.matrix(Sig_z + diag(1e-6, ncol(Sig_z)))
  d_lat <- length(mu_z)

  latent_stats <- list(
    mu    = as.numeric(mu_z),
    Sigma = unname(asplit(Sig_z, 1L)),
    df    = d_lat,
    thr95 = stats::qchisq(0.95, df = d_lat)
  )

  save_metrics_json(base,
                    metrics    = mets_te,
                    conf_q     = list(q90 = q90, q95 = q95),
                    feat_range = feature_minmax(X),
                    latent     = latent_stats)

  saveRDS(list(cv = cv, test = mets_te, n = n),
          paste0(base, "_cv.rds"))

  message(sprintf(
    "      Final fit done: test R2=%.3f, RMSE=%.3f",
    mets_te$R2, mets_te$RMSE
  ))

  list(cv = cv, test = mets_te, n = n)
}

