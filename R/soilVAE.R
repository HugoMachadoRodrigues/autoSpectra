# R/soilVAE.R — soilVAE true Variational Autoencoder (architecture + training)
#
# Architecture (encoder → latent [VAE] → two heads):
#   Input (d_in)
#     Encoder: Dense(256)+BN → Dense(128)+BN → Dense(64)+BN+Dropout → Dense(32)+BN
#     Latent:  z_mean    = Dense(latent_dim, name="z_mean")
#              z_log_var = Dense(latent_dim, name="z_log_var")
#              z         = Sampling()([z_mean, z_log_var])  ← reparam. trick + KL
#     Reconstruction head: Dense(32)+BN → Dense(64)+BN → Dense(128)+BN → Dense(d_in)
#     Prediction head:     Dense(128)+Dropout → Dense(64) → Dense(32) → Dense(1)
#
# KL divergence: injected as add_loss() inside the Sampling Python layer.
# KL annealing:  KLAnnealingCallback ramps beta 0 → target over warmup_epochs,
#                preventing posterior collapse before reconstruction is learned.
# BatchNorm:     on all encoder/decoder Dense layers for OSSL-scale stability.
# L2 reg:        on all Dense layers to prevent overfitting.
# z_mean:        (deterministic) used for Mahalanobis applicability-domain checks.
#
# Loss = MSE(reconstruction) * w1  +  MSE(prediction) * w2  +  beta * KL
# A scaler (z-score on y) is saved alongside each model.
# Conformal calibration quantiles (q90, q95) are saved for uncertainty.
# Latent space statistics (μ, Σ) enable applicability domain checks.

# ---- Python VAE classes (defined once per R session) ------------------------

.vae_env <- new.env(parent = emptyenv())
.vae_env$initialized <- FALSE

#' Initialise VAE Python classes (idempotent)
#'
#' Defines \code{Sampling} and \code{KLAnnealingCallback} in Python via
#' \pkg{reticulate} and stores references in \code{.vae_env}.  Called
#' automatically by \code{build_soilVAE()} and \code{load_soilVAE()}.
#' @keywords internal
.init_vae_python <- function() {
  if (isTRUE(.vae_env$initialized)) return(invisible(NULL))

  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required for the true VAE. ",
         "Install with: install.packages('reticulate')")
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")

  reticulate::py_run_string("
import builtins
import keras
from keras import layers, ops, random as k_random

class Sampling(layers.Layer):
    '''Reparametrisation trick with scaled KL regularisation.

    z = z_mean + exp(0.5 * z_log_var) * epsilon,  epsilon ~ N(0, I)

    The KL divergence is added to the model loss scaled by beta.
    Beta is read from builtins._VAE_BETA, which KLAnnealingCallback
    ramps from 0 to the target value over the warmup period.
    '''
    def call(self, inputs):
        z_mean, z_log_var = inputs
        z_log_var = ops.clip(z_log_var, -10.0, 6.0)
        eps = k_random.normal(shape=ops.shape(z_mean))
        z   = z_mean + ops.exp(0.5 * z_log_var) * eps
        kl  = -0.5 * ops.sum(
            1.0 + z_log_var - ops.square(z_mean) - ops.exp(z_log_var),
            axis=-1
        )
        beta = float(getattr(builtins, '_VAE_BETA', 0.0))
        self.add_loss(beta * ops.mean(kl))
        return z

    def get_config(self):
        return super().get_config()


class KLAnnealingCallback(keras.callbacks.Callback):
    '''Linearly ramps builtins._VAE_BETA from 0 to target_beta.

    Annealing prevents the encoder from being forced into N(0, I)
    before reconstruction loss is learned (posterior collapse).
    '''
    def __init__(self, target_beta, warmup_epochs=30, **kwargs):
        super().__init__(**kwargs)
        self.target_beta   = float(target_beta)
        self.warmup_epochs = int(warmup_epochs)

    def on_train_begin(self, logs=None):
        builtins._VAE_BETA = 0.0

    def on_epoch_begin(self, epoch, logs=None):
        if self.warmup_epochs > 0 and epoch < self.warmup_epochs:
            builtins._VAE_BETA = (
                self.target_beta * (epoch + 1) / self.warmup_epochs
            )
        else:
            builtins._VAE_BETA = self.target_beta
")

  py_main <- reticulate::import_main()
  .vae_env$Sampling            <- py_main$Sampling
  .vae_env$KLAnnealingCallback <- py_main$KLAnnealingCallback
  .vae_env$initialized         <- TRUE
  invisible(NULL)
}

# ---- Model builder ----------------------------------------------------------

#' Build the soilVAE true Variational Autoencoder
#'
#' Requires the \pkg{keras3} and \pkg{reticulate} packages (Keras 3.x).
#'
#' @param d_in Integer; number of input features (preprocessed spectral bands)
#' @param latent_dim Integer; latent space dimension (default 16)
#' @param loss_weights Numeric vector length 2; weights for reconstruction and
#'   prediction losses (default c(0.3, 0.3))
#' @param l2_reg L2 regularisation coefficient on all Dense layers (default 1e-5)
#' @param dropout_rate Dropout rate for the prediction head (default 0.2)
#' @param encoder_dropout Dropout rate inside the encoder (default 0.05)
#' @param learning_rate Initial Adam learning rate (default 1e-3)
#' @return A compiled keras Model object
#' @export
build_soilVAE <- function(d_in,
                          latent_dim      = 16L,
                          loss_weights    = c(0.3, 0.3),
                          l2_reg          = 1e-5,
                          dropout_rate    = 0.2,
                          encoder_dropout = 0.05,
                          learning_rate   = 1e-3) {

  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required. Install with: install.packages('keras3')")

  .init_vae_python()

  reg <- keras3::regularizer_l2(l2 = l2_reg)

  inp <- keras3::layer_input(shape = d_in)

  # ----- Encoder (Dense + BatchNorm pyramid) --------------------------------
  enc <- inp |>
    keras3::layer_dense(256, activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(128, activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(64,  activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dropout(encoder_dropout) |>
    keras3::layer_dense(32,  activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization()

  # ----- Latent space (reparametrisation trick) -----------------------------
  z_mean    <- keras3::layer_dense(enc, latent_dim, name = "z_mean")
  z_log_var <- keras3::layer_dense(enc, latent_dim, name = "z_log_var")
  z         <- .vae_env$Sampling()(list(z_mean, z_log_var))

  # ----- Reconstruction head (decoder pyramid) ------------------------------
  rec <- z |>
    keras3::layer_dense(32,   activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(64,   activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(128,  activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(d_in, activation = "linear", name = "reconstruction")

  # ----- Prediction head ----------------------------------------------------
  pred <- z |>
    keras3::layer_dense(128, activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_dropout(dropout_rate) |>
    keras3::layer_dense(64,  activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_dense(32,  activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_dense(1,   activation = "linear", name = "prediction")

  mdl <- keras3::keras_model(inputs = inp, outputs = list(rec, pred))
  mdl |> keras3::compile(
    optimizer    = keras3::optimizer_adam(learning_rate = learning_rate),
    loss         = list("mse", "mse"),
    loss_weights = as.list(loss_weights)
  )
  mdl
}

# ---- Single-property trainer ------------------------------------------------

#' Train a soilVAE model for a single soil property
#'
#' Performs train / calibration / test split, fits the model with early
#' stopping, learning-rate reduction and KL annealing, computes test metrics
#' and conformal quantiles, and saves the model + scaler + metrics to disk.
#'
#' @param X Numeric matrix (rows = samples, cols = preprocessed bands)
#' @param y Numeric vector of soil property values (length == nrow(X))
#' @param family_id Family identifier string (used for output directory)
#' @param prop Soil property name (used for output file names)
#' @param out_dir Root output directory (default: "models")
#' @param epochs Max training epochs (default 30)
#' @param batch_size Mini-batch size (default 32)
#' @param val_split Validation fraction from training data (default 0.2)
#' @param patience_es Early stopping patience (default 8)
#' @param patience_rl LR reduction patience (default 4)
#' @param latent_dim Latent dimension (default 16)
#' @param kl_beta Target KL weight after annealing (default 4e-5)
#' @param kl_warmup_epochs Epochs to ramp beta from 0 to kl_beta (default 30)
#' @param min_n Minimum samples required to attempt training (default 30)
#' @return Invisible TRUE on success, NULL if skipped
#' @export
train_soilVAE <- function(X, y, family_id, prop,
                          out_dir          = "models",
                          epochs           = 30L,
                          batch_size       = 32L,
                          val_split        = 0.2,
                          patience_es      = 8L,
                          patience_rl      = 4L,
                          latent_dim       = 16L,
                          kl_beta          = 4e-5,
                          kl_warmup_epochs = 30L,
                          min_n            = 30L) {

  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")

  n <- nrow(X)
  if (n < min_n) {
    message("Skip ", family_id, "/", prop, ": need >=", min_n, ", have ", n)
    return(invisible(NULL))
  }

  # Z-score target
  mean_y <- mean(y, na.rm = TRUE)
  sd_y   <- sd(y, na.rm = TRUE)
  if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z <- (y - mean_y) / sd_y

  # Train / calib / test split
  sp    <- split_idx(n, seed = 42)
  X_tr  <- X[sp$train, , drop = FALSE]; y_tr_z <- y_z[sp$train]
  X_ca  <- X[sp$calib, , drop = FALSE]; y_ca   <- y[sp$calib]
  X_te  <- X[sp$test,  , drop = FALSE]; y_te   <- y[sp$test]

  # Build and train
  mdl <- build_soilVAE(ncol(X_tr), latent_dim = latent_dim)
  cbs <- list(
    keras3::callback_early_stopping(
      monitor = "val_loss", patience = patience_es,
      restore_best_weights = TRUE),
    keras3::callback_reduce_lr_on_plateau(
      monitor = "val_loss", patience = patience_rl,
      factor = 0.5, min_lr = 1e-5),
    .vae_env$KLAnnealingCallback(
      target_beta    = kl_beta,
      warmup_epochs  = as.integer(kl_warmup_epochs))
  )
  mdl |> keras3::fit(
    x                = as.matrix(X_tr),
    y                = list(as.matrix(X_tr), y_tr_z),
    epochs           = epochs,
    batch_size       = batch_size,
    validation_split = val_split,
    callbacks        = cbs,
    verbose          = 0
  )

  # --- Save model + scaler ---
  model_dir <- file.path(out_dir, family_id, "models")
  dir_create(model_dir)
  base <- file.path(model_dir, prop)
  keras3::save_model(mdl, paste0(base, ".h5"))
  save_scaler(mean_y, sd_y, path_base = base)

  # --- Conformal calibration (absolute residuals on calib set) ---
  yhat_ca <- .extract_prediction(mdl, X_ca) * sd_y + mean_y
  res_ca  <- abs(y_ca - yhat_ca)
  q90 <- as.numeric(quantile(res_ca, 0.90, na.rm = TRUE))
  q95 <- as.numeric(quantile(res_ca, 0.95, na.rm = TRUE))

  # --- Test metrics ---
  yhat_te <- .extract_prediction(mdl, X_te) * sd_y + mean_y
  mets    <- metrics_from_y(y_te, yhat_te)

  # --- Latent statistics for applicability domain (use z_mean — deterministic) ---
  encoder <- keras3::keras_model(inputs  = mdl$input,
                                 outputs = mdl$get_layer("z_mean")$output)
  Z_tr    <- predict(encoder, as.matrix(X_tr), verbose = 0)
  mu_z    <- colMeans(Z_tr, na.rm = TRUE)
  Sig_z   <- stats::cov(Z_tr, use = "pairwise.complete.obs")
  Sig_z   <- as.matrix(Sig_z + diag(1e-6, ncol(Sig_z)))
  d_lat   <- length(mu_z)
  thr95   <- stats::qchisq(0.95, df = d_lat)

  latent <- list(
    mu    = as.numeric(mu_z),
    Sigma = unname(asplit(Sig_z, 1L)),
    df    = d_lat,
    thr95 = thr95
  )

  # --- Feature range (post-preprocess) ---
  fr <- feature_minmax(X)

  save_metrics_json(base,
                    metrics    = mets,
                    conf_q     = list(q90 = q90, q95 = q95),
                    feat_range = fr,
                    latent     = latent)

  message(sprintf("  [%s/%s] RMSE=%.3f  Bias=%.3f  R2=%.3f  RPIQ=%.2f  CCC=%.3f  n=%d",
                  family_id, prop,
                  mets$RMSE, mets$Bias, mets$R2, mets$RPIQ, mets$CCC, n))
  invisible(TRUE)
}

# ---- Internal: extract the prediction head output from a keras model --------
.extract_prediction <- function(mdl, X) {
  raw <- predict(mdl, as.matrix(X), verbose = 0)
  yhat <- if (is.list(raw)) {
    if (!is.null(names(raw)) && "prediction" %in% names(raw))
      raw[["prediction"]]
    else
      raw[[length(raw)]]
  } else raw
  as.numeric(yhat)
}

# ---- Multi-property OSSL trainer --------------------------------------------

#' Train soilVAE models for all properties of one model family using OSSL data
#'
#' This is the main training entry point. It downloads OSSL data (if needed),
#' applies the family's preprocessing pipeline, and trains one model per
#' soil property, saving outputs to \code{out_dir/<family_id>/models/}.
#'
#' @param family_id Family ID from model_registry, e.g. "OSSL_VisNIR"
#' @param out_dir Root output directory (default: "models")
#' @param cache_dir OSSL data cache directory
#' @param properties Properties to train (default: all in family)
#' @param download_if_missing Download OSSL data if not cached
#' @param ... Additional arguments passed to train_soilVAE()
#' @return Invisible named list of results (TRUE/NULL per property)
#' @export
train_ossl_models <- function(family_id,
                              out_dir              = "models",
                              cache_dir            = ossl_cache_dir(),
                              properties           = NULL,
                              download_if_missing  = TRUE,
                              ...) {
  fam <- get_family(family_id)
  if (fam$source != "ossl")
    stop("family_id '", family_id, "' is not an OSSL family. ",
         "Use train_soilVAE() directly for local families.")

  props <- if (is.null(properties)) fam$properties else properties

  message("=== autoSpectra: Training [", fam$label, "] ===")
  message("  Sensor type : ", fam$sensor_type)
  message("  Grid        : ", length(fam$wavegrid), " channels (",
          min(fam$wavegrid), " - ", max(fam$wavegrid), ")")
  message("  Pipeline    : ", paste(fam$preprocess, collapse = " -> "))
  message("  Properties  : ", length(props))

  # Load and prepare data
  joined_df <- ossl_prepare(
    sensor_type         = fam$sensor_type,
    cache_dir           = cache_dir,
    properties          = props,
    download_if_missing = download_if_missing
  )

  # Extract and resample spectra to family wavegrid
  wl_info <- get_wavelengths(joined_df, id_col = "Soil_ID")
  X_src   <- as.matrix(joined_df[, wl_info$cols, drop = FALSE])
  X_res   <- resample_to_grid(X_src, src_wl = wl_info$wl,
                              target_wl = fam$wavegrid)

  # Apply preprocessing
  message("  Preprocessing ...")
  X_proc <- apply_pipeline(X_res, fam$preprocess)

  # Train one model per property
  results <- list()
  for (prop in props) {
    y <- suppressWarnings(as.numeric(joined_df[[prop]]))
    keep <- is.finite(y) & is.finite(rowSums(X_proc))
    if (sum(keep) < 30) {
      message("  Skip [", prop, "]: only ", sum(keep), " complete rows")
      results[[prop]] <- NULL
      next
    }
    message("  Training [", prop, "] (n=", sum(keep), ") ...")
    results[[prop]] <- try(
      train_soilVAE(X_proc[keep, ], y[keep], family_id, prop,
                    out_dir = out_dir, ...),
      silent = FALSE
    )
  }
  message("=== Done: ", family_id, " ===")
  invisible(results)
}

# ---- Model loader -----------------------------------------------------------

#' Load a saved soilVAE model and its scaler
#'
#' @param family_id Family identifier
#' @param prop Soil property name
#' @param model_dir Root model directory (default "models")
#' @return List with \code{model} (keras) and \code{scaler} (list with mean, sd)
#' @export
load_soilVAE <- function(family_id, prop, model_dir = "models") {
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")

  .init_vae_python()

  base <- file.path(model_dir, family_id, "models", prop)
  h5   <- paste0(base, ".h5")
  if (!file.exists(h5)) stop("Model file not found: ", h5)

  mdl <- keras3::load_model(
    h5,
    custom_objects = list("Sampling" = .vae_env$Sampling)
  )

  scaler <- if (file.exists(paste0(base, "_scaler.rds"))) {
    readRDS(paste0(base, "_scaler.rds"))
  } else {
    list(mean = 0, sd = 1, method = "none")
  }
  list(model = mdl, scaler = scaler)
}
