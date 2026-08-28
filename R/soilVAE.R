# R/soilVAE.R - soilVAE Variational Autoencoder (standard + masked)
#
# Two architectures:
#   build_soilVAE(d_in)        - Standard VAE, full-coverage spectra only
#   build_maskedVAE(d_in)      - Masked VAE; input shape = [SNV-spec | mask] (2*d_in)
#
# Both architectures share:
#   Encoder pyramid : 256 -> 128 -> 64 -> 32 (Dense + BN + Dropout)
#   Latent VAE     : z_mean, z_log_var, Sampling()  (16-D)
#   Decoder pyramid : 32 -> 64 -> 128 -> d_in (reconstruction head)
#   Prediction head : 128 -> 64 -> 32 -> 1 (regression head)
#
# Loss = w_rec * MSE(recon, target_spec) + w_pred * MSE(pred, y) + beta * KL
# beta is annealed from 0 to target via KLAnnealingCallback to prevent
# posterior collapse before reconstruction is learned.
#
# For the masked variant the reconstruction target is the SNV-normalised
# spectrum on COVERED bands only; uncovered bands are masked out of the loss.

# ---- Python VAE classes (defined once per R session) -----------------------

.vae_env <- new.env(parent = emptyenv())
.vae_env$initialized <- FALSE

#' Initialise VAE Python classes (idempotent)
#'
#' Defines the Sampling layer, KL annealing callback, and a masked-MSE
#' reconstruction loss in Python via reticulate.
#'
#' @keywords internal
.init_vae_python <- function() {
  if (isTRUE(.vae_env$initialized)) return(invisible(NULL))

  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required.")
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")

  # NOTE: Apple's tensorflow-metal 1.2.0 produces NaN-collapsed VAE
  # predictions for our custom masked_mse + Sampling layer combo (verified
  # 2026-04-26). For Apple Silicon users the safe configuration is plain
  # tensorflow on CPU (no metal plugin installed). We do not auto-disable
  # GPU here because doing so post-import has no effect; just don't install
  # tensorflow-metal until Apple ships a fix.

  reticulate::py_run_string("
import builtins
import keras
from keras import layers, ops, random as k_random


class Sampling(layers.Layer):
    '''Reparametrisation trick + KL regularisation injected via add_loss.

    z = z_mean + exp(0.5 * z_log_var) * epsilon,  epsilon ~ N(0, I)

    KL is scaled by builtins._VAE_BETA (ramped by KLAnnealingCallback).
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
    '''Linearly ramps builtins._VAE_BETA from 0 to target_beta.'''
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


def masked_mse(y_true, y_pred):
    '''MSE over covered bands only.

    y_true is expected to carry the binary mask in the LAST half of its
    feature axis. Layout: y_true[..., :d_in] = SNV spectrum target,
                          y_true[..., d_in:] = mask (1 covered, 0 missing).
    '''
    n     = ops.shape(y_true)[-1] // 2
    spec  = y_true[..., :n]
    mask  = y_true[..., n:]
    diff2 = ops.square(y_pred - spec) * mask
    denom = ops.maximum(ops.sum(mask, axis=-1), 1.0)
    return ops.sum(diff2, axis=-1) / denom
")

  py_main <- reticulate::import_main()
  .vae_env$Sampling            <- py_main$Sampling
  .vae_env$KLAnnealingCallback <- py_main$KLAnnealingCallback
  .vae_env$masked_mse          <- py_main$masked_mse
  .vae_env$initialized         <- TRUE
  invisible(NULL)
}

# ---- Encoder/decoder builders (shared) -----------------------------------

.build_encoder <- function(inp, latent_dim, l2_reg, encoder_dropout) {
  reg <- keras3::regularizer_l2(l2 = l2_reg)
  # Wider first layer: 1076-1701 spectral inputs deserve more capacity
  # than the bottleneck size before being squeezed into the latent.
  enc <- inp |>
    keras3::layer_dense(512, activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(256, activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(128, activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dropout(encoder_dropout) |>
    keras3::layer_dense(64,  activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization()
  z_mean    <- keras3::layer_dense(enc, latent_dim, name = "z_mean")
  z_log_var <- keras3::layer_dense(enc, latent_dim, name = "z_log_var")
  z         <- .vae_env$Sampling()(list(z_mean, z_log_var))
  list(z = z, z_mean = z_mean, z_log_var = z_log_var)
}

.build_decoder <- function(z, d_out, l2_reg, name = "reconstruction") {
  reg <- keras3::regularizer_l2(l2 = l2_reg)
  z |>
    keras3::layer_dense(64,    activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(128,   activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(256,   activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(d_out, activation = "linear", name = name)
}

.build_pred_head <- function(z, l2_reg, dropout_rate) {
  reg <- keras3::regularizer_l2(l2 = l2_reg)
  # Bigger prediction head: prediction is the primary task,
  # the head can afford more capacity than the decoder.
  z |>
    keras3::layer_dense(256, activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_dropout(dropout_rate) |>
    keras3::layer_dense(128, activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_dropout(dropout_rate / 2) |>
    keras3::layer_dense(64,  activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_dense(1,   activation = "linear", name = "prediction")
}

# ---- Standard soilVAE (full-coverage) ------------------------------------

#' Build a standard soilVAE
#'
#' @param d_in Integer; number of input spectral bands
#' @param latent_dim Integer; latent space dimension (default 16)
#' @param loss_weights Numeric vector length 2 (rec, pred); default c(0.3, 0.3)
#' @param l2_reg L2 regularisation on Dense layers (default 1e-5)
#' @param dropout_rate Dropout rate for the prediction head (default 0.2)
#' @param encoder_dropout Dropout inside encoder (default 0.05)
#' @param learning_rate Initial Adam learning rate (default 1e-3)
#' @return Compiled keras Model object
#' @export
build_soilVAE <- function(d_in,
                          latent_dim      = 32L,
                          loss_weights    = c(0.1, 1.0),
                          l2_reg          = 1e-5,
                          dropout_rate    = 0.2,
                          encoder_dropout = 0.05,
                          learning_rate   = 1e-3) {
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")
  .init_vae_python()

  inp  <- keras3::layer_input(shape = d_in)
  enc  <- .build_encoder(inp, latent_dim, l2_reg, encoder_dropout)
  rec  <- .build_decoder(enc$z, d_in, l2_reg, name = "reconstruction")
  pred <- .build_pred_head(enc$z, l2_reg, dropout_rate)

  mdl <- keras3::keras_model(inputs = inp, outputs = list(rec, pred))
  mdl |> keras3::compile(
    optimizer    = keras3::optimizer_adam(learning_rate = learning_rate),
    loss         = list("mse", "mse"),
    loss_weights = as.list(loss_weights)
  )
  mdl
}

# ---- Masked soilVAE (sensor-agnostic) ------------------------------------

#' Build a masked soilVAE for sensor-agnostic prediction
#'
#' Input layout: a single tensor of shape \code{(2 * d_spec)} where the first
#' \code{d_spec} columns are the SNV-normalised spectrum (uncovered bands
#' zero-filled) and the next \code{d_spec} columns are the binary coverage
#' mask (1 = covered, 0 = missing).
#'
#' Reconstruction loss is masked: only covered bands contribute. The
#' prediction head reads only from the latent z, so its output is invariant
#' to which bands were missing.
#'
#' @param d_spec Integer; number of bands per channel (e.g., 1076 for VisNIR)
#' @param latent_dim Integer; latent dim (default 16)
#' @param loss_weights Numeric vector length 2 (rec, pred); default c(0.3, 0.3)
#' @param l2_reg L2 regularisation (default 1e-5)
#' @param dropout_rate Prediction-head dropout (default 0.2)
#' @param encoder_dropout Encoder dropout (default 0.05)
#' @param learning_rate Initial Adam LR (default 1e-3)
#' @return Compiled keras Model
#' @export
build_maskedVAE <- function(d_spec,
                            latent_dim      = 32L,
                            loss_weights    = c(0.1, 1.0),
                            l2_reg          = 1e-5,
                            dropout_rate    = 0.2,
                            encoder_dropout = 0.05,
                            learning_rate   = 1e-3) {
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")
  .init_vae_python()

  d_in <- 2L * as.integer(d_spec)

  inp  <- keras3::layer_input(shape = d_in)
  enc  <- .build_encoder(inp, latent_dim, l2_reg, encoder_dropout)
  rec  <- .build_decoder(enc$z, as.integer(d_spec), l2_reg, name = "reconstruction")
  pred <- .build_pred_head(enc$z, l2_reg, dropout_rate)

  mdl <- keras3::keras_model(inputs = inp, outputs = list(rec, pred))
  mdl |> keras3::compile(
    optimizer    = keras3::optimizer_adam(learning_rate = learning_rate),
    loss         = list(.vae_env$masked_mse, "mse"),
    loss_weights = as.list(loss_weights)
  )
  mdl
}

# ---- Internal helpers ----------------------------------------------------

#' Pack [spectrum | mask] into a single 2*d-wide matrix for masked-VAE input
#' @keywords internal
.pack_masked_input <- function(spec, mask) {
  if (!identical(dim(spec), dim(mask)))
    stop(".pack_masked_input: dimensions differ.")
  cbind(spec, mask)
}

#' Pack [spectrum_target | mask] for the masked-MSE reconstruction target
#' @keywords internal
.pack_masked_target <- function(spec, mask) {
  cbind(spec, mask)
}

#' Extract the prediction head's output as a numeric vector
#' @keywords internal
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

# ---- Single-property trainer (standard) ----------------------------------

#' Train a (standard) soilVAE model for one soil property
#'
#' @param X Numeric matrix of preprocessed spectra (full coverage)
#' @param y Numeric vector of target values (length == nrow(X))
#' @param family_id Family ID (e.g. "OSSL_VisNIR_Cubist" - not used for VAE)
#' @param prop Soil property name (used for output filenames)
#' @param out_dir Root output directory (default "models")
#' @param epochs Max epochs (default 80)
#' @param batch_size Mini-batch size (default 32)
#' @param val_split Validation fraction (default 0.2)
#' @param patience_es Early-stopping patience (default 10)
#' @param patience_rl LR-reduction patience (default 5)
#' @param latent_dim Latent dim (default 16)
#' @param kl_beta Target KL weight after annealing (default 1e-3)
#' @param kl_warmup_epochs Annealing warmup (default 30)
#' @param min_n Minimum samples to attempt training (default 50)
#' @return Invisible TRUE on success, NULL if skipped
#' @export
train_soilVAE <- function(X, y, family_id, prop,
                          out_dir          = "models",
                          epochs           = 100L,
                          batch_size       = 64L,
                          val_split        = 0.15,
                          patience_es      = 20L,
                          patience_rl      = 7L,
                          latent_dim       = 32L,
                          kl_beta          = 5e-4,
                          kl_warmup_epochs = 30L,
                          min_n            = 50L) {

  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")
  n <- nrow(X)
  if (n < min_n) {
    message("Skip ", family_id, "/", prop, ": need >=", min_n, ", have ", n)
    return(invisible(NULL))
  }

  mean_y <- mean(y, na.rm = TRUE)
  sd_y   <- sd(y, na.rm = TRUE)
  if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z <- (y - mean_y) / sd_y

  sp    <- split_idx(n, seed = 42)
  X_tr  <- X[sp$train, , drop = FALSE]; y_tr_z <- y_z[sp$train]
  X_ca  <- X[sp$calib, , drop = FALSE]; y_ca   <- y[sp$calib]
  X_te  <- X[sp$test,  , drop = FALSE]; y_te   <- y[sp$test]

  mdl <- build_soilVAE(ncol(X_tr), latent_dim = latent_dim)
  cbs <- list(
    keras3::callback_early_stopping(
      monitor = "val_loss", patience = patience_es, restore_best_weights = TRUE),
    keras3::callback_reduce_lr_on_plateau(
      monitor = "val_loss", patience = patience_rl, factor = 0.5, min_lr = 1e-5),
    .vae_env$KLAnnealingCallback(
      target_beta   = kl_beta,
      warmup_epochs = as.integer(kl_warmup_epochs))
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

  model_dir <- file.path(out_dir, family_id, "models")
  dir_create(model_dir)
  base <- file.path(model_dir, prop)
  keras3::save_model(mdl, paste0(base, ".keras"), overwrite = TRUE)
  save_scaler(mean_y, sd_y, path_base = base)

  yhat_ca <- .extract_prediction(mdl, X_ca) * sd_y + mean_y
  res_ca  <- abs(y_ca - yhat_ca)
  q90 <- as.numeric(quantile(res_ca, 0.90, na.rm = TRUE))
  q95 <- as.numeric(quantile(res_ca, 0.95, na.rm = TRUE))

  yhat_te <- .extract_prediction(mdl, X_te) * sd_y + mean_y
  # Metrics in ORIGINAL scale: back-transform log1p so RMSE/R²/RPIQ/CCC are
  # comparable to Safanelli 2025 Tables 4/5 (reported in original units).
  if (prop %in% ossl_l1_log_transform) {
    mets <- metrics_from_y(expm1(y_te), expm1(yhat_te))
  } else {
    mets <- metrics_from_y(y_te, yhat_te)
  }

  encoder <- keras3::keras_model(inputs  = mdl$input,
                                 outputs = mdl$get_layer("z_mean")$output)
  Z_tr    <- predict(encoder, as.matrix(X_tr), verbose = 0)
  mu_z    <- colMeans(Z_tr, na.rm = TRUE)
  Sig_z   <- as.matrix(stats::cov(Z_tr, use = "pairwise.complete.obs") +
                         diag(1e-6, ncol(Z_tr)))
  thr95   <- stats::qchisq(0.95, df = length(mu_z))

  latent <- list(
    mu    = as.numeric(mu_z),
    Sigma = unname(asplit(Sig_z, 1L)),
    df    = length(mu_z),
    thr95 = thr95
  )
  fr <- feature_minmax(X)
  save_metrics_json(base,
                    metrics    = mets,
                    conf_q     = list(q90 = q90, q95 = q95),
                    feat_range = fr,
                    latent     = latent)

  message(sprintf("  [%s/%s] RMSE=%.3f Bias=%.3f R2=%.3f RPIQ=%.2f CCC=%.3f n=%d",
                  family_id, prop, mets$RMSE, mets$Bias,
                  mets$R2, mets$RPIQ, mets$CCC, n))
  invisible(TRUE)
}

# ---- Single-property trainer (masked) ------------------------------------

#' Train a masked soilVAE for one soil property
#'
#' Augments the training set via \code{spectral_dropout()} so the encoder
#' learns to handle partial-coverage inputs typical of real instruments
#' (NeoSpectra, Bruker MPA, PerkinElmer Spectrum One, etc.).
#'
#' @param X Numeric matrix (rows = samples, cols = bands; full-coverage,
#'   already SNV-normalised on the full grid)
#' @param y Numeric vector of target values
#' @param family_id Family ID (used for output directory)
#' @param prop Soil property name
#' @param wavegrid Numeric vector of wavelength positions (length == ncol(X))
#' @param sensor_type "visnir" or "mir" (controls augmentation cut-offs)
#' @param out_dir Output root (default "models")
#' @param epochs Max epochs (default 80)
#' @param batch_size Mini-batch (default 64)
#' @param val_split Validation fraction (default 0.15)
#' @param patience_es Early-stopping patience (default 10)
#' @param patience_rl LR-reduction patience (default 5)
#' @param latent_dim Latent dim (default 16)
#' @param kl_beta Target KL weight (default 1e-3)
#' @param kl_warmup_epochs Annealing epochs (default 30)
#' @param n_aug_copies Augmented copies per sample (default 2)
#' @param min_n Minimum samples (default 50)
#' @return Invisible TRUE on success, NULL on skip
#' @export
train_maskedVAE <- function(X, y, family_id, prop, wavegrid, sensor_type,
                            out_dir          = "models",
                            epochs           = 100L,
                            batch_size       = 64L,
                            val_split        = 0.15,
                            patience_es      = 20L,
                            patience_rl      = 7L,
                            latent_dim       = 32L,
                            kl_beta          = 5e-4,
                            kl_warmup_epochs = 30L,
                            n_aug_copies     = 3L,
                            min_n            = 50L) {
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")
  n <- nrow(X)
  if (n < min_n) {
    message("Skip ", family_id, "/", prop, ": need >=", min_n, ", have ", n)
    return(invisible(NULL))
  }

  mean_y <- mean(y, na.rm = TRUE)
  sd_y   <- sd(y, na.rm = TRUE)
  if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z <- (y - mean_y) / sd_y

  sp    <- split_idx(n, seed = 42)
  X_tr  <- X[sp$train, , drop = FALSE]; y_tr_z <- y_z[sp$train]
  X_ca  <- X[sp$calib, , drop = FALSE]; y_ca   <- y[sp$calib]
  X_te  <- X[sp$test,  , drop = FALSE]; y_te   <- y[sp$test]

  # Augment training set with synthetic partial-coverage copies
  aug   <- spectral_dropout(X_tr, wavegrid = wavegrid,
                            sensor_type = sensor_type,
                            n_copies = n_aug_copies, seed = 42L)
  y_aug <- rep(y_tr_z, times = n_aug_copies + 1L)

  X_in  <- .pack_masked_input(aug$X, aug$mask)
  Y_rec <- .pack_masked_target(aug$X, aug$mask)

  d_spec <- ncol(X_tr)
  mdl <- build_maskedVAE(d_spec, latent_dim = latent_dim)
  cbs <- list(
    keras3::callback_early_stopping(
      monitor = "val_loss", patience = patience_es, restore_best_weights = TRUE),
    keras3::callback_reduce_lr_on_plateau(
      monitor = "val_loss", patience = patience_rl, factor = 0.5, min_lr = 1e-5),
    .vae_env$KLAnnealingCallback(
      target_beta   = kl_beta,
      warmup_epochs = as.integer(kl_warmup_epochs))
  )
  mdl |> keras3::fit(
    x                = X_in,
    y                = list(Y_rec, y_aug),
    epochs           = epochs,
    batch_size       = batch_size,
    validation_split = val_split,
    callbacks        = cbs,
    verbose          = 0
  )

  model_dir <- file.path(out_dir, family_id, "models")
  dir_create(model_dir)
  base <- file.path(model_dir, prop)
  keras3::save_model(mdl, paste0(base, ".keras"), overwrite = TRUE)
  save_scaler(mean_y, sd_y, path_base = base)

  # Calibration / test inputs evaluated with full-coverage masks
  full_mask_ca <- matrix(1L, nrow(X_ca), d_spec)
  full_mask_te <- matrix(1L, nrow(X_te), d_spec)
  X_ca_in <- .pack_masked_input(X_ca, full_mask_ca)
  X_te_in <- .pack_masked_input(X_te, full_mask_te)

  yhat_ca <- .extract_prediction(mdl, X_ca_in) * sd_y + mean_y
  res_ca  <- abs(y_ca - yhat_ca)
  q90 <- as.numeric(quantile(res_ca, 0.90, na.rm = TRUE))
  q95 <- as.numeric(quantile(res_ca, 0.95, na.rm = TRUE))

  yhat_te <- .extract_prediction(mdl, X_te_in) * sd_y + mean_y
  # Metrics in ORIGINAL scale: back-transform log1p so RMSE/R²/RPIQ/CCC match
  # Safanelli 2025 reporting units (Tables 4/5).
  if (prop %in% ossl_l1_log_transform) {
    mets <- metrics_from_y(expm1(y_te), expm1(yhat_te))
  } else {
    mets <- metrics_from_y(y_te, yhat_te)
  }

  encoder <- keras3::keras_model(inputs  = mdl$input,
                                 outputs = mdl$get_layer("z_mean")$output)
  Z_tr    <- predict(encoder,
                     .pack_masked_input(X_tr, matrix(1L, nrow(X_tr), d_spec)),
                     verbose = 0)
  mu_z    <- colMeans(Z_tr, na.rm = TRUE)
  Sig_z   <- as.matrix(stats::cov(Z_tr, use = "pairwise.complete.obs") +
                         diag(1e-6, ncol(Z_tr)))
  thr95   <- stats::qchisq(0.95, df = length(mu_z))

  latent <- list(
    mu    = as.numeric(mu_z),
    Sigma = unname(asplit(Sig_z, 1L)),
    df    = length(mu_z),
    thr95 = thr95,
    masked = TRUE
  )
  fr <- feature_minmax(X)
  save_metrics_json(base,
                    metrics    = mets,
                    conf_q     = list(q90 = q90, q95 = q95),
                    feat_range = fr,
                    latent     = latent)

  message(sprintf("  [%s/%s] RMSE=%.3f Bias=%.3f R2=%.3f RPIQ=%.2f CCC=%.3f n=%d (+aug %d)",
                  family_id, prop, mets$RMSE, mets$Bias,
                  mets$R2, mets$RPIQ, mets$CCC, n, n_aug_copies))
  invisible(TRUE)
}

# ---- Multi-property OSSL trainer (dispatcher) ----------------------------

#' Train soilVAE models for all properties of one model family using OSSL data
#'
#' Dispatches to \code{train_soilVAE} or \code{train_maskedVAE} depending on
#' the family's \code{algorithm} field.
#'
#' @param family_id Family ID from \code{model_registry}
#' @param out_dir Output root (default "models")
#' @param cache_dir OSSL cache directory
#' @param properties Properties to train (default: all in family)
#' @param download_if_missing Auto-download OSSL data if missing
#' @param ... Extra args passed to the per-property trainer
#' @return Invisible named list of per-property results
#' @export
train_ossl_models <- function(family_id,
                              out_dir              = "models",
                              cache_dir            = ossl_cache_dir(),
                              properties           = NULL,
                              download_if_missing  = TRUE,
                              ...) {
  fam <- get_family(family_id)
  if (fam$source != "ossl")
    stop("family_id '", family_id, "' is not an OSSL family.")

  algo  <- fam$algorithm %||% "vae"
  props <- if (is.null(properties)) fam$properties else properties

  message("=== autoSpectra: Training [", fam$label, "] ===")
  message("  Algorithm   : ", algo)
  message("  Sensor type : ", fam$sensor_type)
  message("  Grid        : ", length(fam$wavegrid), " channels (",
          min(fam$wavegrid), " - ", max(fam$wavegrid), ")")
  message("  Pipeline    : ", paste(fam$preprocess, collapse = " -> "))
  message("  Properties  : ", length(props))

  joined_df <- ossl_prepare(
    sensor_type         = fam$sensor_type,
    cache_dir           = cache_dir,
    properties          = props,
    download_if_missing = download_if_missing
  )

  wl_info <- get_wavelengths(joined_df, id_col = "Soil_ID")
  X_src   <- as.matrix(joined_df[, wl_info$cols, drop = FALSE])
  X_res   <- resample_to_grid(X_src, src_wl = wl_info$wl,
                              target_wl = fam$wavegrid)

  message("  Preprocessing ...")
  full_mask <- matrix(1L, nrow(X_res), ncol(X_res))
  X_proc <- apply_pipeline(X_res, fam$preprocess, mask = full_mask)

  results <- list()
  for (prop in props) {
    y <- suppressWarnings(as.numeric(joined_df[[prop]]))

    # Apply Safanelli's log1p transform where appropriate
    if (prop %in% ossl_l1_log_transform) {
      ok <- is.finite(y) & y > -1
      y[!ok] <- NA_real_
      y <- log1p(y)
    }

    keep <- is.finite(y) & is.finite(rowSums(X_proc))
    if (sum(keep) < 50L) {
      message("  Skip [", prop, "]: only ", sum(keep), " complete rows")
      results[[prop]] <- NULL
      next
    }
    message("  Training [", prop, "] (n=", sum(keep), ") ...")
    res <- try({
      if (identical(algo, "vae_masked")) {
        train_maskedVAE(X_proc[keep, ], y[keep], family_id, prop,
                        wavegrid    = fam$wavegrid,
                        sensor_type = fam$sensor_type,
                        out_dir     = out_dir, ...)
      } else {
        train_soilVAE(X_proc[keep, ], y[keep], family_id, prop,
                      out_dir = out_dir, ...)
      }
    }, silent = FALSE)
    results[[prop]] <- res
  }
  message("=== Done: ", family_id, " ===")
  invisible(results)
}

# ---- Model loader --------------------------------------------------------

#' Load a saved soilVAE model and its scaler
#'
#' @param family_id Family identifier
#' @param prop Soil property name
#' @param model_dir Root model directory (default "models")
#' @return List with \code{model} (keras) and \code{scaler}
#' @export
load_soilVAE <- function(family_id, prop, model_dir = "models") {
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")
  .init_vae_python()

  base <- file.path(model_dir, family_id, "models", prop)
  # Prefer .keras (keras3 native); fall back to .h5 for older artefacts
  path <- if (file.exists(paste0(base, ".keras"))) paste0(base, ".keras")
          else if (file.exists(paste0(base, ".h5"))) paste0(base, ".h5")
          else stop("Model file not found for ", family_id, "/", prop)

  custom <- list("Sampling" = .vae_env$Sampling,
                 "masked_mse" = .vae_env$masked_mse)
  mdl <- keras3::load_model(path, custom_objects = custom)

  scaler <- if (file.exists(paste0(base, "_scaler.rds"))) {
    readRDS(paste0(base, "_scaler.rds"))
  } else {
    list(mean = 0, sd = 1, method = "none")
  }
  list(model = mdl, scaler = scaler)
}

# Note: `%||%` operator is defined in R/utils.R (loaded before this file).
