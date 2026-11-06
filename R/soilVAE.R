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

    Outside training the layer returns z_mean. Sampling at inference made
    every prediction random: two identical predict_soil() calls differed by a
    median 1-2%, up to 0.75 %OC at the 95th percentile.
    '''
    def call(self, inputs, training=None):
        z_mean, z_log_var = inputs
        z_log_var = ops.clip(z_log_var, -10.0, 6.0)
        kl  = -0.5 * ops.sum(
            1.0 + z_log_var - ops.square(z_mean) - ops.exp(z_log_var),
            axis=-1
        )
        beta = float(getattr(builtins, '_VAE_BETA', 0.0))
        self.add_loss(beta * ops.mean(kl))
        if not training:
            return z_mean
        eps = k_random.normal(shape=ops.shape(z_mean))
        return z_mean + ops.exp(0.5 * z_log_var) * eps

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

# ---- Training callbacks (shared by every VAE trainer) ---------------------

# Early stopping and LR reduction watch the PREDICTION loss on held-out rows,
# and early stopping waits for the KL warm-up to finish. The total val_loss is
# the wrong signal: its KL term's weight grows every epoch of the warm-up, so
# the total rises while the model improves, and early stopping used to restore
# weights from one of the first epochs (VisNIR oc stopped after ~8 minutes).
.vae_callbacks <- function(patience_es, patience_rl, kl_beta, kl_warmup_epochs) {
  list(
    keras3::callback_early_stopping(
      monitor = "val_prediction_loss", mode = "min", patience = patience_es,
      restore_best_weights = TRUE,
      start_from_epoch = as.integer(kl_warmup_epochs)),
    keras3::callback_reduce_lr_on_plateau(
      monitor = "val_prediction_loss", mode = "min", patience = patience_rl,
      factor = 0.5, min_lr = 1e-5),
    .vae_env$KLAnnealingCallback(
      target_beta   = kl_beta,
      warmup_epochs = as.integer(kl_warmup_epochs))
  )
}

# ---- Internal helpers ----------------------------------------------------

#' Pack the spectrum and its coverage mask side by side into a single 2*d-wide matrix for masked-VAE input
#' @keywords internal
.pack_masked_input <- function(spec, mask) {
  if (!identical(dim(spec), dim(mask)))
    stop(".pack_masked_input: dimensions differ.")
  cbind(spec, mask)
}

#' Pack the target spectrum and its mask side by side for the masked-MSE reconstruction target
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
#' @param groups Optional group label per row; rows sharing one stay in the
#'   same train/calibration/test set
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
                          min_n            = 50L,
                          groups           = NULL) {

  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")
  n <- nrow(X)
  if (n < min_n) {
    message("Skip ", family_id, "/", prop, ": need >=", min_n, ", have ", n)
    return(invisible(NULL))
  }

  sp    <- split_idx(n, seed = 42, groups = groups)
  # Target scaling from training rows only
  mean_y <- mean(y[sp$train], na.rm = TRUE)
  sd_y   <- stats::sd(y[sp$train], na.rm = TRUE)
  if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z <- (y - mean_y) / sd_y

  X_tr  <- X[sp$train, , drop = FALSE]; y_tr_z <- y_z[sp$train]
  X_ca  <- X[sp$calib, , drop = FALSE]; y_ca   <- y[sp$calib]
  X_te  <- X[sp$test,  , drop = FALSE]; y_te   <- y[sp$test]

  mdl <- build_soilVAE(ncol(X_tr), latent_dim = latent_dim)
  cbs <- .vae_callbacks(patience_es, patience_rl, kl_beta, kl_warmup_epochs)
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
  save_scaler(mean_y, sd_y, path_base = base, bounds = target_bounds(y[sp$train]))

  yhat_ca <- .extract_prediction(mdl, X_ca) * sd_y + mean_y
  res_ca  <- abs(y_ca - yhat_ca)
  q90 <- as.numeric(quantile(res_ca, 0.90, na.rm = TRUE))
  q95 <- as.numeric(quantile(res_ca, 0.95, na.rm = TRUE))

  yhat_te <- .extract_prediction(mdl, X_te) * sd_y + mean_y
  # Metrics in ORIGINAL scale: back-transform log1p so RMSE/R²/RPIQ/CCC are
  # comparable to Safanelli 2025 Tables 4/5 (reported in original units).
  mets <- metrics_original_scale(y_te, yhat_te, prop, bounds = target_bounds(y[sp$train]))

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
  save_metrics_json(base,
                    metrics    = mets,
                    conf_q     = list(q90 = q90, q95 = q95),
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
#' \code{X} holds the resampled spectra BEFORE preprocessing, with uncovered
#' bands as \code{NA}. The trainer preprocesses every copy itself with
#' \code{prepare_model_input()}, the call \code{predict_soil()} makes, so a
#' simulated partial-coverage copy is exactly what the model would receive
#' from a real instrument with that range.
#'
#' Besides the usual test metrics (true coverage), the metrics file records
#' \code{metrics_by_coverage}: the same test rows scored as if measured by each
#' simulated instrument range, so partial-coverage performance is measured
#' rather than assumed.
#'
#' @param X Numeric matrix of resampled spectra (rows = samples, cols = bands),
#'   not yet preprocessed; uncovered bands \code{NA}
#' @param y Numeric target vector (modelling scale)
#' @param family_id Family ID (used for output directory)
#' @param prop Soil property name
#' @param wavegrid Numeric vector of wavelength positions (length == ncol(X))
#' @param sensor_type "visnir" or "mir" (controls augmentation ranges)
#' @param out_dir Output root (default "models")
#' @param epochs Max epochs
#' @param batch_size Mini-batch size
#' @param val_split Fraction of training samples held out for early stopping.
#'   They are separate samples, never augmented copies of training samples.
#' @param patience_es Early-stopping patience
#' @param patience_rl LR-reduction patience
#' @param latent_dim Latent dim
#' @param kl_beta Target KL weight
#' @param kl_warmup_epochs Annealing epochs
#' @param n_aug_copies Augmented copies per sample
#' @param min_n Minimum samples (default 50)
#' @param mask Optional coverage mask (same shape as X); derived from the
#'   \code{NA} pattern of \code{X} if omitted
#' @param groups Optional group label per row; rows sharing one stay in the same
#'   train/calibration/test/validation set
#' @param preprocess Preprocessing steps; default the family's
#'   \code{preprocess} field, or per-mask SNV alone for an unknown family
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
                            min_n            = 50L,
                            mask             = NULL,
                            groups           = NULL,
                            preprocess       = NULL) {
  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")
  if (is.null(preprocess))
    preprocess <- if (!is.null(model_registry[[family_id]]))
      model_registry[[family_id]]$preprocess else "SNV_MASKED"
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < min_n) {
    message("Skip ", family_id, "/", prop, ": need >=", min_n, ", have ", n)
    return(invisible(NULL))
  }
  if (is.null(mask))   mask   <- mask_from_coverage(X)
  if (is.null(groups)) groups <- seq_len(n)

  sp <- split_idx(n, seed = 42, groups = groups)
  tr <- sp$train
  # Early-stopping rows are whole samples taken out of training. Keras'
  # validation_split would take the last rows of the augmented matrix, i.e.
  # copies of samples whose other copies are being fitted.
  vf     <- .grouped_partition(length(tr), c(1 - val_split, val_split),
                               groups[tr], seed = 43L)
  fit_ix <- tr[vf[[1]]]
  val_ix <- tr[vf[[2]]]

  # Target scaling from training rows only
  mean_y <- mean(y[tr], na.rm = TRUE)
  sd_y   <- stats::sd(y[tr], na.rm = TRUE)
  if (!is.finite(sd_y) || sd_y == 0) sd_y <- 1
  y_z    <- (y - mean_y) / sd_y
  bounds <- target_bounds(y[tr])

  .augmented <- function(rows, seed) {
    aug <- spectral_dropout(X[rows, , drop = FALSE], wavegrid = wavegrid,
                            sensor_type = sensor_type, n_copies = n_aug_copies,
                            seed = seed, mask = mask[rows, , drop = FALSE],
                            steps = preprocess)
    list(x = .pack_masked_input(aug$X, aug$mask),
         y = list(.pack_masked_target(aug$X, aug$mask),
                  rep(y_z[rows], times = n_aug_copies + 1L)))
  }
  fit_set <- .augmented(fit_ix, 42L)
  val_set <- .augmented(val_ix, 43L)

  d_spec <- ncol(X)
  mdl <- build_maskedVAE(d_spec, latent_dim = latent_dim)
  cbs <- .vae_callbacks(patience_es, patience_rl, kl_beta, kl_warmup_epochs)
  hist <- keras3::fit(
    mdl,
    x               = fit_set$x,
    y               = fit_set$y,
    validation_data = list(val_set$x, val_set$y),
    epochs          = epochs,
    batch_size      = batch_size,
    callbacks       = cbs,
    verbose         = 0
  )
  rm(fit_set, val_set); invisible(gc(verbose = FALSE))
  val_pred   <- as.numeric(unlist(hist$metrics$val_prediction_loss))
  epochs_run <- length(val_pred)
  best_epoch <- if (epochs_run) which.min(val_pred) else NA_integer_

  model_dir <- file.path(out_dir, family_id, "models")
  dir_create(model_dir)
  base <- file.path(model_dir, prop)
  keras3::save_model(mdl, paste0(base, ".keras"), overwrite = TRUE)
  save_scaler(mean_y, sd_y, path_base = base, bounds = bounds)

  # Model input for given rows under a given mask: the predict_soil() path
  .input <- function(rows, mk) {
    .pack_masked_input(.model_input_under_mask(X[rows, , drop = FALSE], mk,
                                               preprocess), mk)
  }
  .predict_rows <- function(rows, mk) {
    .extract_prediction(mdl, .input(rows, mk)) * sd_y + mean_y
  }

  ca <- sp$calib; te <- sp$test
  res_ca <- abs(y[ca] - .predict_rows(ca, mask[ca, , drop = FALSE]))
  q90 <- as.numeric(stats::quantile(res_ca, 0.90, na.rm = TRUE))
  q95 <- as.numeric(stats::quantile(res_ca, 0.95, na.rm = TRUE))

  # Test metrics with each row's real coverage, in original units
  mets <- metrics_original_scale(y[te], .predict_rows(te, mask[te, , drop = FALSE]),
                                 prop, bounds = bounds)

  # The same test rows as each simulated instrument would have measured them
  by_cov <- list()
  for (nm in names(cuts <- .coverage_cuts(sensor_type, wavegrid))) {
    in_cut <- as.integer(wavegrid >= cuts[[nm]][1] & wavegrid <= cuts[[nm]][2])
    mk     <- sweep(mask[te, , drop = FALSE], 2, in_cut, "*")
    ok     <- rowSums(mk) >= 3L
    if (sum(ok) < 10L) next
    m <- metrics_original_scale(y[te][ok], .predict_rows(te[ok], mk[ok, , drop = FALSE]),
                                prop, bounds = bounds)
    m$n_test <- sum(ok)
    by_cov[[nm]] <- m
  }

  # Latent-space applicability domain from training rows at their real
  # coverage. The chi-squared limit assumes a standard Gaussian latent space,
  # which the trained encoder does not produce (it flagged ~13% of OSSL rather
  # than 5%), so the limit used is the empirical 95th percentile.
  encoder <- keras3::keras_model(inputs  = mdl$input,
                                 outputs = mdl$get_layer("z_mean")$output)
  Z_tr    <- predict(encoder, .input(tr, mask[tr, , drop = FALSE]), verbose = 0)
  mu_z    <- colMeans(Z_tr, na.rm = TRUE)
  Sig_z   <- as.matrix(stats::cov(Z_tr, use = "pairwise.complete.obs") +
                         diag(1e-6, ncol(Z_tr)))
  d2_tr   <- mahalanobis_d2(Z_tr, mu_z, Sig_z)

  latent <- list(
    mu          = as.numeric(mu_z),
    Sigma       = unname(asplit(Sig_z, 1L)),
    df          = length(mu_z),
    thr95       = as.numeric(stats::quantile(d2_tr, 0.95, names = FALSE)),
    thr95_method = "empirical_q95_train",
    thr95_chi2  = stats::qchisq(0.95, df = length(mu_z)),
    masked      = TRUE,
    n_train     = length(fit_ix),
    n_val       = length(val_ix),
    n_test      = length(te),
    epochs_run  = epochs_run,
    best_epoch  = best_epoch,
    metrics_by_coverage = by_cov
  )
  save_metrics_json(base,
                    metrics    = mets,
                    conf_q     = list(q90 = q90, q95 = q95),
                    latent     = latent)

  message(sprintf("  [%s/%s] RMSE=%.3f Bias=%.3f R2=%.3f RPIQ=%.2f CCC=%.3f n=%d (+aug %d) epochs=%d best=%d",
                  family_id, prop, mets$RMSE, mets$Bias,
                  mets$R2, mets$RPIQ, mets$CCC, n, n_aug_copies,
                  epochs_run, best_epoch))
  for (nm in names(by_cov))
    message(sprintf("      coverage %-14s R2=%.3f n=%d", nm, by_cov[[nm]]$R2,
                    by_cov[[nm]]$n_test))
  invisible(TRUE)
}

# ---- Multi-property OSSL trainer (dispatcher) ----------------------------

#' Train models for all properties of one model family using OSSL data
#'
#' Dispatches to \code{train_maskedVAE()}, \code{train_cubist_pca()} or
#' \code{train_soilVAE()} depending on the family's \code{algorithm} field,
#' with the same preprocessing \code{predict_soil()} applies.
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

  message("  Preprocessing ...")
  inp <- .family_training_inputs(fam, joined_df)

  results <- list()
  for (prop in props) {
    y    <- apply_target_transform(joined_df[[prop]], prop)
    keep <- is.finite(y) & inp$has_spectra
    if (sum(keep) < 50L) {
      message("  Skip [", prop, "]: only ", sum(keep), " complete rows")
      results[[prop]] <- NULL
      next
    }
    message("  Training [", prop, "] (n=", sum(keep), ") ...")
    res <- try({
      if (identical(algo, "vae_masked")) {
        train_maskedVAE(inp$X_raw[keep, , drop = FALSE], y[keep], family_id, prop,
                        wavegrid    = fam$wavegrid,
                        sensor_type = fam$sensor_type,
                        out_dir     = out_dir,
                        mask        = inp$mask[keep, , drop = FALSE],
                        groups      = inp$groups[keep], ...)
      } else if (identical(algo, "cubist_pca")) {
        train_cubist_pca(inp$X[keep, , drop = FALSE], y[keep], family_id, prop,
                         out_dir = out_dir, groups = inp$groups[keep], ...)
      } else if (identical(algo, "pca_vae")) {
        train_pca_vae(inp$X[keep, , drop = FALSE], y[keep], family_id, prop,
                      out_dir = out_dir, groups = inp$groups[keep], ...)
      } else {
        train_soilVAE(inp$X[keep, , drop = FALSE], y[keep], family_id, prop,
                      out_dir = out_dir, groups = inp$groups[keep], ...)
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
