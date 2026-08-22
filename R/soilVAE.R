# R/soilVAE.R — soilVAE Variational Autoencoder (architecture + training)
#
# Two model builders:
#
#   build_soilVAE(d_in)      — Standard VAE for MIR (consistent spectral coverage)
#     Input (d_in):
#       Encoder: Dense(256)+BN → Dense(128)+BN → Dense(64)+BN+Dropout → Dense(32)+BN
#       Latent:  z_mean, z_log_var → Sampling() (reparameterisation trick)
#       Reconstruction head: Dense(32)+BN → Dense(64)+BN → Dense(128)+BN → Dense(d_in)
#       Prediction head:     Dense(128)+Dropout → Dense(64) → Dense(32) → Dense(1)
#
#   build_maskedVAE(d_spec)  — Sensor-agnostic VAE for VisNIR (variable coverage)
#     Input (2 × d_spec): [SNV-spectrum | binary-mask] concatenation
#       Encoder: Dense(512)+BN → Dense(256)+BN → Dense(128)+BN+Dropout → Dense(64)+BN
#       Latent:  same reparameterisation trick
#       Reconstruction head: outputs d_spec (spectrum only, not mask)
#       Prediction head:     same as above
#
# KL divergence: injected as add_loss() inside the Sampling Python layer.
# KL annealing:  KLAnnealingCallback ramps beta 0 → kl_beta over warmup_epochs.
#                Prevents posterior collapse before reconstruction is learned.
# BatchNorm:     on all encoder/decoder Dense layers for OSSL-scale stability.
# L2 reg:        on all Dense layers to prevent overfitting.
# z_mean:        (deterministic) used for Mahalanobis applicability-domain checks.
#
# Loss = MSE(reconstruction) × w1  +  MSE(prediction) × w2  +  beta × KL
#
# kl_beta = 1e-3: meaningful KL regularisation — forces the encoder towards N(0,I)
# while preserving reconstruction and prediction quality.  Too small (e.g. 4e-5)
# renders the VAE effectively deterministic; too large collapses the latent space.
#
# A scaler (z-score on y) is saved alongside each model.
# Conformal calibration quantiles (q90, q95) are saved for uncertainty intervals.
# Latent space statistics (μ, Σ) enable Mahalanobis applicability-domain checks.

# ---- Python VAE classes (defined once per R session) ------------------------

.vae_env <- new.env(parent = emptyenv())
.vae_env$initialized <- FALSE

.python_support_code <- "import importlib.util\nimport sys\nmods = ('keras', 'tensorflow', 'numpy')\nmissing = [m for m in mods if importlib.util.find_spec(m) is None]\nif missing:\n    sys.stderr.write('missing:' + ','.join(missing))\n    raise SystemExit(1)\nprint('ok')"

.python_has_keras <- function(python) {
  if (is.null(python) || !nzchar(python) || !file.exists(python))
    return(FALSE)

  out <- tryCatch(
    suppressWarnings(system2(
      python,
      c("-c", shQuote(.python_support_code)),
      stdout = TRUE,
      stderr = TRUE
    ),
    ),
    error = function(e) structure(character(0), status = 1L)
  )

  status <- attr(out, "status")
  if (is.null(status) || length(status) == 0L) status <- 0L
  identical(status, 0L)
}

.candidate_python_paths <- function() {
  paths <- character(0)

  env_python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (nzchar(env_python)) paths <- c(paths, env_python)

  if (requireNamespace("reticulate", quietly = TRUE)) {
    venvs <- tryCatch(reticulate::virtualenv_list(), error = function(e) character(0))
    preferred <- c("r-tf311", "r-tensorflow", "r-tf", "r-reticulate")
    ordered <- c(intersect(preferred, venvs), setdiff(venvs, preferred))
    for (venv in ordered) {
      py <- tryCatch(reticulate::virtualenv_python(venv), error = function(e) "")
      if (nzchar(py)) paths <- c(paths, py)
    }
  }

  py_bin <- Sys.which("python")
  if (nzchar(py_bin)) paths <- c(paths, py_bin)

  unique(paths[nzchar(paths)])
}

.autospectra_seed_value <- function(seed) {
  seed <- suppressWarnings(as.integer(seed)[1L])
  if (!is.finite(seed)) return(NULL)
  seed
}

.autospectra_set_random_seed <- function(seed) {
  seed <- .autospectra_seed_value(seed)
  if (is.null(seed))
    return(invisible(NULL))

  set.seed(seed)

  if (!nzchar(Sys.getenv("PYTHONHASHSEED", unset = "")))
    Sys.setenv(PYTHONHASHSEED = as.character(seed))
  if (!nzchar(Sys.getenv("TF_DETERMINISTIC_OPS", unset = "")))
    Sys.setenv(TF_DETERMINISTIC_OPS = "1")
  if (!nzchar(Sys.getenv("TF_CUDNN_DETERMINISTIC", unset = "")))
    Sys.setenv(TF_CUDNN_DETERMINISTIC = "1")

  if (requireNamespace("keras3", quietly = TRUE)) {
    try(keras3::set_random_seed(seed), silent = TRUE)
  }

  if (requireNamespace("reticulate", quietly = TRUE) &&
      reticulate::py_available(initialize = FALSE)) {
    try(reticulate::py_run_string(sprintf("
import os
import random
import numpy as np
import tensorflow as tf
seed = %d
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)
tf.random.set_seed(seed)
try:
    tf.config.experimental.enable_op_determinism()
except Exception:
    pass
", seed)), silent = TRUE)
  }

  .vae_env$seed <- seed
  invisible(seed)
}

.autospectra_configure_python <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE))
    return(invisible(NULL))

  if (!nzchar(Sys.getenv("TF_CPP_MIN_LOG_LEVEL", unset = "")))
    Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
  if (!nzchar(Sys.getenv("TF_DETERMINISTIC_OPS", unset = "")))
    Sys.setenv(TF_DETERMINISTIC_OPS = "1")
  if (!nzchar(Sys.getenv("TF_CUDNN_DETERMINISTIC", unset = "")))
    Sys.setenv(TF_CUDNN_DETERMINISTIC = "1")
  if (!nzchar(Sys.getenv("PYTHONHASHSEED", unset = "")) &&
      !is.null(.vae_env$seed))
    Sys.setenv(PYTHONHASHSEED = as.character(.vae_env$seed))

  if (reticulate::py_available(initialize = FALSE)) {
    .vae_env$python_path <- tryCatch(reticulate::py_config()$python, error = function(e) NA_character_)
    ok <- tryCatch(
      reticulate::py_module_available("keras") &&
        reticulate::py_module_available("tensorflow") &&
        reticulate::py_module_available("numpy"),
      error = function(e) FALSE
    )
    if (ok)
      return(invisible(.vae_env$python_path))

    stop(
      "reticulate is already initialized with Python interpreter '",
      .vae_env$python_path,
      "' but it does not provide 'keras', 'tensorflow', and 'numpy'. ",
      "Restart the R session or set RETICULATE_PYTHON before loading the model code."
    )
  }

  candidates <- .candidate_python_paths()
  for (python in candidates) {
    if (!.python_has_keras(python)) next
    reticulate::use_python(python, required = TRUE)
    .vae_env$python_path <- python
    return(invisible(python))
  }

  stop(
    "Could not find a Python environment with 'keras', 'tensorflow', and 'numpy'. ",
    "Set RETICULATE_PYTHON to a compatible interpreter or install them in one of the detected virtualenvs."
  )
}

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

  .autospectra_configure_python()

  reticulate::py_run_string("
import os
import contextlib
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '3')
import logging
try:
    from absl import logging as absl_logging
    absl_logging.set_verbosity(absl_logging.ERROR)
except Exception:
    absl_logging = None
import tensorflow as tf
import numpy as np
try:
    tf.config.experimental.enable_op_determinism()
except Exception:
    pass
tf.get_logger().setLevel('ERROR')
logging.getLogger('tensorflow').setLevel(logging.ERROR)
logging.getLogger('absl').setLevel(logging.ERROR)
import builtins
import keras
from keras import layers, ops, random as k_random

@contextlib.contextmanager
def silence_stderr_fd():
    saved_stderr = os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved_stderr, 2)
        os.close(saved_stderr)
        os.close(devnull)

def _to_numpy(obj, preserve_sequence=False):
    if obj is None:
        return None
    if isinstance(obj, tuple):
        if preserve_sequence:
            return tuple(_to_numpy(item, preserve_sequence=False) for item in obj)
        return tuple(_to_numpy(item, preserve_sequence=False) for item in obj)
    if isinstance(obj, list):
        if preserve_sequence:
            return [_to_numpy(item, preserve_sequence=False) for item in obj]
        try:
            arr = np.asarray(obj)
            if arr.dtype != object:
                return arr
        except Exception:
            pass
        return [_to_numpy(item, preserve_sequence=False) for item in obj]
    return np.asarray(obj)

def keras_fit_quiet(model, x, y=None, y_recon=None, y_pred=None, **kwargs):
    x = _to_numpy(x, preserve_sequence=isinstance(x, (list, tuple)))
    if y is None and y_recon is not None and y_pred is not None:
        y = [_to_numpy(y_recon), _to_numpy(y_pred)]
    elif isinstance(y, (list, tuple)):
        y = _to_numpy(y, preserve_sequence=True)
    elif y is not None:
        y = _to_numpy(y)
    if 'sample_weight' in kwargs and kwargs['sample_weight'] is not None:
        kwargs['sample_weight'] = _to_numpy(
            kwargs['sample_weight'],
            preserve_sequence=isinstance(kwargs['sample_weight'], (list, tuple))
        )
    with silence_stderr_fd():
        return model.fit(x=x, y=y, **kwargs)

def keras_predict_quiet(model, *args, **kwargs):
    if args:
        args = (_to_numpy(args[0], preserve_sequence=isinstance(args[0], (list, tuple))),) + tuple(args[1:])
    with silence_stderr_fd():
        return model.predict(*args, **kwargs)

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
  .vae_env$keras_fit_quiet     <- py_main$keras_fit_quiet
  .vae_env$keras_predict_quiet <- py_main$keras_predict_quiet
  .vae_env$initialized         <- TRUE
  .autospectra_set_random_seed(.vae_env$seed)
  invisible(NULL)
}

.keras_fit_quiet <- function(model, x, y = NULL, y_recon = NULL, y_pred = NULL, ...) {
  .init_vae_python()
  dots <- list(...)
  int_args <- c(
    "epochs", "batch_size", "verbose", "initial_epoch",
    "steps_per_epoch", "validation_steps"
  )
  for (nm in intersect(names(dots), int_args)) {
    dots[[nm]] <- as.integer(dots[[nm]])
  }
  if ("validation_split" %in% names(dots)) {
    dots[["validation_split"]] <- as.numeric(dots[["validation_split"]])
  }
  to_keras_array <- function(obj) {
    if (is.list(obj)) {
      return(lapply(obj, to_keras_array))
    }
    if (is.matrix(obj)) return(obj)
    if (is.data.frame(obj)) return(as.matrix(obj))
    if (is.atomic(obj)) return(obj)
    obj
  }
  if (is.null(y) && !is.null(y_recon) && !is.null(y_pred)) {
    y <- list(y_recon, y_pred)
  }
  if (is.list(y) && length(y) == 2L && !is.list(y[[1L]])) {
    args <- c(list(
      model = model,
      x = to_keras_array(x),
      y_recon = to_keras_array(y[[1L]]),
      y_pred = to_keras_array(y[[2L]])
    ), dots)
    return(do.call(.vae_env$keras_fit_quiet, args))
  }
  args <- c(list(
    model = model,
    x = to_keras_array(x),
    y = to_keras_array(y)
  ), dots)
  do.call(.vae_env$keras_fit_quiet, args)
}

.keras_predict_quiet <- function(model, X, verbose = 0L) {
  .init_vae_python()
  if (is.list(X)) {
    X <- lapply(X, function(x) if (is.matrix(x)) x else as.matrix(x))
  } else {
    X <- as.matrix(X)
  }
  .vae_env$keras_predict_quiet(
    model,
    X,
    verbose = as.integer(verbose)
  )
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

# ---- Masked VAE builder -----------------------------------------------------

#' Build the Masked soilVAE for sensor-agnostic training
#'
#' The encoder receives the concatenation of the SNV-normalised spectrum and a
#' binary validity mask (\code{2 × d_spec} inputs), so the network can
#' distinguish genuinely zero absorbance from absent spectral bands.
#' The reconstruction head outputs \code{d_spec} values (spectrum only).
#' Training targets: reconstruction = masked spectrum (zeros for absent bands);
#' prediction = z-scored soil property.
#'
#' Training call pattern:
#' \preformatted{
#'   fit(x = cbind(X_snv, mask),       # (n, 2*d_spec)
#'       y = list(X_snv, y_z), ...)    # recon + pred targets
#' }
#'
#' @param d_spec        Integer; number of spectral bands in the canonical grid.
#' @param latent_dim    Integer; latent space dimension (default 16).
#' @param loss_weights  Numeric vector length 2; reconstruction and prediction
#'   weights (default c(0.3, 0.3)).
#' @param l2_reg        L2 regularisation coefficient (default 1e-5).
#' @param dropout_rate  Dropout rate in the prediction head (default 0.2).
#' @param encoder_dropout Dropout inside the encoder (default 0.05).
#' @param learning_rate Initial Adam learning rate (default 1e-3).
#' @return A compiled keras Model with input shape \code{(2 × d_spec)}.
#' @export
build_maskedVAE <- function(d_spec,
                             latent_dim      = 16L,
                             loss_weights    = c(0.3, 0.3),
                             l2_reg          = 1e-5,
                             dropout_rate    = 0.2,
                             encoder_dropout = 0.05,
                             learning_rate   = 1e-3) {

  if (!requireNamespace("keras3", quietly = TRUE))
    stop("Package 'keras3' is required.")

  .init_vae_python()

  d_spec <- as.integer(d_spec)
  d_in   <- 2L * d_spec          # [SNV-spectrum | binary-mask] concatenated
  reg    <- keras3::regularizer_l2(l2 = l2_reg)

  inp <- keras3::layer_input(shape = d_in, name = "spec_mask")

  # ---- Encoder: wider first layer to handle the doubled input ---------------
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

  # ---- Latent space (reparametrisation trick) --------------------------------
  z_mean    <- keras3::layer_dense(enc, latent_dim, name = "z_mean")
  z_log_var <- keras3::layer_dense(enc, latent_dim, name = "z_log_var")
  z         <- .vae_env$Sampling()(list(z_mean, z_log_var))

  # ---- Reconstruction head: outputs d_spec (spectrum only, not mask) ---------
  rec <- z |>
    keras3::layer_dense(64,     activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(128,    activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(256,    activation = "relu", kernel_regularizer = reg) |>
    keras3::layer_batch_normalization() |>
    keras3::layer_dense(d_spec, activation = "linear", name = "reconstruction")

  # ---- Prediction head -------------------------------------------------------
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
                          kl_beta          = 1e-3,
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
  .keras_fit_quiet(
    mdl,
    x                = as.matrix(X_tr),
    y_recon          = as.matrix(X_tr),
    y_pred           = y_tr_z,
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
  keras3::save_model(mdl, .default_model_path(base), overwrite = TRUE)
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
  Z_tr    <- .keras_predict_quiet(encoder, X_tr, verbose = 0)
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
  raw <- .keras_predict_quiet(mdl, X, verbose = 0)
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
  model_path <- .resolve_model_path(base)

  mdl <- keras3::load_model(
    model_path,
    compile = FALSE,
    custom_objects = list("Sampling" = .vae_env$Sampling)
  )

  scaler <- if (file.exists(paste0(base, "_scaler.rds"))) {
    readRDS(paste0(base, "_scaler.rds"))
  } else {
    list(mean = 0, sd = 1, method = "none", target_transform = "none")
  }
  list(model = mdl, scaler = scaler)
}
