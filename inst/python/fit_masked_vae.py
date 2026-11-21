"""Fit a masked soilVAE in its own Python process (used by train_maskedVAE()).

R cannot load PyTorch in the same process: data.table and torch each bring an
OpenMP runtime and the second one aborts or segfaults the session. R writes the
training matrices to a private temporary folder, this script fits the model
(on the Apple GPU via the Keras torch backend when available) and saves it as a
backend-independent .keras file that R then loads with its TensorFlow backend.

Usage: python fit_masked_vae.py <job_dir>
       python fit_masked_vae.py --describe <d_spec> <latent_dim>
         prints the layer list as JSON (used by the tests)

<job_dir>/job.json      settings written by R
<job_dir>/fit_x.f32     fit input (rows x 2d, float32, column-major); also the
                        reconstruction target, which is identical
<job_dir>/fit_y.f32     standardised target (rows)
<job_dir>/val_x.f32, val_y.f32
Writes:
<job_dir>/heartbeat     touched every few batches (R's watchdog reads its time)
<job_dir>/progress.log  one line per epoch
<job_dir>/model.weights.h5, check.f32 (predictions for the first 256 val rows),
                        history.json, done
The architecture must match .build_encoder/.build_decoder/.build_pred_head in
R/soilVAE.R; tests/testthat/test-fit-process.R compares the two.
"""
import builtins
import json
import os
import sys
import time

job_dir = None if sys.argv[1] == "--describe" else sys.argv[1]
job = {}
if job_dir is not None:
    with open(os.path.join(job_dir, "job.json")) as fh:
        job = json.load(fh)

os.environ.setdefault("KERAS_BACKEND", job.get("backend", "torch"))
import numpy as np
import keras
from keras import layers, ops, random as k_random


def _read(name, ncol):
    a = np.fromfile(os.path.join(job_dir, name), dtype="<f4")
    if ncol == 1:
        return a.reshape(-1, 1)
    return a.reshape(ncol, -1).T  # R writes column-major


builtins._VAE_BETA_VAR = keras.Variable(0.0, trainable=False, dtype="float32",
                                        name="vae_kl_beta")


def _set_vae_beta(value):
    builtins._VAE_BETA = float(value)
    builtins._VAE_BETA_VAR.assign(float(value))


class Sampling(layers.Layer):
    def call(self, inputs, training=None):
        z_mean, z_log_var = inputs
        z_log_var = ops.clip(z_log_var, -10.0, 6.0)
        kl = -0.5 * ops.sum(
            1.0 + z_log_var - ops.square(z_mean) - ops.exp(z_log_var), axis=-1)
        beta = ops.cast(builtins._VAE_BETA_VAR, kl.dtype)
        self.add_loss(beta * ops.mean(kl))
        if not training:
            return z_mean
        eps = k_random.normal(shape=ops.shape(z_mean))
        return z_mean + ops.exp(0.5 * z_log_var) * eps

    def get_config(self):
        return super().get_config()


def masked_mse(y_true, y_pred):
    n = ops.shape(y_true)[-1] // 2
    spec = y_true[..., :n]
    mask = y_true[..., n:]
    diff2 = ops.square(y_pred - spec) * mask
    denom = ops.maximum(ops.sum(mask, axis=-1), 1.0)
    return ops.sum(diff2, axis=-1) / denom


class KLAnnealing(keras.callbacks.Callback):
    def __init__(self, target_beta, warmup_epochs):
        super().__init__()
        self.target_beta = float(target_beta)
        self.warmup_epochs = int(warmup_epochs)

    def on_train_begin(self, logs=None):
        _set_vae_beta(0.0)

    def on_epoch_begin(self, epoch, logs=None):
        if self.warmup_epochs > 0 and epoch < self.warmup_epochs:
            _set_vae_beta(self.target_beta * (epoch + 1) / self.warmup_epochs)
        else:
            _set_vae_beta(self.target_beta)


class Heartbeat(keras.callbacks.Callback):
    def __init__(self, every=10):
        super().__init__()
        self.every = every
        self.path = os.path.join(job_dir, "heartbeat")
        self.t0 = time.time()

    def _beat(self):
        with open(self.path, "w") as fh:
            fh.write(str(time.time()))

    def on_train_begin(self, logs=None):
        self._beat()

    def on_train_batch_end(self, batch, logs=None):
        if batch % self.every == 0:
            self._beat()

    def on_test_batch_end(self, batch, logs=None):
        if batch % self.every == 0:
            self._beat()

    def on_epoch_end(self, epoch, logs=None):
        self._beat()
        logs = logs or {}
        with open(os.path.join(job_dir, "progress.log"), "a") as fh:
            fh.write("epoch %d val_prediction_loss=%.4f %.0fs\n" % (
                epoch + 1, logs.get("val_prediction_loss", float("nan")),
                time.time() - self.t0))


def build(d_spec, latent_dim, l2_reg=1e-5, dropout_rate=0.2,
          encoder_dropout=0.05, loss_weights=(0.1, 1.0), learning_rate=1e-3):
    reg = lambda: keras.regularizers.L2(l2_reg)
    dense = lambda u, **kw: layers.Dense(u, activation="relu",
                                         kernel_regularizer=reg(), **kw)
    inp = keras.Input(shape=(2 * d_spec,))
    h = layers.BatchNormalization()(dense(512)(inp))
    h = layers.BatchNormalization()(dense(256)(h))
    h = layers.BatchNormalization()(dense(128)(h))
    h = layers.Dropout(encoder_dropout)(h)
    h = layers.BatchNormalization()(dense(64)(h))
    z_mean = layers.Dense(latent_dim, name="z_mean")(h)
    z_log_var = layers.Dense(latent_dim, name="z_log_var")(h)
    z = Sampling()([z_mean, z_log_var])
    r = layers.BatchNormalization()(dense(64)(z))
    r = layers.BatchNormalization()(dense(128)(r))
    r = layers.BatchNormalization()(dense(256)(r))
    rec = layers.Dense(d_spec, activation="linear", name="reconstruction")(r)
    p = dense(256)(z)
    p = layers.Dropout(dropout_rate)(p)
    p = dense(128)(p)
    p = layers.Dropout(dropout_rate / 2)(p)
    p = dense(64)(p)
    pred = layers.Dense(1, activation="linear", name="prediction")(p)
    mdl = keras.Model(inp, [rec, pred])
    mdl.compile(optimizer=keras.optimizers.Adam(learning_rate=learning_rate),
                loss=[masked_mse, "mse"], loss_weights=list(loss_weights))
    return mdl


def main():
    d = int(job["d_spec"])
    keras.utils.set_random_seed(int(job.get("seed", 42)))
    fit_x = _read("fit_x.f32", 2 * d)
    fit_y = _read("fit_y.f32", 1)
    val_x = _read("val_x.f32", 2 * d)
    val_y = _read("val_y.f32", 1)
    mdl = build(d, int(job["latent_dim"]))
    cbs = [
        keras.callbacks.EarlyStopping(
            monitor="val_prediction_loss", mode="min",
            patience=int(job["patience_es"]), restore_best_weights=True,
            start_from_epoch=int(job["kl_warmup_epochs"])),
        keras.callbacks.ReduceLROnPlateau(
            monitor="val_prediction_loss", mode="min",
            patience=int(job["patience_rl"]), factor=0.5, min_lr=1e-5),
        KLAnnealing(job["kl_beta"], job["kl_warmup_epochs"]),
        Heartbeat(),
    ]
    hist = mdl.fit(fit_x, [fit_x, fit_y],
                   validation_data=(val_x, [val_x, val_y]),
                   epochs=int(job["epochs"]), batch_size=int(job["batch_size"]),
                   callbacks=cbs, verbose=0)
    # Weights only: a whole .keras file names the torch optimizer, and loading
    # it would import torch into the R session. R rebuilds the same network
    # with build_maskedVAE() and checks it reproduces these predictions.
    mdl.save_weights(os.path.join(job_dir, "model.weights.h5"))
    check = mdl.predict(val_x[:256], verbose=0)[1].ravel().astype("<f4")
    check.tofile(os.path.join(job_dir, "check.f32"))
    with open(os.path.join(job_dir, "history.json"), "w") as fh:
        json.dump({k: [float(v) for v in vals]
                   for k, vals in hist.history.items()}, fh)
    with open(os.path.join(job_dir, "done"), "w") as fh:
        dev = str(mdl.weights[0].value.device) if keras.backend.backend() == "torch" else "cpu"
        fh.write("%s (%s)" % (keras.backend.backend(), dev))


def describe(d_spec, latent_dim):
    mdl = build(d_spec, latent_dim)
    print(json.dumps([[l.__class__.__name__, int(l.count_params())]
                      for l in mdl.layers]))


if __name__ == "__main__":
    if job_dir is None:
        describe(int(sys.argv[2]), int(sys.argv[3]))
    else:
        main()
