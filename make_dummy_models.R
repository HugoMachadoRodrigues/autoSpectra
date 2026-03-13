# make_dummy_models.R
# Creates simple .rds models for each family/property so the app can predict.

# --- Minimal copy of your registry (edit if you changed it in app) ------------
default_registry <- list(
  ASD_DRY_full = list(
    id="ASD_DRY_full", wavegrid=350:2500,
    properties=c("B","Ca","Cu","S","Active Carbon","Al","BD","Fe","K","Mg","Mn",
                 "Organic Matter","P","Ph","Soil Protein","Respiration","SOC",
                 "Clay","Sand","Silt","Total C","Total N","Zn")
  ),
  Agnostic_Combined = list(
    id="Agnostic_Combined", wavegrid=1350:2500,
    properties=c("Active Carbon","Organic Matter","SOC","Total C")
  ),
  NeoSpectra_DRY = list(
    id="NeoSpectra_DRY", wavegrid=1350:2500,
    properties=c("Active Carbon","Organic Matter","SOC","Total C")
  ),
  NaturaSpec_DRY = list(
    id="NaturaSpec_DRY", wavegrid=350:2500,
    properties=c("Active Carbon","Organic Matter","SOC","Total C")
  ),
  NeoSpectra_Moisture_1ML = list(
    id="NeoSpectra_Moisture_1ML", wavegrid=1350:2500,
    properties=c("Active Carbon","Organic Matter","SOC","Total C")
  ),
  NeoSpectra_Moisture_3ML = list(
    id="NeoSpectra_Moisture_3ML", wavegrid=1350:2500,
    properties=c("Active Carbon","Organic Matter","SOC","Total C")
  ),
  NaturaSpec_Moisture_1ML = list(
    id="NaturaSpec_Moisture_1ML", wavegrid=350:2500,
    properties=c("Active Carbon","Organic Matter","SOC","Total C")
  ),
  NaturaSpec_Moisture_3ML = list(
    id="NaturaSpec_Moisture_3ML", wavegrid=350:2500,
    properties=c("Active Carbon","Organic Matter","SOC","Total C")
  )
)

# --- Tiny helpers -------------------------------------------------------------
dir_create <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

# Smooth-ish synthetic spectra with a few Gaussian bumps
sim_spectra <- function(n, wl) {
  set.seed(42)
  k <- 4L
  centers <- sample(wl, k, replace = FALSE)
  widths  <- sample(80:200, k, replace = TRUE)
  amps    <- runif(k, 0.5, 1.5)
  M <- matrix(0, n, length(wl))
  for (j in seq_len(k)) {
    M <- M + amps[j] * exp(-((rep(wl, each=n) - centers[j])^2)/(2*widths[j]^2))
  }
  # add small noise
  M <- M + matrix(rnorm(n*length(wl), sd = 0.02), n, length(wl))
  colnames(M) <- as.character(wl)
  M
}

# Pick ~30 wavelengths across the grid to keep models tiny & stable
pick_features <- function(wl, k = 30L) {
  if (length(wl) <= k) return(as.character(wl))
  idx <- round(seq(1, length(wl), length.out = k))
  as.character(wl[idx])
}

# Fit very small linear model on the selected wavelengths
fit_small_lm <- function(Xdf, y) {
  # Use a formula with only those ~30 columns; avoids p>>n issues
  nm <- names(Xdf)
  frm <- as.formula(paste("y ~", paste(sprintf("`%s`", nm), collapse = " + ")))
  lm(frm, data = cbind(Xdf, y = y))
}

# --- Main: build models for each family ---------------------------------------
build_family <- function(fam) {
  wl <- fam$wavegrid
  props <- fam$properties
  message("Building dummy models for family: ", fam$id,
          "  (grid ", min(wl), "-", max(wl), " nm, props: ", length(props), ")")
  
  outdir <- file.path("models", fam$id, "models")
  dir_create(outdir)
  
  n <- 80L                        # synthetic samples
  X <- sim_spectra(n, wl)         # n x length(wl)
  feats <- pick_features(wl, 30L) # ~30 wavelength columns
  Xdf <- as.data.frame(X[, feats, drop = FALSE], check.names = FALSE)
  
  # For each property, create a pseudo-target as a weighted sum + noise
  for (prop in props) {
    set.seed(abs(crc32 <- sum(utf8ToInt(prop))) %% 1e6)
    w <- rnorm(ncol(Xdf))
    y <- as.numeric(scale(as.matrix(Xdf) %*% w + rnorm(n, sd = 0.05)))
    
    mdl <- fit_small_lm(Xdf, y)
    
    # Save model
    saveRDS(mdl, file = file.path(outdir, paste0(prop, ".rds")))
  }
  invisible(TRUE)
}

# Build for all families in the list
invisible(lapply(default_registry, build_family))
cat("Done. Dummy models written under ./models/<FAMILY>/models/<Property>.rds\n")

