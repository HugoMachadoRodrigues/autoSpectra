# R/registry.R - Sensor / model-family registry for autoSpectra
#
# Four official families (two spectral domains x two algorithms):
#   OSSL_VisNIR        - masked soilVAE on 350-2500 nm (sensor-agnostic)
#   OSSL_MIR           - masked soilVAE on 600-4000 cm-1 (sensor-agnostic)
#   OSSL_VisNIR_Cubist - Cubist + PCA(120) + SNV (Safanelli 2025 replication)
#   OSSL_MIR_Cubist    - Cubist + PCA(120) + SNV (Safanelli 2025 replication)
#
# Spectral grids match Safanelli et al. (2025) PLoS ONE 20(1):e0296545:
#   VisNIR: 350-2500 nm at 2 nm   -> 1076 bands
#   MIR   : 600-4000 cm-1 at 2 cm-1 -> 1701 bands
#
# Preprocessing pipeline matches Safanelli (SNV alone). The VAE families
# additionally use a binary mask channel and per-mask SNV so partial-coverage
# instruments (e.g., NeoSpectra 1350-2550 nm) work without a separate model.

# ---- OSSL Level-1 soil properties ----------------------------------------

#' All OSSL Level-1 harmonized soil property variable names (34 targets)
#' @export
ossl_l1_properties <- c(
  "oc", "c.tot", "n.tot",
  "clay.tot", "silt.tot", "sand.tot",
  "ph.h2o", "ph.cacl2",
  "bd", "cec", "caco3",
  "al.ext", "al.ox", "al.dith",
  "ca.ext", "k.ext", "mg.ext", "na.ext",
  "p.ext",
  "fe.ext", "fe.dith", "fe.ox",
  "mn.ext", "zn.ext", "cu.ext", "b.ext",
  "s.tot", "s.ext",
  "acidity", "ec", "aggstb",
  "awc.33.1500kPa", "wr.33kPa", "wr.1500kPa"
)

#' Properties that Safanelli log1p-transforms before fitting
#'
#' Per Safanelli et al. 2025 (Methods, p.10): "Except for clay, silt, sand,
#' pH H2O, and pH CaCl2, all other soil properties were natural-log
#' transformed (with offset = 1, log1p() R function) to improve the
#' prediction performance of soil properties with highly skewed distribution."
#' @export
ossl_l1_log_transform <- setdiff(
  ossl_l1_properties,
  c("clay.tot", "silt.tot", "sand.tot", "ph.h2o", "ph.cacl2")
)

#' Fancy display labels for OSSL L1 properties
#' @export
ossl_l1_labels <- c(
  oc              = "Organic Carbon (OC, %)",
  c.tot           = "Total Carbon (%)",
  n.tot           = "Total Nitrogen (%)",
  clay.tot        = "Clay (%)",
  silt.tot        = "Silt (%)",
  sand.tot        = "Sand (%)",
  ph.h2o          = "pH (H\u2082O)",
  ph.cacl2        = "pH (CaCl\u2082)",
  bd              = "Bulk Density (g/cm\u00b3)",
  cec             = "CEC (cmolc/kg)",
  caco3           = "CaCO\u2083 (%)",
  al.ext          = "Al extractable (mg/kg)",
  al.ox           = "Al amorphous (w.pct)",
  al.dith         = "Al crystalline (w.pct)",
  ca.ext          = "Ca extractable (mg/kg)",
  k.ext           = "K extractable (mg/kg)",
  mg.ext          = "Mg extractable (mg/kg)",
  na.ext          = "Na extractable (mg/kg)",
  p.ext           = "P extractable (mg/kg)",
  fe.ext          = "Fe extractable (mg/kg)",
  fe.dith         = "Fe crystalline (w.pct)",
  fe.ox           = "Fe amorphous (w.pct)",
  mn.ext          = "Mn extractable (mg/kg)",
  zn.ext          = "Zn extractable (mg/kg)",
  cu.ext          = "Cu extractable (mg/kg)",
  b.ext           = "B extractable (mg/kg)",
  s.tot           = "Total Sulfur (%)",
  s.ext           = "S extractable (mg/kg)",
  acidity         = "Exchangeable Acidity (cmolc/kg)",
  ec              = "Electrical Conductivity (dS/m)",
  aggstb          = "Aggregate Stability (%)",
  awc.33.1500kPa  = "Available Water Content (w.frac)",
  wr.33kPa        = "Water Retention 33 kPa (%)",
  wr.1500kPa      = "Water Retention 1500 kPa (%)"
)

# ---- Instrument metadata -------------------------------------------------

#' VisNIR instruments contributing to OSSL v1.2
#' @export
ossl_visnir_instruments <- c(
  "ASD FieldSpec 3 / 4 (Malvern Panalytical) - 350-2500 nm",
  "Foss XDS Rapid Content Analyzer - 400-2500 nm",
  "Bruker MPA FT-NIR - 800-2500 nm",
  "PerkinElmer Spectrum One FT-NIR - 1000-2500 nm",
  "Si-Ware NeoSpectra Handheld - 1350-2550 nm",
  "Other diffuse-reflectance VisNIR (>=350 nm)"
)

#' MIR instruments contributing to OSSL v1.2
#' @export
ossl_mir_instruments <- c(
  "Bruker ALPHA FTIR (ATR diamond crystal)",
  "Bruker Tensor 27 FTIR (DRIFTS)",
  "Bruker Vertex 70 (Pike Autodiff)",
  "PerkinElmer Spectrum Two FT-MIR",
  "Thermo Nicolet FT-MIR (various)",
  "Other FTIR instruments (600-4000 cm-1)"
)

# ---- Model family registry -----------------------------------------------

# Common preprocessing pipelines (match Safanelli 2025)
.SAFANELLI_VISNIR_PP <- c("ABSORBANCE", "SNV")
.SAFANELLI_MIR_PP    <- c("SNV")

# VAE families add SG smoothing before SNV to denoise the input;
# SNV is then computed per-sample using only the masked (covered) bands.
.VAE_VISNIR_PP       <- c("ABSORBANCE", "SG_SMOOTH(11,2)", "SNV_MASKED")
.VAE_MIR_PP          <- c("SG_SMOOTH(11,2)", "SNV_MASKED")

#' Official autoSpectra model family registry
#'
#' Four sensor-agnostic families trained on the full OSSL v1.2 corpus.
#' Two model types per spectral domain:
#'   - soilVAE (masked, sensor-agnostic) - production models
#'   - Cubist + PCA(120) (Safanelli 2025 replication) - benchmark models
#'
#' All families share the same OSSL spectral grids:
#'   VisNIR: 350-2500 nm at 2 nm
#'   MIR   : 600-4000 cm-1 at 2 cm-1
#'
#' @export
model_registry <- list(

  OSSL_VisNIR = list(
    id              = "OSSL_VisNIR",
    label           = "OSSL VisNIR \u2014 masked soilVAE (350-2500 nm, any instrument)",
    sensor_type     = "visnir",
    algorithm       = "vae_masked",
    sensors_allowed = NULL,
    moisture_levels = "agnostic",
    wavegrid        = seq(350, 2500, by = 2),
    preprocess      = .VAE_VISNIR_PP,
    properties      = ossl_l1_properties,
    source          = "ossl",
    ossl_version    = "v1.2",
    ossl_level      = "L1",
    instruments     = ossl_visnir_instruments,
    citation        = "Safanelli et al. (2025) doi:10.1371/journal.pone.0296545"
  ),

  OSSL_MIR = list(
    id              = "OSSL_MIR",
    label           = "OSSL MIR \u2014 masked soilVAE (600-4000 cm\u207b\u00b9, any instrument)",
    sensor_type     = "mir",
    algorithm       = "vae_masked",
    sensors_allowed = NULL,
    moisture_levels = "agnostic",
    wavegrid        = seq(600, 4000, by = 2),
    preprocess      = .VAE_MIR_PP,
    properties      = ossl_l1_properties,
    source          = "ossl",
    ossl_version    = "v1.2",
    ossl_level      = "L1",
    instruments     = ossl_mir_instruments,
    citation        = "Safanelli et al. (2025) doi:10.1371/journal.pone.0296545"
  ),

  OSSL_VisNIR_Cubist = list(
    id              = "OSSL_VisNIR_Cubist",
    label           = "OSSL VisNIR \u2014 Cubist + PCA(120) [Safanelli 2025 replication]",
    sensor_type     = "visnir",
    algorithm       = "cubist_pca",
    sensors_allowed = NULL,
    moisture_levels = "agnostic",
    wavegrid        = seq(350, 2500, by = 2),
    preprocess      = .SAFANELLI_VISNIR_PP,
    properties      = ossl_l1_properties,
    source          = "ossl",
    ossl_version    = "v1.2",
    ossl_level      = "L1",
    instruments     = ossl_visnir_instruments,
    citation        = "Safanelli et al. (2025) doi:10.1371/journal.pone.0296545"
  ),

  OSSL_MIR_Cubist = list(
    id              = "OSSL_MIR_Cubist",
    label           = "OSSL MIR \u2014 Cubist + PCA(120) [Safanelli 2025 replication]",
    sensor_type     = "mir",
    algorithm       = "cubist_pca",
    sensors_allowed = NULL,
    moisture_levels = "agnostic",
    wavegrid        = seq(600, 4000, by = 2),
    preprocess      = .SAFANELLI_MIR_PP,
    properties      = ossl_l1_properties,
    source          = "ossl",
    ossl_version    = "v1.2",
    ossl_level      = "L1",
    instruments     = ossl_mir_instruments,
    citation        = "Safanelli et al. (2025) doi:10.1371/journal.pone.0296545"
  )
)

#' Lookup a model family from the registry
#'
#' @param id Family ID string
#' @return Named list describing the family
#' @export
get_family <- function(id) {
  fam <- model_registry[[id]]
  if (is.null(fam)) stop("Unknown model family: '", id,
    "'. Available: ", paste(names(model_registry), collapse = ", "))
  fam
}

#' Test whether a family matches a sensor type
#'
#' @param fam Family list from \code{model_registry}
#' @param sensor_type Character: \code{"visnir"} or \code{"mir"}
#' @param moisture Ignored; retained for backward compatibility
#' @return Logical
#' @export
family_matches <- function(fam, sensor_type, moisture = "agnostic") {
  identical(fam$sensor_type, sensor_type)
}

#' Get a display label for a soil property key
#'
#' @param key Soil property key string (e.g. \code{"oc"})
#' @return Human-readable label string
#' @export
property_label <- function(key) {
  lbl <- ossl_l1_labels[key]
  if (!is.na(lbl)) return(unname(lbl))
  key
}
