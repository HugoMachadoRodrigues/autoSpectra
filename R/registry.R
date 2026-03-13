# R/registry.R — Sensor / model-family registry for autoSpectra
#
# Only two official families:
#   OSSL_VisNIR — all VisNIR instruments in OSSL v1.2, 400–2500 nm
#   OSSL_MIR   — all MIR instruments in OSSL v1.2, 600–4000 cm-1
#
# NOTE on VisNIR range: The OSSL v1.2 corpus has 64k samples with coverage
# from 400 nm, but only 24k with coverage from 350 nm (instruments lacking the
# 350-400 nm range, e.g., Foss XDS, many benchtop NIR).  Using 400-2500 nm
# triples the training corpus while losing negligible soil-relevant signal
# (the 350-400 nm Vis region contributes little to soil property prediction).
#
# Both are sensor-agnostic: trained on the full OSSL corpus (30+ datasets,
# dozens of instruments). The two-step SG first-derivative preprocessing
# eliminates instrument-specific baseline and multiplicative-scatter
# differences without discarding spectral information.

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
  "ASD FieldSpec 3 / 4 (Malvern Panalytical)",
  "Foss XDS Rapid Content Analyzer",
  "Bruker MPA FT-NIR",
  "PerkinElmer Spectrum One FT-NIR",
  "Vis-NIR (400-2500nm, other instruments)"
)

#' MIR instruments contributing to OSSL v1.2
#' @export
ossl_mir_instruments <- c(
  "Bruker ALPHA FTIR (ATR diamond crystal)",
  "Bruker Tensor 27 FTIR (DRIFTS)",
  "PerkinElmer Spectrum Two FT-MIR",
  "Thermo Nicolet FT-MIR (various)",
  "MIR (600-4000 cm-1, other instruments)"
)

# ---- Model family registry -----------------------------------------------

#' Official autoSpectra model family registry
#'
#' Two sensor-agnostic families trained on the full OSSL v1.2 corpus.
#' Use \code{OSSL_VisNIR} for any diffuse-reflectance VisNIR instrument
#' (400-2500 nm) and \code{OSSL_MIR} for any FTIR/ATR/DRIFTS MIR instrument
#' (600-4000 cm-1).
#'
#' @export
model_registry <- list(

  OSSL_VisNIR = list(
    id              = "OSSL_VisNIR",
    label           = "OSSL VisNIR \u2014 Agnostic (all instruments, 400-2500 nm)",
    sensor_type     = "visnir",
    sensors_allowed = NULL,      # accepts any VisNIR instrument via resampling
    moisture_levels = "agnostic",
    wavegrid        = seq(400, 2500, by = 2),   # 1051 channels at 2 nm
    preprocess      = c("ABSORBANCE", "SG_SMOOTH(11,2)", "SG_DERIV(11,2,1)"),
    properties      = ossl_l1_properties,
    source          = "ossl",
    ossl_version    = "v1.2",
    ossl_level      = "L1",
    instruments     = ossl_visnir_instruments,
    citation        = "Safanelli et al. (2023) doi:10.5194/essd-15-3829-2023"
  ),

  OSSL_MIR = list(
    id              = "OSSL_MIR",
    label           = "OSSL MIR \u2014 Agnostic (all instruments, 600-4000 cm\u207b\u00b9)",
    sensor_type     = "mir",
    sensors_allowed = NULL,      # accepts any MIR instrument via resampling
    moisture_levels = "agnostic",
    wavegrid        = seq(600, 4000, by = 2),   # 1701 channels at 2 cm-1
    preprocess      = c("SG_SMOOTH(11,2)", "SG_DERIV(11,2,1)"),
    properties      = ossl_l1_properties,
    source          = "ossl",
    ossl_version    = "v1.2",
    ossl_level      = "L1",
    instruments     = ossl_mir_instruments,
    citation        = "Safanelli et al. (2023) doi:10.5194/essd-15-3829-2023"
  )
)

#' Lookup a model family from the registry
#'
#' @param id Family ID string: \code{"OSSL_VisNIR"} or \code{"OSSL_MIR"}
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
