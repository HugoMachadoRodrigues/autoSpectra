# R/registry.R — Sensor / model-family registry for autoSpectra
#
# Each entry defines:
#   id              : unique string key (matches models/ subdirectory name)
#   label           : human-readable label for UI
#   sensor_type     : "visnir" | "mir"
#   sensors_allowed : character vector of specific sensor names (for local data)
#   moisture_levels : character vector ("DRY","1ML","3ML","agnostic")
#   wavegrid        : numeric vector of target wavelength/wavenumber positions (nm or cm-1)
#   preprocess      : character vector of pipeline step strings (see preprocess.R)
#   properties      : character vector of soil property keys this model predicts

# ---- OSSL Level-1 soil properties ----------------------------------------

#' All OSSL Level-1 harmonized soil property variable names
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

# ---- Local (legacy) soil properties --------------------------------------

local_soil_properties <- c(
  "soil_texture_sand", "soil_texture_silt", "soil_texture_clay",
  "organic_matter", "soc", "total_c", "total_n",
  "active_carbon", "ph", "p", "k",
  "mg", "fe", "mn", "zn",
  "al", "Ca", "Cu", "S",
  "B", "pred_soil_protein", "respiration", "bd_ws"
)

local_fancy_labels <- c(
  soil_texture_sand  = "Texture \u2014 Sand (%)",
  soil_texture_silt  = "Texture \u2014 Silt (%)",
  soil_texture_clay  = "Texture \u2014 Clay (%)",
  organic_matter     = "Organic Matter (%)",
  soc                = "Soil Organic Carbon (SOC, %)",
  total_c            = "Total C (%)",
  total_n            = "Total N (%)",
  active_carbon      = "Active Carbon (mg/kg)",
  ph                 = "pH",
  p                  = "Phosphorus (P, mg/kg)",
  k                  = "Potassium (K, mg/kg)",
  mg                 = "Magnesium (Mg, mg/kg)",
  fe                 = "Iron (Fe, mg/kg)",
  mn                 = "Manganese (Mn, mg/kg)",
  zn                 = "Zinc (Zn, mg/kg)",
  al                 = "Aluminum (Al, mg/kg)",
  Ca                 = "Calcium (Ca, mg/kg)",
  Cu                 = "Copper (Cu, mg/kg)",
  S                  = "Sulfur (S, mg/kg)",
  B                  = "Boron (B, mg/kg)",
  pred_soil_protein  = "Soil Protein (pred., mg/kg)",
  respiration        = "Respiration (\u00b5g CO\u2082-C g\u207b\u00b9 d\u207b\u00b9)",
  bd_ws              = "Bulk Density (g/cm\u00b3)"
)

# ---- Model family registry -----------------------------------------------

#' Full model family registry (local + OSSL families)
#'
#' Each family is a named list describing its sensor type, target wavegrid,
#' preprocessing pipeline, and predicted soil properties.
#'
#' @export
model_registry <- list(

  # --- Local / user-trained families (backward compatible) -----------------
  ASD_DRY = list(
    id              = "ASD_DRY",
    label           = "ASD \u2014 DRY (local, 23 props)",
    sensor_type     = "visnir",
    sensors_allowed = "ASD",
    moisture_levels = "DRY",
    wavegrid        = 350:2500,
    preprocess      = c("ABSORBANCE", "SG(11,2,1)"),
    properties      = local_soil_properties,
    source          = "local"
  ),

  NeoSpectra_DRY = list(
    id              = "NeoSpectra_DRY",
    label           = "NeoSpectra \u2014 DRY (local, 23 props)",
    sensor_type     = "visnir",
    sensors_allowed = "NeoSpectra",
    moisture_levels = "DRY",
    wavegrid        = 1350:2500,
    preprocess      = c("ABSORBANCE", "SG(11,2,1)"),
    properties      = local_soil_properties,
    source          = "local"
  ),

  NaturaSpec_DRY = list(
    id              = "NaturaSpec_DRY",
    label           = "NaturaSpec \u2014 DRY (local, 23 props)",
    sensor_type     = "visnir",
    sensors_allowed = "NaturaSpec",
    moisture_levels = "DRY",
    wavegrid        = 350:2500,
    preprocess      = c("ABSORBANCE", "SG(11,2,1)"),
    properties      = local_soil_properties,
    source          = "local"
  ),

  Agnostic_DRY = list(
    id              = "Agnostic_DRY",
    label           = "Agnostic \u2014 DRY (ASD + NaturaSpec + NeoSpectra, local)",
    sensor_type     = "visnir",
    sensors_allowed = c("ASD", "NaturaSpec", "NeoSpectra"),
    moisture_levels = "agnostic",
    wavegrid        = 1350:2500,
    preprocess      = c("ABSORBANCE", "SG(11,2,1)"),
    properties      = local_soil_properties,
    source          = "local"
  ),

  Agnostic_Moisture = list(
    id              = "Agnostic_Moisture",
    label           = "Agnostic \u2014 DRY+1ML+3ML (ASD + NaturaSpec + NeoSpectra, local)",
    sensor_type     = "visnir",
    sensors_allowed = c("ASD", "NaturaSpec", "NeoSpectra"),
    moisture_levels = "agnostic",
    wavegrid        = 1350:2500,
    preprocess      = c("ABSORBANCE", "SG(11,2,1)"),
    properties      = local_soil_properties,
    source          = "local"
  ),

  # --- OSSL-trained agnostic families (two-step SG) ------------------------
  OSSL_VisNIR = list(
    id              = "OSSL_VisNIR",
    label           = "OSSL \u2014 VisNIR Agnostic (all sensors, v1.2, 350-2500 nm)",
    sensor_type     = "visnir",
    sensors_allowed = NULL,   # accepts any VisNIR instrument
    moisture_levels = "agnostic",
    wavegrid        = seq(350, 2500, by = 2),   # 1076 channels
    preprocess      = c("ABSORBANCE", "SG_SMOOTH(11,2)", "SG_DERIV(11,2,1)"),
    properties      = ossl_l1_properties,
    source          = "ossl",
    ossl_version    = "v1.2",
    ossl_level      = "L1"
  ),

  OSSL_MIR = list(
    id              = "OSSL_MIR",
    label           = "OSSL \u2014 MIR Agnostic (all sensors, v1.2, 600-4000 cm\u207b\u00b9)",
    sensor_type     = "mir",
    sensors_allowed = NULL,   # accepts any MIR instrument
    moisture_levels = "agnostic",
    wavegrid        = seq(600, 4000, by = 2),   # 1701 channels (wavenumbers)
    preprocess      = c("SG_SMOOTH(11,2)", "SG_DERIV(11,2,1)"),  # already absorbance
    properties      = ossl_l1_properties,
    source          = "ossl",
    ossl_version    = "v1.2",
    ossl_level      = "L1"
  )
)

#' Lookup a model family from the registry
#'
#' @param id Family ID string (e.g., "OSSL_VisNIR")
#' @return List describing the family, or an error if not found
#' @export
get_family <- function(id) {
  fam <- model_registry[[id]]
  if (is.null(fam)) stop("Unknown model family: ", id)
  fam
}

#' Test whether a family matches a given sensor and moisture combination
#'
#' @param fam Family list from model_registry
#' @param sensor Character string for the sensor name
#' @param moisture Character string for the moisture mode
#' @return Logical
#' @export
family_matches <- function(fam, sensor, moisture) {
  ok_sensor  <- is.null(fam$sensors_allowed) || sensor %in% fam$sensors_allowed
  ok_moisture <- moisture %in% fam$moisture_levels || "agnostic" %in% fam$moisture_levels
  ok_sensor && ok_moisture
}

#' Get a display label for a soil property key
#'
#' Checks OSSL L1 labels first, then local labels, then returns the key itself.
#'
#' @param key Soil property key string
#' @return Human-readable label string
#' @export
property_label <- function(key) {
  lbl <- ossl_l1_labels[key]
  if (!is.na(lbl)) return(unname(lbl))
  lbl <- local_fancy_labels[key]
  if (!is.na(lbl)) return(unname(lbl))
  key
}
