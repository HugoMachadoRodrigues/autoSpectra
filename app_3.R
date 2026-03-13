# app.R — autoSpectra (Horizontal Sophisticated Layout)

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(shinyWidgets)
  library(readxl)
  library(readr)
  library(DT)
  library(ggplot2)
  library(prospectr)  # SG filter
  library(writexl)
})

options(shiny.fullstacktrace = TRUE, shiny.maxRequestSize = 200*1024^2)

# Read metrics JSON sidecar
read_metrics_json <- function(family_id, prop) {
  fp <- file.path("models", family_id, "models", paste0(prop, "_metrics.json"))
  if (file.exists(fp)) jsonlite::read_json(fp, simplifyVector = TRUE) else NULL
}

# Robustly coerce latent$Sigma from _metrics.json into a square matrix
coerce_sigma <- function(Sigma_json, d_expected = NULL) {
  # Case A: already matrix/data.frame
  if (is.matrix(Sigma_json)) return(Sigma_json)
  if (is.data.frame(Sigma_json)) return(as.matrix(Sigma_json))
  
  # Case B: list of rows (what we saved originally as asplit(..., 1L))
  if (is.list(Sigma_json)) {
    # ensure each element is a numeric vector
    rows <- lapply(Sigma_json, function(x) as.numeric(unlist(x, use.names = FALSE)))
    # rbind safely
    mat <- do.call(rbind, rows)
    # if it's 1 x n (single row), try to reshape later
    if (nrow(mat) == 1L && !is.null(d_expected) && length(mat) == d_expected * d_expected) {
      mat <- matrix(as.numeric(mat), nrow = d_expected, byrow = TRUE)
    }
    return(as.matrix(mat))
  }
  
  # Case C: numeric vector (flattened)
  vec <- as.numeric(Sigma_json)
  if (!is.null(d_expected) && length(vec) == d_expected * d_expected) {
    return(matrix(vec, nrow = d_expected, byrow = TRUE))
  }
  
  stop("Could not coerce latent covariance 'Sigma' to a square matrix.")
}

# Robust inverse via chol; small ridge in case Σ is near-singular
inv_pd <- function(Sigma, ridge = 1e-6) {
  S <- as.matrix(Sigma)
  d <- ncol(S); S <- S + diag(ridge, d)
  chol2inv(chol(S))
}

# % of wavelengths inside training min/max after preprocessing
applicability_pct <- function(x_row, feat_range_df) {
  common <- intersect(names(x_row), feat_range_df$band)
  if (!length(common)) return(NA_real_)
  fr <- feat_range_df[match(common, feat_range_df$band), ]
  vals <- as.numeric(x_row[common])
  inside <- (vals >= fr$min) & (vals <= fr$max)
  100 * mean(inside, na.rm = TRUE)
}

# ------------------------------------------------------------------------------
# THEME & CSS
# ------------------------------------------------------------------------------
theme <- bs_theme(version = 5, bootswatch = "lux",
                  base_font = font_google("Inter"),
                  heading_font = font_google("Inter Tight"),
                  primary = "#2B6CB0", secondary = "#6B46C1")

app_css <- HTML("
  .app-header {
    position: sticky; top: 0; z-index: 1030;
    background: #ffffffE6; backdrop-filter: blur(6px);
    border-bottom: 1px solid rgba(0,0,0,0.06);
  }
  .toolbar {
  display: flex; flex-wrap: wrap;
  gap: 12px; align-items: center;
  padding: 10px 14px;
}
.toolbar .inputs, .toolbar .actions {
  display: flex; flex-wrap: wrap;
  gap: 10px; align-items: center;
}
  .brand {
    display: flex; align-items: center; gap: 12px;
  }
  .brand .title {
    font-weight: 700; font-size: 1.25rem; letter-spacing: .2px;
  }
  .pill {
    background: #f8f9fa; border-radius: 999px; padding: 10px 14px;
    box-shadow: inset 0 0 0 1px rgba(0,0,0,0.06);
  }
  .toolbar .inputs {
    display: grid; grid-auto-flow: column; gap: 10px; align-items: center;
  }
  .toolbar .actions {
    display: grid; grid-auto-flow: column; gap: 10px; align-items: center;
  }
  .cardish {
    background: #fff; border-radius: 16px; padding: 16px;
    box-shadow: 0 4px 18px rgba(0,0,0,0.05), 0 1px 2px rgba(0,0,0,0.05);
    border: 1px solid rgba(0,0,0,0.05);
  }
  .muted { color: #6c757d; }
  .tab-title { margin: 4px 0 14px 0; font-weight: 600; }
  .spacer { height: 10px; }
  .form-label { margin-bottom: 4px; }
  .btn-primary { border-radius: 999px; }
  .btn-success { border-radius: 999px; }
  .shiny-output-error { color: #B00020; }
")

# ------------------------------------------------------------------------------
# VALIDATOR
# ------------------------------------------------------------------------------
vneed <- function(condition, message = "Validation failed.") {
  msg <- tryCatch({
    m <- paste(as.character(message), collapse = " ")
    if (!nzchar(m)) "Validation failed." else m
  }, error = function(e) "Validation failed.")
  shiny::validate(shiny::need(condition, msg))
}

safe_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), x, NA_real_)
}

# ------------------------------------------------------------------------------
# TARGETS & FANCY LABELS
# ------------------------------------------------------------------------------
soil_properties <- c(
  "soil_texture_sand","soil_texture_silt","soil_texture_clay",
  "organic_matter","soc","total_c","total_n","active_carbon",
  "ph","p","k","mg","fe","mn","zn","al","Ca","Cu","S","B",
  "pred_soil_protein","respiration","bd_ws"
)

property_order <- c(
  "soil_texture_sand","soil_texture_silt","soil_texture_clay",
  "organic_matter","soc","total_c","total_n","active_carbon",
  "ph","p","k","mg","fe","mn","zn","al","Ca","Cu","S","B",
  "pred_soil_protein","respiration","bd_ws"
)

fancy_labels <- c(
  soil_texture_sand  = "Texture — Sand (%)",
  soil_texture_silt  = "Texture — Silt (%)",
  soil_texture_clay  = "Texture — Clay (%)",
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
  respiration        = "Respiration (µg CO₂-C g⁻¹ d⁻¹)",
  bd_ws              = "Bulk Density (g/cm³)"
)
get_label <- function(key) {
  lbl <- unname(fancy_labels[key])
  if (is.na(lbl)) key else lbl
}

# ------------------------------------------------------------------------------
# MODEL REGISTRY
# ------------------------------------------------------------------------------
default_registry <- list(
  ASD_DRY = list(
    id = "ASD_DRY",
    label = "ASD — DRY",
    sensors_allowed = c("ASD"),
    moisture_levels = c("DRY"),
    properties = soil_properties,
    wavegrid = 350:2500,
    preprocess = c("ABSORBANCE","SG(11,2,1)")
  ),
  NeoSpectra_DRY = list(
    id = "NeoSpectra_DRY",
    label = "NeoSpectra — DRY",
    sensors_allowed = c("NeoSpectra"),
    moisture_levels = c("DRY"),
    properties = soil_properties,
    wavegrid = 1350:2500,
    preprocess = c("ABSORBANCE","SG(11,2,1)")
  ),
  NaturaSpec_DRY = list(
    id = "NaturaSpec_DRY",
    label = "NaturaSpec — DRY",
    sensors_allowed = c("NaturaSpec"),
    moisture_levels = c("DRY"),
    properties = soil_properties,
    wavegrid = 350:2500,
    preprocess = c("ABSORBANCE","SG(11,2,1)")
  ),
  Agnostic_DRY = list(
    id = "Agnostic_DRY",
    label = "Agnostic — DRY (All Sensors)",
    sensors_allowed = c("ASD","NaturaSpec","NeoSpectra"),
    moisture_levels = c("agnostic"),
    properties = soil_properties,
    wavegrid = 1350:2500,
    preprocess = c("ABSORBANCE","SG(11,2,1)")
  ),
  Agnostic_Moisture = list(
    id = "Agnostic_Moisture",
    label = "Agnostic — DRY + 1ML + 3ML (All Sensors)",
    sensors_allowed = c("ASD","NaturaSpec","NeoSpectra"),
    moisture_levels = c("agnostic"),
    properties = soil_properties,
    wavegrid = 1350:2500,
    preprocess = c("ABSORBANCE","SG(11,2,1)")
  )
)

family_matches <- function(fam, sensor, moisture) {
  ok_sensor <- sensor %in% fam$sensors_allowed
  ok_moist  <- moisture %in% fam$moisture_levels || ("agnostic" %in% fam$moisture_levels)
  ok_sensor && ok_moist
}

# ------------------------------------------------------------------------------
# WAVELENGTHS & PREPROCESSING
# ------------------------------------------------------------------------------
get_wavelengths <- function(df) {
  cols <- setdiff(names(df), "Soil_ID")
  cols <- as.character(cols)
  wl   <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", cols)))
  idx  <- !is.na(wl)
  list(wl = wl[idx], cols = cols[idx])
}

resample_to_grid <- function(df, src_wl, target_wl) {
  M <- as.matrix(df)
  out <- matrix(NA_real_, nrow = nrow(M), ncol = length(target_wl))
  for (i in seq_len(nrow(M))) {
    out[i, ] <- approx(x = src_wl, y = M[i, ], xout = target_wl, rule = 1, ties = mean)$y
  }
  colnames(out) <- as.character(target_wl)
  out
}

# Reflectance -> Absorbance (-ln R); auto-scale percentages to [0,1]
reflectance_to_absorbance <- function(M, base10 = FALSE) {
  Mnum <- as.matrix(M)
  if (is.finite(max(Mnum, na.rm = TRUE)) && max(Mnum, na.rm = TRUE) > 2) {
    Mnum <- Mnum / 100
  }
  eps <- 1e-6
  Mnum[Mnum < eps] <- eps
  if (base10) -log10(Mnum) else -log(Mnum)
}

apply_snv <- function(M) {
  M_centered <- sweep(M, 1, rowMeans(M, na.rm = TRUE), "-")
  denom <- sqrt(rowMeans(M_centered^2, na.rm = TRUE))
  denom[denom == 0 | is.na(denom)] <- 1
  sweep(M_centered, 1, denom, "/")
}

parse_sg <- function(st) {
  inside <- substr(st, 4, nchar(st) - 1)   # "m,p" or "m,p,d"
  parts  <- strsplit(inside, ",", fixed = TRUE)[[1]]
  parts  <- trimws(parts)
  m <- suppressWarnings(as.integer(parts[1])); if (is.na(m)) m <- 11
  p <- suppressWarnings(as.integer(parts[2])); if (is.na(p)) p <- 2
  d <- if (length(parts) >= 3) suppressWarnings(as.integer(parts[3])) else 0
  if (is.na(d)) d <- 0
  list(m = m, p = p, d = d)
}
sg_safe_params <- function(m, p, ncols) {
  if (!is.finite(m) || m < 3) m <- 3
  if (!is.finite(p) || p < 0) p <- 2
  if (m %% 2 == 0) m <- m + 1
  min_req <- p + 1
  if (m < min_req) m <- min_req + ifelse(min_req %% 2 == 0, 1, 0)
  if (ncols < 3) return(NULL)
  if (m > ncols) m <- ncols - ifelse(ncols %% 2 == 0, 1, 0)
  if (m < 3) return(NULL)
  list(m = m, p = p)
}
apply_sg_matrix <- function(M, m, p, d = 0) {
  ncols <- ncol(M)
  pars  <- sg_safe_params(m, p, ncols)
  if (is.null(pars)) return(M)
  m <- pars$m; p <- pars$p
  res <- t(apply(M, 1, function(r) {
    if (all(is.na(r))) return(r)
    tryCatch(prospectr::savitzkyGolay(r, m = m, p = p, w = d),
             error = function(e) r)
  }))
  colnames(res) <- colnames(M)
  res
}

apply_pipeline <- function(M, steps, absorbance_base10 = FALSE) {
  out <- M
  if (length(steps) == 0) return(out)
  for (st in steps) {
    if (identical(st, "ABSORBANCE")) {
      out <- reflectance_to_absorbance(out, base10 = absorbance_base10)
    } else if (startsWith(st, "SG(")) {
      sg <- parse_sg(st)
      out <- apply_sg_matrix(out, m = sg$m, p = sg$p, d = sg$d)
    } else if (identical(st, "SNV")) {
      out <- apply_snv(out)
    }
  }
  out
}

# ------------------------------------------------------------------------------
# UI
# ------------------------------------------------------------------------------
ui <- page_fluid(
  theme = theme,
  tags$head(tags$style(app_css),
            tags$link(rel="icon", type="image/png", href="logo.png")),
  div(class="app-header",
      div(class="toolbar container-fluid",
          div(class="brand",
              img(src="logo.png", height="210px"),
              span(class="title", "")
          ),
          div(class="inputs",
              div(class="pill",
                  pickerInput("sensor", NULL,
                              choices = c("ASD","NeoSpectra","NaturaSpec"),
                              selected = "ASD", width = "180px")
              ),
              div(class="pill",
                  pickerInput("moisture", NULL,
                              choices = c("DRY","1ML","3ML","agnostic"),
                              selected = "DRY", width="140px")
              ),
              div(class="pill",
                  pickerInput("family", NULL,
                              choices = setNames(names(default_registry),
                                                 vapply(default_registry, `[[`, "", "label")),
                              selected = "ASD_DRY", width="260px")
              ),
              div(class="pill",
                  uiOutput("props_ui", inline=TRUE)
              )
          ),
          div(class="actions",
              fileInput("file", NULL, accept = c(".xlsx",".xls",".csv"), buttonLabel = "Upload data", placeholder = "Excel/CSV"),
              actionButton("preview", "Preview", class="btn btn-primary"),
              actionButton("run", "Predict", class="btn btn-success"),
              checkboxInput("disable_pp", "Raw", value = FALSE, width = "80px")
          )
      )
  ),
  div(class="container my-3",
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Analysis",
                 fluidRow(
                   column(12,
                          div(class="cardish",
                              h5(class="tab-title", "Data Preview"),
                              uiOutput("sheet_ui"),
                              div(class="spacer"),
                              verbatimTextOutput("info", placeholder = TRUE),
                              DTOutput("head_dt")
                          )
                   )
                 ),
                 div(class="spacer"),
                 fluidRow(
                   column(6,
                          div(class="cardish",
                              h5(class="tab-title", "Spectrum"),
                              uiOutput("soil_pick_ui"),
                              plotOutput("spec_plot", height = "320px"),
                              span(class="muted", textOutput("coverage"))
                          )
                   ),
                   column(6,
                          div(class="cardish",
                              h5(class="tab-title", "Predictions"),
                              DTOutput("pred_dt"),
                              div(class="spacer"),
                              DTOutput("model_info_dt"),
                              div(class="spacer"),
                              downloadButton("dl_metrics", "Download model metrics"),
                              div(class="spacer"),
                              downloadButton("dl_preds", "Download predictions")
                          )
                   )
                 )
        ),
        tabPanel("Help",
                 div(class="cardish",
                     h4("About autoSpectra"),
                     p("autoSpectra predicts soil properties from spectral readings collected with different sensors (ASD, NaturaSpec, NeoSpectra)."),
                     p("Models are trained with pretreatments (Absorbance + Savitzky–Golay smoothing) and saved separately for:"),
                     tags$ul(
                       tags$li("ASD (350–2500 nm, DRY samples)"),
                       tags$li("NaturaSpec (350–2500 nm, DRY and moist 1ML/3ML samples)"),
                       tags$li("NeoSpectra (1350–2500 nm, DRY and moist 1ML/3ML samples)"),
                       tags$li("Agnostic families (intersection grid 1350–2500 nm, combining sensors)")
                     ),
                     h4("Predicted Soil Properties"),
                     p("Depending on the family, models can predict texture (sand, silt, clay), soil chemistry (pH, C, N, P, K, Mg, Fe, Mn, Zn, Al, Ca, Cu, S, B), organic indicators (Organic Matter, SOC, Active Carbon, Total C, Total N), and biological proxies (Soil Protein, Respiration, Bulk Density)."),
                     h4("Interpreting Prediction Outputs"),
                     tags$ul(
                       tags$li(strong("MC sd:"), " Monte Carlo dropout standard deviation. Captures model uncertainty due to limited training data or complex relationships. Larger MC sd means the model is less certain about the prediction."),
                       tags$li(strong("PI95 low / PI95 high:"), " 95% conformal prediction interval. The range within which the true value is expected to lie with 95% probability, based on calibration residuals."),
                       tags$li(strong("Latent App. %:"), " Latent applicability score. Derived from the Mahalanobis distance in the latent space. High values (>80%) mean the sample is similar to the training set. Low values (<20%) suggest extrapolation or out-of-distribution input.")
                     ),
                     h4("Practical Use"),
                     p("For each soil sample, examine not just the predicted values, but also the uncertainty (MC sd and PI95) and the applicability %. Predictions with low applicability or high MC sd should be interpreted cautiously.")
                 )
        )
      )
  )
)

# ------------------------------------------------------------------------------
# SERVER
# ------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Filter family choices by sensor & moisture
  observe({
    req(input$sensor, input$moisture)
    filt_ids <- names(default_registry)[vapply(
      default_registry,
      function(f) family_matches(f, as.character(input$sensor), as.character(input$moisture)),
      logical(1)
    )]
    if (length(filt_ids) == 0) filt_ids <- names(default_registry)
    updatePickerInput(session, "family",
                      choices = setNames(filt_ids, vapply(default_registry[filt_ids], `[[`, "", "label")),
                      selected = if (input$family %in% filt_ids) input$family else filt_ids[1])
  })
  
  fam <- reactive({ default_registry[[ as.character(input$family) ]] })
  
  # Properties picker with fancy labels, horizontal friendly
  output$props_ui <- renderUI({
    f <- fam(); req(f)
    ordered <- intersect(property_order, f$properties)
    choices <- setNames(ordered, vapply(ordered, get_label, FUN.VALUE = character(1)))
    pickerInput("props", NULL, choices = choices, selected = ordered, multiple = TRUE,
                options = pickerOptions(actionsBox = TRUE, liveSearch = TRUE, dropupAuto = FALSE, size = 10),
                width = "320px")
  })
  
  # Sheet selector for Excel
  output$sheet_ui <- renderUI({
    req(input$file)
    ext <- tolower(tools::file_ext(input$file$name))
    if (ext %in% c("xlsx","xls")) {
      path <- normalizePath(as.character(input$file$datapath), winslash = "/")
      shs  <- as.character(readxl::excel_sheets(path))
      selectInput("sheet", "Sheet", choices = shs, selected = shs[1], width = "240px")
    }
  })
  
  # Read data on Preview
  raw <- eventReactive(input$preview, {
    req(input$file)
    ext  <- tolower(tools::file_ext(input$file$name))
    path <- normalizePath(as.character(input$file$datapath), winslash = "/")
    if (ext %in% c("xlsx","xls")) {
      sh  <- if (isTruthy(input$sheet)) as.character(input$sheet) else 1
      df  <- readxl::read_excel(path, sheet = sh)
    } else {
      df  <- readr::read_csv(path, show_col_types = FALSE)
    }
    df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
    vneed("Soil_ID" %in% names(df), "Column 'Soil_ID' not found. Please include a Soil_ID column.")
    df
  }, ignoreInit = TRUE)
  
  # Preview info
  output$info <- renderPrint({
    req(raw())
    df <- raw(); f <- fam()
    wl <- get_wavelengths(df)
    cat("Rows x Cols:", nrow(df), "x", ncol(df), "\n")
    if (length(wl$wl) > 0) {
      cat("Detected spectral columns:", length(wl$wl),
          sprintf(" (%.0f-%.0f nm)\n", min(wl$wl), max(wl$wl)))
    } else {
      cat("Detected spectral columns: 0\n")
    }
    cat("Family grid:", length(f$wavegrid),
        sprintf(" (%.0f-%.0f nm)\n", min(f$wavegrid), max(f$wavegrid)))
    cat("Preprocess:", paste(f$preprocess, collapse = " → "),
        if (isTRUE(input$disable_pp)) " [DISABLED]" else "", "\n")
  })
  
  output$head_dt <- renderDT({
    req(raw())
    DT::datatable(head(raw(), 8), options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  # Spectrum viewer
  output$soil_pick_ui <- renderUI({
    req(raw())
    selectInput("soil_pick", "Sample", choices = raw()$Soil_ID, width = "240px")
  })
  
  output$spec_plot <- renderPlot({
    req(raw(), input$soil_pick)
    df <- raw(); f <- fam(); wl <- get_wavelengths(df)
    row <- df[df$Soil_ID == input$soil_pick, , drop = FALSE]
    vals <- as.numeric(row[1, wl$cols, drop = TRUE])
    dd <- data.frame(wl = wl$wl, refl = vals)
    ggplot(dd, aes(wl, refl)) +
      geom_line() +
      geom_vline(xintercept = range(f$wavegrid), linetype = "dashed") +
      labs(x = "Wavelength (nm)", y = "Response",
           subtitle = sprintf("Model grid %.0f-%.0f nm", min(f$wavegrid), max(f$wavegrid))) +
      theme_minimal(base_size = 12)
  })
  
  output$coverage <- renderText({
    req(raw())
    df <- raw(); f <- fam(); wl <- get_wavelengths(df)
    overlap <- sum(wl$wl %in% f$wavegrid)
    sprintf("Header overlap with model grid: %.1f%% (%d/%d). Missing bands are linearly interpolated.",
            100*overlap/length(f$wavegrid), overlap, length(f$wavegrid))
  })
  
  output$model_info <- renderPrint({
    req(fam(), input$props)
    # Pick the first selected property for summary (or make a dedicated picker)
    prop <- input$props[[1]]
    mets <- read_metrics_json(fam()$id, prop)
    if (is.null(mets)) {
      cat("Model metrics not found for", prop, "\nExpected:", file.path("models", fam()$id, "models", paste0(prop, "_metrics.json")))
      return(invisible(NULL))
    }
    cat("Model:", fam()$label, "\nProperty:", get_label(prop), "\n\n")
    cat(sprintf("Test RMSE : %.3f\n", as.numeric(mets$metrics$RMSE)))
    cat(sprintf("Test R²   : %.3f\n", as.numeric(mets$metrics$R2)))
    cat(sprintf("Test RPIQ : %.3f\n", as.numeric(mets$metrics$RPIQ)))
    if (!is.null(mets$conformal$q95)) {
      cat(sprintf("\n95%% Conformal absolute error (q95): ±%.3f (prediction ± q95)\n", as.numeric(mets$conformal$q95)))
    }
    if (!is.null(mets$feature_range)) {
      cat("\nApplicability: % of bands within training min–max is reported per sample in the table.\n")
    }
  })
  
  metrics_tbl <- reactive({
    req(fam(), input$props)
    f_id   <- fam()$id
    props  <- input$props
    if (!length(props)) return(NULL)
    
    # Build one row per property
    rows <- lapply(props, function(prop) {
      m <- read_metrics_json(f_id, prop)
      
      data.frame(
        Property        = get_label(prop),
        Family          = fam()$label,
        Grid_nm         = sprintf("%.0f–%.0f", min(fam()$wavegrid), max(fam()$wavegrid)),
        Features        = length(fam()$wavegrid),
        `Test RMSE`     = if (!is.null(m)) safe_num(m$metrics$RMSE) else NA,
        `Test R²`       = if (!is.null(m)) safe_num(m$metrics$R2)   else NA,
        `Test RPIQ`     = if (!is.null(m)) safe_num(m$metrics$RPIQ) else NA,
        `q95 (abs err)` = if (!is.null(m)) safe_num(m$conformal$q95) else NA,
        `Latent df`     = if (!is.null(m)) safe_num(m$latent$df)      else NA,
        `Latent χ²(0.95)` = if (!is.null(m)) safe_num(m$latent$thr95) else NA,
        check.names = FALSE
      )
    })
    
    # Order by your preferred property order
    df <- do.call(rbind, rows)
    ord <- match(df$Property, vapply(intersect(property_order, input$props), get_label, ""))
    df[order(ord), , drop = FALSE]
  })
  
  # --------------------- PREDICTIONS ---------------------
  preds <- eventReactive(input$run, {
    req(raw())
    df <- raw(); f <- fam(); wl <- get_wavelengths(df)
    vneed(length(wl$wl) > 20, "No spectral columns detected.")
    
    # resample to family grid
    X_src <- as.matrix(df[, wl$cols, drop = FALSE])
    X_res <- resample_to_grid(as.data.frame(X_src), src_wl = wl$wl, target_wl = f$wavegrid)
    rownames(X_res) <- df$Soil_ID
    
    # preprocessing (Absorbance -> SG(11,2,1)) unless disabled
    X_proc <- if (isTRUE(input$disable_pp)) X_res else apply_pipeline(X_res, f$preprocess)
    
    mdl_dir <- file.path("models", f$id, "models")
    props_wanted <- input$props
    vneed(length(props_wanted) > 0, "Select at least one property.")
    
    # check models exist
    paths_found <- vapply(props_wanted, function(p) {
      fp <- file.path(mdl_dir, paste0(p, ".h5"))
      if (file.exists(fp)) fp else NA_character_
    }, character(1))
    vneed(!all(is.na(paths_found)),
          paste0("No .h5 model files found for '", f$id, "'. Expected at: ", mdl_dir))
    
    out <- data.frame(Soil_ID = df$Soil_ID, stringsAsFactors = FALSE)
    
    for (prop in props_wanted) {
      fp_h5  <- file.path(mdl_dir, paste0(prop, ".h5"))
      if (!file.exists(fp_h5)) {
        out[[prop]] <- NA_real_
        out[[paste0(prop, "_MCsd")]] <- NA_real_
        out[[paste0(prop, "_PI95_low")]]  <- NA_real_
        out[[paste0(prop, "_PI95_high")]] <- NA_real_
        out[[paste0(prop, "_LatentApp_%")]] <- NA_real_
        next
      }
      
      if (!requireNamespace("keras", quietly = TRUE)) {
        vneed(FALSE, paste0("Keras model found for ", prop, " but 'keras' is not installed."))
      }
      mdl <- keras::load_model_hdf5(fp_h5, compile = FALSE)
      
      expected <- tryCatch(as.integer(mdl$inputs[[1]]$shape[[2]]), error = function(e) NA_integer_)
      found    <- ncol(X_proc)
      if (!is.na(expected) && expected != found) {
        vneed(FALSE, paste0(
          "Model input size mismatch for ", prop, ". Expects ", expected,
          " features; current family produced ", found, "."
        ))
      }
      
      # scaler
      base <- file.path(mdl_dir, prop)
      scaler_rds <- paste0(base, "_scaler.rds")
      if (file.exists(scaler_rds)) {
        sc <- readRDS(scaler_rds)
        mu_y <- if (is.null(sc$mean)) 0 else sc$mean
        sd_y <- if (is.null(sc$sd) || sc$sd == 0) 1 else sc$sd
      } else { mu_y <- 0; sd_y <- 1 }
      
      # metrics (for q95 and latent μ/Σ)
      mets <- read_metrics_json(f$id, prop)
      
      # ---------- MC Dropout inference (TF2 eager-friendly) ----------
      # ---------- MC Dropout inference (TF2 eager-friendly) ----------
      X_in <- as.matrix(X_proc)
      mc_T <- 30L
      
      # Build an encoder to grab the latent layer
      encoder <- keras::keras_model(
        inputs  = mdl$input,
        outputs = mdl$get_layer("latent")$output
      )
      
      mc_pred_mat <- matrix(NA_real_, nrow = nrow(X_in), ncol = mc_T)
      mc_lat_sum  <- NULL
      
      for (t in seq_len(mc_T)) {
        # Call the full model with dropout active
        out_t <- mdl(X_in, training = TRUE)           # list(rec, pred)
        
        # prediction is the SECOND output of your AE
        # Convert eager tensor -> R numeric
        yhat_t <- tryCatch(
          as.numeric(reticulate::py_to_r(out_t[[2]])),
          error = function(e) as.numeric(as.array(out_t[[2]]))  # fallback
        )
        mc_pred_mat[, t] <- yhat_t
        
        # Latent with dropout active
        Z_t <- tryCatch(
          reticulate::py_to_r(encoder(X_in, training = TRUE)),
          error = function(e) as.array(encoder(X_in, training = TRUE))
        )
        # ensure matrix shape
        Z_t <- as.matrix(Z_t)
        if (is.null(mc_lat_sum)) mc_lat_sum <- Z_t else mc_lat_sum <- mc_lat_sum + Z_t
      }
      
      yhat_scaled_mean <- rowMeans(mc_pred_mat)
      yhat_scaled_sd   <- apply(mc_pred_mat, 1, stats::sd)
      
      # back-transform to original units
      yhat <- yhat_scaled_mean * sd_y + mu_y
      mc_sd <- yhat_scaled_sd   * sd_y
      
      out[[prop]] <- yhat
      out[[paste0(prop, "_MCsd")]] <- mc_sd
      
      # ---------- Conformal PI95 (absolute residual q95 in metrics)
      if (!is.null(mets) && !is.null(mets$conformal$q95)) {
        q95 <- as.numeric(mets$conformal$q95)
        out[[paste0(prop, "_PI95_low")]]  <- yhat - q95
        out[[paste0(prop, "_PI95_high")]] <- yhat + q95
      } else {
        out[[paste0(prop, "_PI95_low")]]  <- NA_real_
        out[[paste0(prop, "_PI95_high")]] <- NA_real_
      }
      
      # ---------- Latent applicability (%)
      if (!is.null(mets) && !is.null(mets$latent$mu) && !is.null(mets$latent$Sigma)) {
        Z_mean <- mc_lat_sum / mc_T
        mu_z   <- as.numeric(mets$latent$mu)
        # rebuild Sigma matrix from JSON rows
        # rebuild Sigma matrix robustly
        Sigma_z <- coerce_sigma(mets$latent$Sigma, d_expected = length(mu_z))
        
        # ensure square & symmetric-ish
        if (nrow(Sigma_z) != length(mu_z) || ncol(Sigma_z) != length(mu_z)) {
          vneed(FALSE, paste0("Latent Σ dimension mismatch for ", prop, ": got ",
                              nrow(Sigma_z), "x", ncol(Sigma_z), ", expected ",
                              length(mu_z), "x", length(mu_z)))
        }
        # small symmetrization + ridge
        Sigma_z <- (Sigma_z + t(Sigma_z)) / 2
        Sinv <- tryCatch(inv_pd(Sigma_z, ridge = 1e-6),
                         error = function(e) {
                           if (requireNamespace("MASS", quietly = TRUE)) MASS::ginv(Sigma_z) else solve(Sigma_z)
                         })
        dlat <- length(mu_z)
        
        D2 <- apply(Z_mean, 1, function(z) {
          v <- z - mu_z
          as.numeric(t(v) %*% Sinv %*% v)
        })
        # map chi-square CDF to an intuitive 0..100 score
        latent_app <- pmax(0, 100 * (1 - stats::pchisq(D2, df = dlat)))
        out[[paste0(prop, "_LatentApp_%")]] <- latent_app
      } else {
        out[[paste0(prop, "_LatentApp_%")]] <- NA_real_
      }
    }
    
    out
  }, ignoreInit = TRUE)
  
  
  # Pretty table (2 decimals + fancy headers)
  preds_pretty <- reactive({
    req(preds())
    src <- preds()
    
    # round numerics
    for (cl in setdiff(names(src), "Soil_ID")) {
      src[[cl]] <- round(as.numeric(src[[cl]]), 2)
    }
    
    # rename columns (keys -> fancy labels)
    new_names <- names(src)
    for (k in soil_properties) {
      lbl <- get_label(k)
      new_names <- sub(paste0("^", k, "$"), lbl, new_names)
      new_names <- sub(paste0("^", k, "_PI95_low$"),     paste0(lbl, " (PI95 low)"),      new_names)
      new_names <- sub(paste0("^", k, "_PI95_high$"),    paste0(lbl, " (PI95 high)"),     new_names)
      new_names <- sub(paste0("^", k, "_MCsd$"),         paste0(lbl, " (MC sd)"),         new_names)
      new_names <- sub(paste0("^", k, "_LatentApp_%$"),  paste0(lbl, " (Latent app. %)"), new_names)
    }
    names(src) <- new_names
    src
  })
  output$model_info_dt <- renderDT({
    req(metrics_tbl())
    df <- metrics_tbl()
    
    # Which columns are numeric?
    num_cols <- c("Features","Test RMSE","Test R²","Test RPIQ","q95 (abs err)","Latent df","Latent χ²(0.95)")
    num_cols <- intersect(names(df), num_cols)
    
    DT::datatable(
      df,
      options = list(scrollX = TRUE, pageLength = 10),
      rownames = FALSE
    ) |>
      DT::formatRound(columns = setdiff(num_cols, c("Features","Latent df")), digits = 2) |>
      DT::formatRound(columns = intersect(num_cols, c("Features","Latent df")), digits = 0)
  })
  
  output$dl_metrics <- downloadHandler(
    filename = function() paste0("autoSpectra_model_metrics_", format(Sys.time(), "%Y%m%d_%H%M"), ".xlsx"),
    content  = function(file) {
      req(metrics_tbl())
      writexl::write_xlsx(metrics_tbl(), file)
    }
  )
  
  output$pred_dt <- renderDT({
    req(preds_pretty())
    df <- preds_pretty()
    num_cols <- setdiff(names(df), "Soil_ID")
    DT::datatable(df, options = list(scrollX = TRUE), rownames = FALSE) |>
      DT::formatRound(columns = num_cols, digits = 2)
  })
  
  output$dl_preds <- downloadHandler(
    filename = function() paste0("autoSpectra_predictions_", format(Sys.time(), "%Y%m%d_%H%M"), ".xlsx"),
    content  = function(file) { writexl::write_xlsx(preds_pretty(), file) }
  )
}

shinyApp(ui, server)
