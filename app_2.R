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
              img(src="logo.png", height="190px"),
              span(class="title", "autoSpectra")
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
                   downloadButton("dl_preds", "Download predictions")
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
  
  # --------------------- PREDICTIONS ---------------------
  preds <- eventReactive(input$run, {
    req(raw())
    df <- raw(); f <- fam(); wl <- get_wavelengths(df)
    vneed(length(wl$wl) > 20, "No spectral columns detected.")
    
    # Resample to family grid
    X_src <- as.matrix(df[, wl$cols, drop = FALSE])
    X_res <- resample_to_grid(as.data.frame(X_src), src_wl = wl$wl, target_wl = f$wavegrid)
    rownames(X_res) <- df$Soil_ID
    
    # Preprocessing (Absorbance -> SG(11,2,1)) unless disabled
    X_proc <- if (isTRUE(input$disable_pp)) X_res else apply_pipeline(X_res, f$preprocess)
    
    mdl_dir <- file.path("models", f$id, "models")
    props_wanted <- input$props
    vneed(length(props_wanted) > 0, "Select at least one property.")
    
    # Check files exist
    paths_found <- vapply(props_wanted, function(p) {
      fp <- file.path(mdl_dir, paste0(p, ".h5"))
      if (file.exists(fp)) fp else NA_character_
    }, character(1))
    vneed(!all(is.na(paths_found)),
          paste0("No .h5 model files found for '", f$id, "'. Expected at: ", mdl_dir))
    
    out <- data.frame(Soil_ID = df$Soil_ID, stringsAsFactors = FALSE)
    
    for (prop in props_wanted) {
      fp_h5  <- file.path(mdl_dir, paste0(prop, ".h5"))
      if (file.exists(fp_h5)) {
        if (!requireNamespace("keras", quietly = TRUE)) {
          vneed(FALSE, paste0("Keras model found for ", prop, " but 'keras' is not installed."))
        }
        mdl <- keras::load_model_hdf5(fp_h5, compile = FALSE)
        
        expected <- tryCatch(as.integer(mdl$inputs[[1]]$shape[[2]]), error = function(e) NA_integer_)
        found    <- ncol(X_proc)
        if (!is.na(expected) && expected != found) {
          vneed(FALSE, paste0(
            "Model input size mismatch for ", prop, ". Model expects ", expected,
            " features; current family produced ", found, "."
          ))
        }
        
        X_in <- as.matrix(X_proc)
        pred_out <- predict(mdl, X_in, verbose = 0)
        if (is.list(pred_out)) {
          if (!is.null(names(pred_out)) && "prediction" %in% names(pred_out)) {
            yhat_scaled <- pred_out[["prediction"]]
          } else {
            yhat_scaled <- pred_out[[length(pred_out)]]
          }
        } else {
          yhat_scaled <- pred_out
        }
        yhat_scaled <- as.numeric(yhat_scaled)
        
        base <- file.path(mdl_dir, prop)
        scaler_rds <- paste0(base, "_scaler.rds")
        if (file.exists(scaler_rds)) {
          sc <- readRDS(scaler_rds)
          mu <- if (is.null(sc$mean)) 0 else sc$mean
          sg <- if (is.null(sc$sd)   || sc$sd == 0) 1 else sc$sd
          yhat <- yhat_scaled * sg + mu
        } else {
          yhat <- yhat_scaled
        }
        out[[prop]] <- yhat
      } else {
        out[[prop]] <- NA_real_
      }
    }
    out
  }, ignoreInit = TRUE)
  
  # Pretty table (2 decimals + fancy headers)
  preds_pretty <- reactive({
    req(preds())
    src <- preds()
    num_cols <- setdiff(names(src), "Soil_ID")
    for (cl in num_cols) src[[cl]] <- round(as.numeric(src[[cl]]), 2)
    pretty_names <- c("Soil_ID", vapply(num_cols, get_label, FUN.VALUE = character(1)))
    names(src) <- pretty_names
    src
  })
  
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
