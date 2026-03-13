# app.R — autoSpectra (real-models edition)

suppressPackageStartupMessages({
  library(shiny)
  library(shinyWidgets)
  library(readxl)
  library(readr)
  library(DT)
  library(ggplot2)
  library(prospectr)  # SG
  library(writexl)
})

options(shiny.fullstacktrace = TRUE, shiny.maxRequestSize = 200*1024^2)

# Order the properties the way humans read them
property_order <- c(
  "soil_texture_sand","soil_texture_silt","soil_texture_clay",
  "organic_matter","soc","total_c","total_n","active_carbon",
  "ph","p","k","mg","fe","mn","zn","al","Ca","Cu","S","B",
  "pred_soil_protein","respiration","bd_ws"
)

# Fancy display labels (left = key, right = label)
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

# ensure labels exist for all keys (fallback to key if missing)
get_label <- function(key) if (!is.na(fancy_labels[[key]])) fancy_labels[[key]] else key

# ----------------------------- Safe validator ---------------------------------
vneed <- function(condition, message = "Validation failed.") {
  msg <- tryCatch({
    m <- paste(as.character(message), collapse = " ")
    if (!nzchar(m)) "Validation failed." else m
  }, error = function(e) "Validation failed.")
  shiny::validate(shiny::need(condition, msg))
}

# ----------------------------- Targets ----------------------------------------
soil_properties <- c(
  "soil_texture_sand","soil_texture_silt","soil_texture_clay",
  "organic_matter","soc","total_c","total_n",
  "active_carbon","ph","p","k",
  "mg","fe","mn","zn",
  "al","Ca","Cu","S",
  "B","pred_soil_protein","respiration","bd_ws"
)

# ----------------------------- Registry ---------------------------------------
default_registry <- list(
  ASD_DRY = list(
    id = "ASD_DRY",
    label = "ASD — DRY (23 props)",
    sensors_allowed = c("ASD"),
    moisture_levels = c("DRY"),
    properties = soil_properties,
    wavegrid = 350:2500,
    preprocess = c("ABSORBANCE","SG(11,2,1)")
  ),
  NeoSpectra_DRY = list(
    id = "NeoSpectra_DRY",
    label = "NeoSpectra — DRY (23 props)",
    sensors_allowed = c("NeoSpectra"),
    moisture_levels = c("DRY"),
    properties = soil_properties,
    wavegrid = 1350:2500,
    preprocess = c("ABSORBANCE","SG(11,2,1)")
  ),
  NaturaSpec_DRY = list(
    id = "NaturaSpec_DRY",
    label = "NaturaSpec — DRY (23 props)",
    sensors_allowed = c("NaturaSpec"),
    moisture_levels = c("DRY"),
    properties = soil_properties,
    wavegrid = 350:2500,
    preprocess = c("ABSORBANCE","SG(11,2,1)")
  ),
  Agnostic_DRY = list(
    id = "Agnostic_DRY",
    label = "Agnostic — DRY (ASD + NaturaSpec + NeoSpectra, 23 props)",
    sensors_allowed = c("ASD","NaturaSpec","NeoSpectra"),
    moisture_levels = c("agnostic"),
    properties = soil_properties,
    wavegrid = 1350:2500,  # intersection grid
    preprocess = c("ABSORBANCE","SG(11,2,1)")
  ),
  Agnostic_Moisture = list(
    id = "Agnostic_Moisture",
    label = "Agnostic — DRY+1ML+3ML (ASD + NaturaSpec + NeoSpectra, 23 props)",
    sensors_allowed = c("ASD","NaturaSpec","NeoSpectra"),
    moisture_levels = c("agnostic"),
    properties = soil_properties,
    wavegrid = 1350:2500,  # intersection grid
    preprocess = c("ABSORBANCE","SG(11,2,1)")
  )
)

family_matches <- function(fam, sensor, moisture) {
  ok_sensor <- sensor %in% fam$sensors_allowed
  ok_moist  <- moisture %in% fam$moisture_levels || ("agnostic" %in% fam$moisture_levels)
  ok_sensor && ok_moist
}

# ------------------------ Wavelength & preprocessing --------------------------
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

# ----- Reflectance -> Absorbance (A = -ln R; auto-scale % to [0,1]) -----------
reflectance_to_absorbance <- function(M, base10 = FALSE) {
  Mnum <- as.matrix(M)
  if (is.finite(max(Mnum, na.rm = TRUE)) && max(Mnum, na.rm = TRUE) > 2) {
    Mnum <- Mnum / 100
  }
  eps <- 1e-6
  Mnum[Mnum < eps] <- eps
  if (base10) -log10(Mnum) else -log(Mnum)
}

# ----- SNV (not used now, but kept available) ---------------------------------
apply_snv <- function(M) {
  M_centered <- sweep(M, 1, rowMeans(M, na.rm = TRUE), "-")
  denom <- sqrt(rowMeans(M_centered^2, na.rm = TRUE))
  denom[denom == 0 | is.na(denom)] <- 1
  sweep(M_centered, 1, denom, "/")
}

# ----- SG (robust; supports derivative d via "SG(m,p,d)") ---------------------
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

# --------------------------------- UI -----------------------------------------
ui <- fluidPage(
  tags$head(
    tags$link(rel="icon", type="image/png", href="logo.png")
  ),
  
  titlePanel(
    div(
      style = "display: flex; align-items: center;",
      img(src = "logo.png", height = "160px", style = "margin-right: 20px;"),
      span("autoSpectra — Soil Spectral Prediction", style = "font-size: 24px; font-weight: bold;")
    )
  ),
  sidebarLayout(
    sidebarPanel(
      pickerInput("sensor", "Sensor",
                  choices = c("ASD","NeoSpectra","NaturaSpec"), selected = "ASD"),
      pickerInput("moisture", "Moisture mode",
                  choices = c("DRY","1ML","3ML","agnostic"), selected = "DRY"),
      pickerInput("family", "Model family",
                  choices = setNames(names(default_registry),
                                     vapply(default_registry, `[[`, "", "label")),
                  selected = "ASD_DRY"),
      checkboxInput("disable_pp", "Disable preprocessing (debug)", value = FALSE),
      tags$hr(),
      fileInput("file", "Upload Excel (.xlsx/.xls) or CSV", accept = c(".xlsx",".xls",".csv")),
      uiOutput("sheet_ui"),
      textInput("soil_col", "Soil ID column", value = "Soil_ID"),
      actionButton("preview", "Preview", class = "btn btn-primary"),
      tags$hr(),
      uiOutput("props_ui"),
      actionButton("run", "Predict", class = "btn btn-success"),
      tags$hr(),
      downloadButton("dl_preds", "Download predictions")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Preview",
                 verbatimTextOutput("info"),
                 DTOutput("head_dt")),
        tabPanel("Spectrum",
                 uiOutput("soil_pick_ui"),
                 plotOutput("spec_plot", height = "320px"),
                 verbatimTextOutput("coverage")),
        tabPanel("Predictions",
                 DTOutput("pred_dt"))
      )
    )
  )
)

# -------------------------------- Server --------------------------------------
server <- function(input, output, session) {
  
  # filter family picker by sensor + moisture
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
  
  output$props_ui <- renderUI({
    f <- fam(); req(f)
    # order choices by our preferred order
    ordered <- intersect(property_order, f$properties)
    # map keys -> fancy labels (display), values remain keys (important!)
    choices <- setNames(ordered, vapply(ordered, get_label, FUN.VALUE = character(1)))
    pickerInput("props", "Properties", choices = choices, selected = ordered, multiple = TRUE)
  })
  
  # sheet selector for Excel
  output$sheet_ui <- renderUI({
    req(input$file)
    ext <- tolower(tools::file_ext(input$file$name))
    if (ext %in% c("xlsx","xls")) {
      path <- normalizePath(as.character(input$file$datapath), winslash = "/")
      shs  <- as.character(readxl::excel_sheets(path))
      selectInput("sheet", "Sheet", choices = shs, selected = shs[1])
    }
  })
  
  # read data on Preview
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
    soil_col <- as.character(input$soil_col)
    vneed(soil_col %in% names(df), paste0("Column '", soil_col, "' not found."))
    names(df)[names(df) == soil_col] <- "Soil_ID"
    df
  }, ignoreInit = TRUE)
  
  # preview info
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
    cat("Preprocess:", paste(f$preprocess, collapse = " -> "),
        if (isTRUE(input$disable_pp)) " [DISABLED]" else "", "\n")
  })
  
  output$head_dt <- renderDT({
    req(raw())
    DT::datatable(head(raw(), 5), options = list(scrollX = TRUE))
  })
  
  # spectrum viewer
  output$soil_pick_ui <- renderUI({
    req(raw())
    selectInput("soil_pick", "Sample", choices = raw()$Soil_ID)
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
  
  # predictions
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
    if (length(props_wanted) == 0) vneed(FALSE, "Select at least one property.")
    
    # sanity: do we have any files?
    paths_found <- vapply(props_wanted, function(p) {
      fp <- file.path(mdl_dir, paste0(p, ".h5"))
      if (file.exists(fp)) fp else NA_character_
    }, character(1))
    if (all(is.na(paths_found))) {
      vneed(FALSE, paste0(
        "No .h5 model files found for family '", f$id, "'. Expected at: ", mdl_dir,
        ". Place your trained .h5 files there (one per property)."
      ))
    }
    
    out <- data.frame(Soil_ID = df$Soil_ID, stringsAsFactors = FALSE)
    
    # ---- Post-process for display: round to 2 decimals & fancy column names
    disp <- out
    num_cols <- setdiff(names(disp), "Soil_ID")
    for (cl in num_cols) {
      # round safely and keep as numeric (so DT can still sort)
      disp[[cl]] <- round(as.numeric(disp[[cl]]), 2)
    }
    # rename columns to fancy labels
    nice_names <- c("Soil_ID", vapply(num_cols, get_label, FUN.VALUE = character(1)))
    names(disp) <- nice_names
    disp
    
    for (prop in props_wanted) {
      fp_h5  <- file.path(mdl_dir, paste0(prop, ".h5"))
      if (file.exists(fp_h5)) {
        if (!requireNamespace("keras", quietly = TRUE)) {
          vneed(FALSE, paste0("Keras model found for ", prop, " but 'keras' is not installed."))
        }
        mdl <- keras::load_model_hdf5(fp_h5, compile = FALSE)
        
        # shape check
        expected <- tryCatch(as.integer(mdl$inputs[[1]]$shape[[2]]), error = function(e) NA_integer_)
        found    <- ncol(X_proc)
        if (!is.na(expected) && expected != found) {
          vneed(FALSE, paste0(
            "Model input size mismatch for ", prop, ". Model expects ", expected,
            " features, current family produced ", found, ".\n",
            "Tip: pick the family whose wavegrid matches the model's training grid."
          ))
        }
        
        # run predict
        X_in <- as.matrix(X_proc)
        pred_out <- predict(mdl, X_in, verbose = 0)
        
        # pick the prediction head (second output), or by name if available
        if (is.list(pred_out)) {
          if (!is.null(names(pred_out)) && "prediction" %in% names(pred_out)) {
            yhat_scaled <- pred_out[["prediction"]]
          } else {
            # assume last element is the prediction head
            yhat_scaled <- pred_out[[length(pred_out)]]
          }
        } else {
          # single-output fallback
          yhat_scaled <- pred_out
        }
        yhat_scaled <- as.numeric(yhat_scaled)
        
        # back-transform with scaler if present
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
  
  # Pretty table for UI & download: round to 2 decimals + fancy headers
  preds_pretty <- reactive({
    req(preds())
    src <- preds()
    
    # Round to 2 decimals (keep numeric)
    num_cols <- setdiff(names(src), "Soil_ID")
    for (cl in num_cols) src[[cl]] <- round(as.numeric(src[[cl]]), 2)
    
    # Rename columns to fancy labels
    fancy_names <- c("Soil_ID", vapply(num_cols, get_label, FUN.VALUE = character(1)))
    names(src) <- fancy_names
    
    src
  })
  
  output$pred_dt <- renderDT({
    req(preds_pretty())
    df <- preds_pretty()
    # all numeric (non-ID) columns:
    num_cols <- setdiff(names(df), "Soil_ID")
    DT::datatable(df, options = list(scrollX = TRUE), rownames = FALSE) |>
      DT::formatRound(columns = num_cols, digits = 2)
  })
  
  output$dl_preds <- downloadHandler(
    filename = function() paste0("autoSpectra_predictions_", format(Sys.time(), "%Y%m%d_%H%M"), ".xlsx"),
    content  = function(file) {
      writexl::write_xlsx(preds_pretty(), file)
    }
  )
}

shinyApp(ui, server)
