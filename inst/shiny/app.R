# inst/shiny/app.R — autoSpectra Shiny interface
# Supports: VisNIR + MIR spectra | local + OSSL model families
# Preprocessing: ABSORBANCE -> SG_SMOOTH -> SG_DERIV (two-step)

suppressPackageStartupMessages({
  library(shiny)
  library(shinyWidgets)
  library(readxl)
  library(readr)
  library(DT)
  library(ggplot2)
  library(prospectr)
  library(writexl)
})

# Load autoSpectra if installed; otherwise source from parent directory
if (requireNamespace("autoSpectra", quietly = TRUE)) {
  library(autoSpectra)
} else {
  # Development mode: source R files from package source tree
  pkg_root <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", "..")),
    error = function(e) getwd()
  )
  for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE))
    source(f)
  model_dir_default <- file.path(pkg_root, "models")
}

# Resolve model directory (can be set via option)
model_dir <- getOption("autoSpectra.model_dir",
                        default = if (exists("model_dir_default")) model_dir_default else "models")

# ---- Registry snapshot for UI ----------------------------------------
reg         <- model_registry
reg_ids     <- names(reg)
reg_labels  <- vapply(reg, `[[`, "", "label")

# Sensor type display
sensor_type_label <- function(fam) {
  if (identical(fam$sensor_type, "mir")) "MIR (cm\u207b\u00b9)" else "VisNIR (nm)"
}

# ---- UI --------------------------------------------------------------
ui <- fluidPage(
  tags$head(
    tags$link(rel = "icon", type = "image/png", href = "logo.png"),
    tags$style(HTML("
      .sidebar { padding-top: 10px; }
      .section-header { font-weight: bold; color: #2c5f8a; margin-top: 10px; }
      .ossl-tag { background: #e8f4fd; border: 1px solid #90caf9;
                  border-radius: 4px; padding: 2px 6px; font-size: 0.8em; }
    "))
  ),

  titlePanel(
    div(
      style = "display:flex; align-items:center; gap:16px;",
      img(src = "logo.png", height = "80px"),
      div(
        span("autoSpectra", style = "font-size:26px; font-weight:bold;"),
        br(),
        span("Soil Spectral Modeling \u2014 Prediction \u2014 Visualization",
             style = "font-size:13px; color:#555;")
      )
    )
  ),

  sidebarLayout(
    sidebarPanel(width = 3,

      # --- Sensor & moisture ---
      div(class = "section-header", "1. Sensor / Mode"),
      pickerInput("sensor_type", "Spectrum type",
                  choices = c("VisNIR (nm)" = "visnir", "MIR (cm\u207b\u00b9)" = "mir"),
                  selected = "visnir"),
      pickerInput("sensor", "Sensor (for local models)",
                  choices = c("ASD", "NeoSpectra", "NaturaSpec", "Other"),
                  selected = "ASD"),
      pickerInput("moisture", "Moisture",
                  choices = c("DRY", "1ML", "3ML", "agnostic"),
                  selected = "DRY"),

      tags$hr(),

      # --- Model family ---
      div(class = "section-header", "2. Model family"),
      uiOutput("family_ui"),

      tags$hr(),

      # --- Upload ---
      div(class = "section-header", "3. Upload spectra"),
      fileInput("file", NULL,
                accept = c(".xlsx", ".xls", ".csv"),
                placeholder = "Excel / CSV"),
      uiOutput("sheet_ui"),
      textInput("soil_col", "Sample ID column", value = "Soil_ID"),
      actionButton("preview_btn", "Preview", class = "btn-primary btn-sm"),

      tags$hr(),

      # --- Properties ---
      div(class = "section-header", "4. Properties to predict"),
      uiOutput("props_ui"),

      tags$hr(),

      # --- Options ---
      checkboxInput("disable_pp", "Disable preprocessing (debug)", FALSE),

      tags$hr(),
      actionButton("run_btn", "\u25b6 Predict", class = "btn-success"),
      tags$hr(),
      downloadButton("dl_preds", "Download predictions (.xlsx)")
    ),

    mainPanel(width = 9,
      tabsetPanel(id = "tabs",

        tabPanel("Preview",
          verbatimTextOutput("info_txt"),
          DTOutput("head_dt")
        ),

        tabPanel("Spectrum viewer",
          uiOutput("sample_pick_ui"),
          plotOutput("spec_plot", height = "320px"),
          verbatimTextOutput("coverage_txt")
        ),

        tabPanel("Mean spectrum",
          plotOutput("mean_spec_plot", height = "360px")
        ),

        tabPanel("Predictions",
          DTOutput("pred_dt")
        ),

        tabPanel("About / Help",
          includeMarkdown(
            system.file("shiny/HELP.md", package = "autoSpectra",
                        mustWork = FALSE) %||%
            textConnection("## autoSpectra\nLoad spectra, select a model family, and press Predict.")
          )
        )
      )
    )
  )
)

# ---- `%||%` null-coalescing operator ----------------------------------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && nzchar(a)) a else b

# ---- Server -----------------------------------------------------------
server <- function(input, output, session) {

  # --- Reactive: filtered family list ---
  avail_families <- reactive({
    req(input$sensor_type)
    ids <- reg_ids[vapply(reg, function(f) {
      type_ok <- identical(f$sensor_type, input$sensor_type)
      sens_ok <- is.null(f$sensors_allowed) ||
        input$sensor %in% f$sensors_allowed
      moi_ok  <- is.null(input$moisture) ||
        input$moisture %in% f$moisture_levels ||
        "agnostic" %in% f$moisture_levels
      type_ok  # always show all matching type; filter narrows via sensor + moi
    }, logical(1))]
    if (length(ids) == 0) ids <- reg_ids
    setNames(ids, vapply(reg[ids], `[[`, "", "label"))
  })

  output$family_ui <- renderUI({
    choices <- avail_families()
    pickerInput("family", NULL,
                choices  = choices,
                selected = choices[1],
                options  = pickerOptions(liveSearch = TRUE))
  })

  fam <- reactive({
    req(input$family)
    model_registry[[as.character(input$family)]]
  })

  # --- Properties UI ---
  output$props_ui <- renderUI({
    f <- fam(); req(f)
    props <- f$properties
    labels <- vapply(props, property_label, character(1))
    pickerInput("props", NULL,
                choices  = setNames(props, labels),
                selected = props,
                multiple = TRUE,
                options  = pickerOptions(
                  actionsBox = TRUE,
                  liveSearch = TRUE,
                  selectedTextFormat = "count > 3"
                ))
  })

  # --- Sheet selector for Excel ---
  output$sheet_ui <- renderUI({
    req(input$file)
    ext <- tolower(tools::file_ext(input$file$name))
    if (ext %in% c("xlsx", "xls")) {
      path <- normalizePath(input$file$datapath, winslash = "/")
      shs  <- readxl::excel_sheets(path)
      selectInput("sheet", "Sheet", choices = shs, selected = shs[1])
    }
  })

  # --- Read data on Preview ---
  raw <- eventReactive(input$preview_btn, {
    req(input$file)
    ext  <- tolower(tools::file_ext(input$file$name))
    path <- normalizePath(input$file$datapath, winslash = "/")
    df   <- if (ext %in% c("xlsx", "xls")) {
      sh <- if (isTruthy(input$sheet)) input$sheet else 1
      readxl::read_excel(path, sheet = sh)
    } else {
      readr::read_csv(path, show_col_types = FALSE)
    }
    df       <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
    soil_col <- as.character(input$soil_col)
    vneed(soil_col %in% names(df),
          paste0("Column '", soil_col, "' not found. ",
                 "Available: ", paste(names(df)[1:min(8, ncol(df))], collapse = ", ")))
    names(df)[names(df) == soil_col] <- "Soil_ID"
    df
  }, ignoreInit = TRUE)

  # --- Preview info ---
  output$info_txt <- renderPrint({
    req(raw()); df <- raw(); f <- fam()
    wl <- get_wavelengths(df)
    cat(sprintf("Samples:  %d\n", nrow(df)))
    if (length(wl$wl) > 0)
      cat(sprintf("Spectral: %d bands  (%.0f \u2013 %.0f)\n",
                  length(wl$wl), min(wl$wl), max(wl$wl)))
    cat(sprintf("Family:   %s\n", f$label))
    cat(sprintf("Type:     %s\n", sensor_type_label(f)))
    cat(sprintf("Grid:     %d bands  (%.0f \u2013 %.0f)\n",
                length(f$wavegrid), min(f$wavegrid), max(f$wavegrid)))
    cat(sprintf("Pipeline: %s%s\n",
                paste(f$preprocess, collapse = " \u2192 "),
                if (isTRUE(input$disable_pp)) "  [DISABLED]" else ""))
  })

  output$head_dt <- renderDT({
    req(raw())
    DT::datatable(head(raw(), 5), options = list(scrollX = TRUE))
  })

  # --- Spectrum viewer ---
  output$sample_pick_ui <- renderUI({
    req(raw())
    selectInput("sample_pick", "Sample", choices = raw()$Soil_ID)
  })

  output$spec_plot <- renderPlot({
    req(raw(), input$sample_pick)
    df <- raw(); f <- fam()
    row <- df[df$Soil_ID == input$sample_pick, , drop = FALSE]
    wl  <- get_wavelengths(row)
    vals <- as.numeric(row[1, wl$cols, drop = TRUE])
    dd   <- data.frame(wl = wl$wl, y = vals)

    x_lab <- if (identical(f$sensor_type, "mir"))
      "Wavenumber (cm\u207b\u00b9)" else "Wavelength (nm)"

    ggplot(dd, aes(wl, y)) +
      geom_line(colour = "#2166ac", linewidth = 0.6) +
      geom_vline(xintercept = range(f$wavegrid),
                 linetype = "dashed", colour = "grey50") +
      labs(x = x_lab, y = "Response",
           subtitle = sprintf("Model grid: %.0f \u2013 %.0f  (%d bands)",
                              min(f$wavegrid), max(f$wavegrid),
                              length(f$wavegrid))) +
      theme_minimal(base_size = 12)
  })

  output$coverage_txt <- renderText({
    req(raw()); df <- raw(); f <- fam()
    wl <- get_wavelengths(df)
    ov <- sum(round(wl$wl) %in% round(f$wavegrid))
    sprintf("Overlap with model grid: %.1f%%  (%d / %d bands). Missing bands interpolated.",
            100 * ov / length(f$wavegrid), ov, length(f$wavegrid))
  })

  # --- Mean spectrum ---
  output$mean_spec_plot <- renderPlot({
    req(raw()); df <- raw(); f <- fam()
    wl_info <- get_wavelengths(df)
    M <- as.matrix(df[, wl_info$cols, drop = FALSE])
    x_lab <- if (identical(f$sensor_type, "mir"))
      "Wavenumber (cm\u207b\u00b9)" else "Wavelength (nm)"
    plot_mean_spectrum(M, wl = wl_info$wl, xlab = x_lab,
                       title = paste("Mean \u00b1 SD \u2014", nrow(df), "samples"))
  })

  # --- Predictions ---
  preds <- eventReactive(input$run_btn, {
    req(raw())
    df    <- raw()
    f     <- fam()
    props <- as.character(input$props)
    vneed(length(props) > 0, "Select at least one property.")

    wl_info <- get_wavelengths(df)
    vneed(length(wl_info$wl) > 20, "No spectral columns detected (need > 20 numeric columns).")

    # Resample
    X_src <- as.matrix(df[, wl_info$cols, drop = FALSE])
    X_res <- resample_to_grid(X_src, wl_info$wl, f$wavegrid)
    rownames(X_res) <- df$Soil_ID

    # Preprocess
    X_proc <- if (isTRUE(input$disable_pp)) X_res else apply_pipeline(X_res, f$preprocess)

    mdl_dir <- file.path(model_dir, f$id, "models")
    any_found <- any(vapply(props, function(p)
      file.exists(file.path(mdl_dir, paste0(p, ".h5"))), logical(1)))
    vneed(any_found, paste0(
      "No model files found in: ", mdl_dir, "\n",
      "Train models with train_ossl_models('", f$id, "') or train_soilVAE()."))

    out <- data.frame(Soil_ID = df$Soil_ID, stringsAsFactors = FALSE)

    for (prop in props) {
      fp <- file.path(mdl_dir, paste0(prop, ".h5"))
      if (!file.exists(fp)) { out[[prop]] <- NA_real_; next }
      if (!requireNamespace("keras", quietly = TRUE)) {
        vneed(FALSE, "Package 'keras' is required for predictions.")
      }
      info <- load_soilVAE(f$id, prop, model_dir)
      mdl  <- info$model; sc <- info$scaler

      exp_d <- tryCatch(as.integer(mdl$inputs[[1]]$shape[[2]]),
                        error = function(e) NA_integer_)
      if (!is.na(exp_d) && exp_d != ncol(X_proc))
        vneed(FALSE, paste0("Input size mismatch for [", prop, "]: model=", exp_d,
                            " vs. family=", ncol(X_proc)))

      yhat_z    <- .extract_prediction(mdl, X_proc)
      mu        <- if (is.null(sc$mean)) 0 else sc$mean
      sg        <- if (is.null(sc$sd) || sc$sd == 0) 1 else sc$sd
      out[[prop]] <- yhat_z * sg + mu
    }
    out
  }, ignoreInit = TRUE)

  preds_pretty <- reactive({
    req(preds())
    format_predictions(preds())
  })

  output$pred_dt <- renderDT({
    req(preds_pretty())
    df <- preds_pretty()
    num_cols <- setdiff(names(df), "Soil_ID")
    DT::datatable(df, rownames = FALSE,
                  options = list(scrollX = TRUE, pageLength = 20)) |>
      DT::formatRound(columns = num_cols, digits = 2)
  })

  output$dl_preds <- downloadHandler(
    filename = function()
      paste0("autoSpectra_", input$family, "_",
             format(Sys.time(), "%Y%m%d_%H%M"), ".xlsx"),
    content = function(file) writexl::write_xlsx(preds_pretty(), file)
  )
}

shinyApp(ui, server)
