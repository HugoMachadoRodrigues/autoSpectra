# inst/shiny/app.R - autoSpectra Shiny interface
#
# Four OSSL models available:
#   OSSL_VisNIR        - masked soilVAE (sensor-agnostic, 350-2500 nm)
#   OSSL_MIR           - masked soilVAE (sensor-agnostic, 600-4000 cm-1)
#   OSSL_VisNIR_Cubist - Cubist + PCA(120) + SNV (Safanelli 2025)
#   OSSL_MIR_Cubist    - Cubist + PCA(120) + SNV (Safanelli 2025)

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

if (requireNamespace("autoSpectra", quietly = TRUE)) {
  library(autoSpectra)
} else {
  pkg_root <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", "..")),
    error = function(e) getwd()
  )
  for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                        full.names = TRUE))
    source(f)
}

model_dir <- getOption("autoSpectra.model_dir", default = "models")

# ---- Family choices ------------------------------------------------------
FAMILIES <- list(
  "VisNIR - masked soilVAE (sensor-agnostic)"      = "OSSL_VisNIR",
  "MIR - masked soilVAE (sensor-agnostic)"         = "OSSL_MIR",
  "VisNIR - Cubist + PCA(120) [Safanelli 2025]"    = "OSSL_VisNIR_Cubist",
  "MIR - Cubist + PCA(120) [Safanelli 2025]"       = "OSSL_MIR_Cubist"
)
FAMILY_PROPS <- lapply(FAMILIES, function(fid) get_family(fid)$properties)

# ---- UI ------------------------------------------------------------------
ui <- fluidPage(
  tags$head(
    tags$link(rel = "icon", type = "image/png", href = "logo.png"),
    tags$style(HTML("
      .sidebar        { padding-top: 10px; }
      .section-header { font-weight: bold; color: #2c5f8a; margin-top: 10px; }
      .model-tag      { background: #e8f5e9; border: 1px solid #81c784;
                        border-radius: 4px; padding: 2px 8px;
                        font-size: 0.82em; }
      .ossl-badge     { background: #e3f2fd; border: 1px solid #90caf9;
                        border-radius: 4px; padding: 2px 8px;
                        font-size: 0.82em; }
      .algo-vae       { background: #fff3e0; border: 1px solid #ffb74d;
                        border-radius: 4px; padding: 2px 8px;
                        font-size: 0.82em; }
      .algo-cubist    { background: #f3e5f5; border: 1px solid #ba68c8;
                        border-radius: 4px; padding: 2px 8px;
                        font-size: 0.82em; }
      .warn-ood       { color: #b71c1c; font-weight: bold; }
    "))
  ),

  titlePanel(
    div(
      style = "display:flex; align-items:center; gap:16px;",
      img(src = "logo.png", height = "80px"),
      div(
        span("autoSpectra", style = "font-size:26px; font-weight:bold;"),
        br(),
        span("Soil Spectral Prediction \u2014 OSSL v1.2 Sensor-Agnostic Models",
             style = "font-size:13px; color:#555;")
      )
    )
  ),

  sidebarLayout(
    sidebarPanel(width = 3,

      div(class = "section-header", "1. Spectral Model"),
      pickerInput(
        "family_id", NULL,
        choices  = FAMILIES,
        selected = "OSSL_VisNIR",
        options  = pickerOptions(style = "btn-outline-success")
      ),
      uiOutput("family_info_ui"),

      tags$hr(),

      div(class = "section-header", "2. Soil Properties"),
      uiOutput("props_ui"),
      actionButton("select_all_btn",  "All",  class = "btn-xs btn-outline-secondary"),
      actionButton("select_none_btn", "None", class = "btn-xs btn-outline-secondary"),

      tags$hr(),

      div(class = "section-header", "3. Upload Spectra"),
      fileInput("file", NULL,
                accept = c(".xlsx", ".xls", ".csv"),
                placeholder = "Excel / CSV"),
      uiOutput("sheet_ui"),
      textInput("soil_col", "Sample ID column", value = "Soil_ID"),
      actionButton("preview_btn", "\U0001f441 Preview",
                   class = "btn-primary btn-sm w-100"),

      tags$hr(),

      div(class = "section-header", "4. Predict"),
      checkboxInput("check_domain", "Show applicability domain", value = FALSE),
      checkboxInput("show_pi",      "Prediction intervals (Cubist only)", value = FALSE),
      actionButton("predict_btn", "\U0001f9ea Predict",
                   class = "btn-success btn-sm w-100"),

      tags$hr(),

      div(class = "section-header", "5. Export"),
      downloadButton("dl_excel", "Download Excel",
                     class = "btn-outline-primary btn-sm w-100")
    ),

    mainPanel(width = 9,
      tabsetPanel(id = "tabs",
        tabPanel("Preview",
          br(),
          uiOutput("preview_summary"),
          br(),
          tableOutput("preview_head")
        ),
        tabPanel("Spectrum Viewer",
          br(),
          uiOutput("sample_picker_ui"),
          plotOutput("spec_plot", height = "380px")
        ),
        tabPanel("Mean Spectrum",
          br(),
          plotOutput("mean_spec_plot", height = "420px")
        ),
        tabPanel("Predictions",
          br(),
          uiOutput("pred_status"),
          DT::dataTableOutput("pred_table"),
          br(),
          uiOutput("domain_ui")
        )
      )
    )
  )
)

# ---- Server --------------------------------------------------------------
server <- function(input, output, session) {

  fam <- reactive({ get_family(input$family_id) })

  output$family_info_ui <- renderUI({
    f <- fam()
    n_bands <- length(f$wavegrid)
    rng <- paste0(min(f$wavegrid), "\u2013", max(f$wavegrid),
                  if (f$sensor_type == "mir") " cm\u207b\u00b9" else " nm")
    algo_class <- if (identical(f$algorithm, "cubist_pca")) "algo-cubist"
                   else "algo-vae"
    algo_label <- if (identical(f$algorithm, "cubist_pca"))
                    "Cubist + PCA(120)" else "masked soilVAE"
    tagList(
      tags$small(
        span(class = "ossl-badge", "OSSL v1.2"), " ",
        span(class = algo_class, algo_label), br(),
        n_bands, " bands | ", rng
      )
    )
  })

  output$props_ui <- renderUI({
    props  <- fam()$properties
    labels <- setNames(ossl_l1_labels[props], props)
    labels[is.na(labels)] <- props[is.na(labels)]
    pickerInput(
      "properties", NULL,
      choices  = setNames(props, labels),
      selected = props,
      multiple = TRUE,
      options  = pickerOptions(
        actionsBox         = FALSE,
        liveSearch         = TRUE,
        selectedTextFormat = "count > 3",
        countSelectedText  = "{0} properties selected"
      )
    )
  })

  observeEvent(input$select_all_btn, {
    updatePickerInput(session, "properties", selected = fam()$properties)
  })
  observeEvent(input$select_none_btn, {
    updatePickerInput(session, "properties", selected = character(0))
  })

  raw_df <- reactiveVal(NULL)

  sheets_available <- reactive({
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    if (tolower(ext) %in% c("xlsx", "xls"))
      readxl::excel_sheets(input$file$datapath)
    else character(0)
  })

  output$sheet_ui <- renderUI({
    sh <- sheets_available()
    if (length(sh) > 0)
      pickerInput("sheet", "Sheet", choices = sh, selected = sh[1])
  })

  read_file <- reactive({
    req(input$file)
    ext <- tolower(tools::file_ext(input$file$name))
    if (ext %in% c("xlsx", "xls")) {
      sh <- if (!is.null(input$sheet)) input$sheet else 1
      readxl::read_excel(input$file$datapath, sheet = sh)
    } else {
      readr::read_csv(input$file$datapath, show_col_types = FALSE)
    }
  })

  observeEvent(input$preview_btn, {
    df <- tryCatch(as.data.frame(read_file()), error = function(e) NULL)
    if (is.null(df)) {
      showNotification("Could not read file.", type = "error")
      return()
    }
    id_col <- input$soil_col
    if (!id_col %in% names(df)) {
      candidates <- grep("id|sample|soil", names(df), ignore.case = TRUE, value = TRUE)
      if (length(candidates) > 0) {
        id_col <- candidates[1]
        updateTextInput(session, "soil_col", value = id_col)
        showNotification(paste("Auto-detected ID column:", id_col), type = "message")
      } else {
        showNotification(paste("Column", input$soil_col, "not found."), type = "warning")
        return()
      }
    }
    names(df)[names(df) == id_col] <- "Soil_ID"
    raw_df(df)
    updateTabsetPanel(session, "tabs", selected = "Preview")
    showNotification("File loaded.", type = "message", duration = 2)
  })

  output$preview_summary <- renderUI({
    df <- raw_df(); req(df)
    wl_info <- get_wavelengths(df, id_col = "Soil_ID")
    f  <- fam()
    overlap <- sum(wl_info$wl >= min(f$wavegrid) & wl_info$wl <= max(f$wavegrid))
    tagList(
      tags$p(
        strong("Samples: "), nrow(df), " | ",
        strong("Spectral bands: "), length(wl_info$wl), " | ",
        strong("Overlap with model grid: "), overlap, " bands"
      ),
      if (overlap < 100)
        tags$p(class = "warn-ood",
               "\u26a0 Low spectral overlap with the selected model. ",
               "Check that the file's wavelength range matches the selected domain.")
    )
  })

  output$preview_head <- renderTable({
    df <- raw_df(); req(df)
    wl_info   <- get_wavelengths(df, id_col = "Soil_ID")
    meta_cols <- setdiff(names(df), wl_info$cols)
    spec_cols <- head(wl_info$cols, 5)
    head(df[, c(meta_cols[1:min(3, length(meta_cols))], spec_cols)], 6)
  }, striped = TRUE, hover = TRUE, digits = 4)

  output$sample_picker_ui <- renderUI({
    df <- raw_df(); req(df)
    pickerInput("sample_id", "Sample", choices = df[["Soil_ID"]],
                selected = df[["Soil_ID"]][1],
                options  = pickerOptions(liveSearch = TRUE))
  })

  output$spec_plot <- renderPlot({
    df <- raw_df(); req(df, input$sample_id)
    f  <- fam()
    sub_df <- df[df$Soil_ID == input$sample_id, , drop = FALSE]
    plot_spectra(sub_df, family = f,
                 title = paste("Spectrum:", input$sample_id))
  })

  output$mean_spec_plot <- renderPlot({
    df <- raw_df(); req(df)
    plot_mean_spectrum(df, family = fam())
  })

  pred_df    <- reactiveVal(NULL)
  domain_df  <- reactiveVal(NULL)

  observeEvent(input$predict_btn, {
    df <- raw_df()
    vneed(!is.null(df), "Please upload and preview a file first.")
    props <- input$properties
    vneed(length(props) > 0, "Please select at least one property.")

    f         <- fam()
    is_cubist <- identical(f$algorithm, "cubist_pca")
    want_pi   <- isTRUE(input$show_pi) && is_cubist

    withProgress(message = "Predicting ...", value = 0, {
      result <- tryCatch({
        incProgress(0.3, detail = "Running model ...")
        p <- predict_soil(df, family_id = input$family_id,
                          properties = props, model_dir = model_dir,
                          with_pi    = want_pi,
                          confidence = 0.68)
        incProgress(0.5)
        p
      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8)
        NULL
      })
      pred_df(result)
      incProgress(1)
    })

    if (!is.null(pred_df())) {
      if (input$check_domain && length(props) >= 1) {
        ref_prop <- props[1]
        dom <- tryCatch(
          predict_applicability(df, input$family_id, ref_prop,
                                model_dir = model_dir),
          error = function(e) NULL
        )
        domain_df(dom)
      }
      updateTabsetPanel(session, "tabs", selected = "Predictions")
      pred_cols <- grep("_lo$|_hi$", names(pred_df()), invert = TRUE,
                        value = TRUE)
      pred_cols <- setdiff(pred_cols, "Soil_ID")
      n_models <- sum(vapply(pred_df()[, pred_cols, drop = FALSE],
                             function(x) !all(is.na(x)), logical(1)))
      showNotification(
        paste0("Done. ", n_models, "/", length(props),
               " properties predicted."),
        type = "message", duration = 4)
    }
  })

  output$pred_status <- renderUI({
    p <- pred_df(); req(p)
    f <- fam()
    algo_class <- if (identical(f$algorithm, "cubist_pca")) "algo-cubist"
                   else "algo-vae"
    pred_cols  <- grep("_lo$|_hi$", names(p), invert = TRUE, value = TRUE)
    pred_cols  <- setdiff(pred_cols, "Soil_ID")
    n_ok  <- sum(vapply(p[, pred_cols, drop = FALSE],
                        function(x) !all(is.na(x)), logical(1)))
    n_all <- length(pred_cols)
    tagList(
      tags$p(
        span(class = "model-tag", paste(n_ok, "/", n_all, "properties")),
        " | ",
        span(class = algo_class, fam()$algorithm %||% "vae"),
        " | ",
        strong(fam()$label)
      )
    )
  })

  output$pred_table <- DT::renderDataTable({
    p <- pred_df(); req(p)
    DT::datatable(
      format_predictions(p),
      rownames  = FALSE,
      selection = "none",
      options   = list(
        pageLength = 15, scrollX = TRUE,
        dom = "Bfrtip",
        buttons = list("csv", "excel")
      ),
      extensions = "Buttons"
    )
  })

  output$domain_ui <- renderUI({
    dom <- domain_df(); req(dom)
    method <- attr(dom, "method") %||% "unknown"
    n_in   <- sum(dom$in_domain, na.rm = TRUE)
    n_tot  <- nrow(dom)
    pct    <- round(100 * n_in / n_tot, 1)
    subtitle <- if (identical(method, "q_statistic"))
      "Q-statistic on PCA residuals (Safanelli 2025)"
    else
      "Mahalanobis distance in latent space"
    tagList(
      tags$h5("Applicability Domain"),
      tags$p(em(subtitle)),
      tags$p(
        strong(n_in), " / ", n_tot, " samples within domain (",
        pct, "%)"
      ),
      renderPlot(plot_applicability(dom), height = 300)
    )
  })

  output$dl_excel <- downloadHandler(
    filename = function() {
      paste0("autoSpectra_", input$family_id, "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      p <- pred_df()
      req(p)
      writexl::write_xlsx(format_predictions(p), path = file)
    }
  )
}

shinyApp(ui = ui, server = server)
