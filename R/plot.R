# R/plot.R — Spectral visualization functions for autoSpectra

#' Plot one or more spectra from a data frame
#'
#' @param df Data frame with Soil_ID column and numeric spectral columns
#' @param sample_ids Character vector of Soil_IDs to plot. If NULL, plots all.
#' @param id_col Name of the sample ID column (default "Soil_ID")
#' @param xlab X-axis label; auto-detected if NULL ("Wavelength (nm)" or
#'   "Wavenumber (cm-1)")
#' @param ylab Y-axis label (default "Response")
#' @param colour_by Name of a non-spectral column to colour lines by; NULL for
#'   a single colour
#' @param family Family list (from model_registry) used to draw the model grid
#'   range as dashed vertical lines. Pass NULL to skip.
#' @param title Plot title (default: "Spectra")
#' @param alpha Line transparency (default 0.7)
#' @return A ggplot2 object
#' @export
plot_spectra <- function(df,
                         sample_ids = NULL,
                         id_col     = "Soil_ID",
                         xlab       = NULL,
                         ylab       = "Response",
                         colour_by  = NULL,
                         family     = NULL,
                         title      = "Spectra",
                         alpha      = 0.7) {

  if (!is.null(sample_ids))
    df <- df[df[[id_col]] %in% sample_ids, , drop = FALSE]

  wl_info <- get_wavelengths(df, id_col = id_col)
  if (length(wl_info$wl) == 0)
    stop("No numeric spectral columns found in df.")

  # Build long format for ggplot
  spec_mat <- as.matrix(df[, wl_info$cols, drop = FALSE])
  ids      <- df[[id_col]]

  long <- data.frame(
    Soil_ID = rep(ids, times = length(wl_info$wl)),
    wl      = rep(wl_info$wl, each = length(ids)),
    value   = as.vector(spec_mat),
    stringsAsFactors = FALSE
  )

  # Optionally merge a colour variable
  if (!is.null(colour_by) && colour_by %in% names(df)) {
    meta <- df[, c(id_col, colour_by), drop = FALSE]
    names(meta)[1] <- "Soil_ID"
    long <- merge(long, meta, by = "Soil_ID", all.x = TRUE)
  }

  # Auto x-axis label
  if (is.null(xlab)) {
    xlab <- if (max(wl_info$wl) > 5000) "Wavenumber (cm\u207b\u00b9)" else "Wavelength (nm)"
  }

  # Base plot
  if (!is.null(colour_by) && colour_by %in% names(long)) {
    p <- ggplot2::ggplot(long, ggplot2::aes(
      x = wl, y = value,
      group = Soil_ID,
      colour = .data[[colour_by]]
    ))
  } else {
    p <- ggplot2::ggplot(long, ggplot2::aes(
      x = wl, y = value, group = Soil_ID
    ))
  }

  p <- p +
    ggplot2::geom_line(alpha = alpha, linewidth = 0.4) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_minimal(base_size = 12)

  # Draw model grid boundaries
  if (!is.null(family) && !is.null(family$wavegrid)) {
    p <- p + ggplot2::geom_vline(
      xintercept = range(family$wavegrid),
      linetype = "dashed", colour = "steelblue", linewidth = 0.6
    )
  }

  p
}

#' Plot predicted vs observed values for a soil property
#'
#' @param observed Numeric vector of observed values
#' @param predicted Numeric vector of predicted values
#' @param prop_label Display label for the property (e.g., "Organic Carbon (%)")
#' @param show_metrics Logical; add RMSE, R2, RPIQ annotation (default TRUE)
#' @return A ggplot2 object
#' @export
plot_predictions <- function(observed, predicted,
                             prop_label    = "Property",
                             show_metrics  = TRUE) {
  df_p <- data.frame(obs = observed, pred = predicted)
  df_p <- df_p[is.finite(df_p$obs) & is.finite(df_p$pred), ]

  lims <- range(c(df_p$obs, df_p$pred), na.rm = TRUE)

  p <- ggplot2::ggplot(df_p, ggplot2::aes(x = obs, y = pred)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5, colour = "#2166ac") +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey40") +
    ggplot2::coord_equal(xlim = lims, ylim = lims) +
    ggplot2::labs(x = paste("Observed:", prop_label),
                  y = paste("Predicted:", prop_label)) +
    ggplot2::theme_minimal(base_size = 12)

  if (show_metrics) {
    m <- metrics_from_y(df_p$obs, df_p$pred)
    ann <- sprintf("RMSE = %.3f\nR\u00b2 = %.3f\nRPIQ = %.2f\nn = %d",
                   m$RMSE, m$R2, m$RPIQ, nrow(df_p))
    p <- p + ggplot2::annotate("text",
      x = lims[1] + 0.05 * diff(lims),
      y = lims[2] - 0.05 * diff(lims),
      label = ann, hjust = 0, vjust = 1, size = 3.5,
      colour = "grey20"
    )
  }
  p
}

#' Plot latent-space applicability domain scores
#'
#' Mahalanobis distances are shown as a horizontal bar chart; samples outside
#' the 95% threshold are highlighted in red.
#'
#' @param app_df Data frame from predict_applicability()
#' @param title Plot title
#' @return A ggplot2 object
#' @export
plot_applicability <- function(app_df, title = "Applicability Domain") {
  app_df <- app_df[order(app_df$mahal_dist, decreasing = TRUE), ]
  app_df$Soil_ID <- factor(app_df$Soil_ID, levels = app_df$Soil_ID)

  ggplot2::ggplot(app_df, ggplot2::aes(
    x = mahal_dist, y = Soil_ID,
    fill = in_domain
  )) +
    ggplot2::geom_col(show.legend = TRUE) +
    ggplot2::geom_vline(xintercept = app_df$thr95[1],
                        linetype = "dashed", colour = "red", linewidth = 0.7) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#4dac26", "FALSE" = "#d01c8b"),
      labels = c("TRUE" = "In domain", "FALSE" = "Out of domain"),
      name   = NULL
    ) +
    ggplot2::labs(x = "Mahalanobis distance (latent space)",
                  y = NULL, title = title,
                  caption = "Dashed line = 95% chi-squared threshold") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")
}

#' Plot a mean spectrum with ±1 SD ribbon from a spectral matrix
#'
#' @param M Numeric matrix (rows = samples, cols = wavelengths)
#' @param wl Numeric vector of wavelength positions (same length as ncol(M))
#' @param xlab X-axis label (default "Wavelength (nm)")
#' @param ylab Y-axis label (default "Mean Response")
#' @param title Plot title
#' @param colour Line/ribbon colour (default "#2166ac")
#' @return A ggplot2 object
#' @export
plot_mean_spectrum <- function(M, wl,
                               xlab   = "Wavelength (nm)",
                               ylab   = "Mean Response",
                               title  = "Mean \u00b1 SD spectrum",
                               colour = "#2166ac") {
  mu  <- colMeans(M, na.rm = TRUE)
  sig <- apply(M, 2, stats::sd, na.rm = TRUE)
  df  <- data.frame(wl = wl, mu = mu, lo = mu - sig, hi = mu + sig)

  ggplot2::ggplot(df, ggplot2::aes(x = wl)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         fill = colour, alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = mu), colour = colour, linewidth = 0.7) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_minimal(base_size = 12)
}
