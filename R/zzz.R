# R/zzz.R — package-level declarations

# Suppress R CMD check NOTEs for ggplot2 NSE column names
# These are column names created inside functions via data.frame() and passed
# to ggplot2::aes() using non-standard evaluation.
utils::globalVariables(c(
  # plot_spectra
  "wl", "value", "Soil_ID", ".data",
  # plot_predictions
  "obs", "pred",
  # plot_applicability
  "mahal_dist", "in_domain",
  # plot_mean_spectrum
  "lo", "hi"
))
