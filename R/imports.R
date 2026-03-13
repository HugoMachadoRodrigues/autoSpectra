# R/imports.R — Centralised importFrom declarations for autoSpectra
#
# roxygen2 reads these tags and writes the corresponding importFrom() lines
# to NAMESPACE. Keeping them here prevents R CMD check WARNINGs about
# "no visible binding" and "namespace not declared in DESCRIPTION".

#' @importFrom MASS ginv
#' @importFrom stats approx cov IQR predict qchisq quantile sd
#' @importFrom shinyWidgets pickerInput pickerOptions updatePickerInput
#' @importFrom utils globalVariables
NULL
