#' Spatial correlation used in the geoBayes package
#'
#' This hidden variable contains a choice of correlation functions
#' that can be fit with this package. The function can be chosen in
#' the \code{corrfcn} input of the relevant function. This variable
#' cannot be overwritten.
#' @name geoBayes_correlation
#' @export
.geoBayes_corrfcn <- data.frame(
  corrfcn   = c("matern", "spherical", "powerexponential"),
  needkappa = c(    TRUE,       FALSE,               TRUE),
  validklb  = c(       0,          NA,                  0),
  validkub  = c(     Inf,          NA,                  2),
  stringsAsFactors = FALSE)

lockBinding(".geoBayes_corrfcn", environment())

.geoBayes_correlation <- function (corrfcn)
{
  i <- pmatch(corrfcn, .geoBayes_corrfcn$corrfcn)
  if (is.na(i)) stop ("Cannot deduce correlation.")
  i
}



## Checks if the value of kappa is valid for the given correlation
## function. Returns the value of kappa.
.geoBayes_getkappa <- function (kappa, icf) {
  if (length(icf) != 1) stop ("Multiple correlations inputted.")
  if (!is.numeric(icf)) icf <- .geoBayes_correlation(icf)
  if (.geoBayes_corrfcn$needkappa[icf]) {
    kappa <- as.double(kappa)
    if (length(kappa) < 1) stop ("No kappa value inputted.")
    if (is.finite(.geoBayes_corrfcn$validklb[icf]) &&
        any(kappa <= .geoBayes_corrfcn$validklb[icf])) {
      stop (paste0("Input kappa cannot be less than or equal to ",
                   .geoBayes_corrfcn$validklb[icf], "."))
    }
    if (is.finite(.geoBayes_corrfcn$validkub[icf]) &&
        any(kappa > .geoBayes_corrfcn$validkub[icf])) {
      stop (paste0("Input kappa cannot be greater than ",
                   .geoBayes_corrfcn$validkub[icf], "."))
    }
    kappa
  } else {
    0
  }
}
