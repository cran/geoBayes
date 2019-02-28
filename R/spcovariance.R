##' @param ... Further arguments. Not currently in use.
##' @export
##' @name spcovariance
spcovariance <- function (...) {
  UseMethod("spcovariance")
}

##' Calculates the spatial variance-covariance matrix for a selection
##' of correlation functions.
##'
##' @title Spatial variance-covariance matrix
##' @param formula A formula of the form \code{~ Xcoord + Ycoord}
##'   specifying the sampled locations.
##' @param data An optional data frame containing the variables in the
##'   model.
##' @param subset An optional set of indices. The covariance will be
##'   calculated for those coordinates only.
##' @param corrfcn The correlation function to use.
##' @param ssq The partial sill parameter.
##' @param phi The spatial range parameter.
##' @param omg The relative nugget parameter.
##' @param kappa The spatial smoothness parameter.
##' @param longlat How to compute the distance between locations. If
##'   \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##'   distance. See \code{\link[sp]{spDists}}.
##' @return For a formula input, a variance-covariance matrix. For a
##'   numeric input, an object of the the same dimensions as its first
##'   input. 
##' @export
##' @importFrom sp spDists
##' @useDynLib geoBayes spcorr
##' @name spcovariance
spcovariance.formula <- function (formula, data, subset, corrfcn, ssq, 
                                  phi, omg, kappa, longlat = FALSE, ...)
{
  ## Correlation function
  icf <- .geoBayes_correlation(corrfcn)
  corrfcn <- .geoBayes_corrfcn$corrfcn[icf]
  if (!.geoBayes_corrfcn$needkappa[icf]) kappa <- 0

  ## Design matrix and data
  if (missing(data)) data <- environment(formula)
  mfc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(mfc), 0L)
  mfc <- mfc[c(1L, m)]
  atsample <- update(formula, NULL ~ . + 0) # No response and no intercept
  mfc$formula <- atsample
  mfc$drop.unused.levels <- TRUE
  mfc$na.action <- "na.pass"
  mfc[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfc, parent.frame())
  loc <- as.matrix(mf)
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations")
  if (corrfcn == "spherical" && NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }
  nloc <- dim(loc)[1]
  labels <- rownames(loc)
  dm <- sp::spDists(loc, longlat = longlat)
  dmup <- dm[upper.tri(dm)]
  spcorr <- spcovariance.numeric(dmup, corrfcn, ssq=ssq, phi=phi,
                                 omg=0, kappa=kappa)
  Sg <- diag(ssq*(1 + omg), nloc, nloc)
  Sg[upper.tri(Sg)] <- spcorr
  Sg[lower.tri(Sg)] <- t(Sg)[lower.tri(Sg)]
  if (is.null(labels)) {
    dimnames(Sg) <- list(seq_len(nloc), seq_len(nloc))
  } else {
      dimnames(Sg) <- list(labels, labels)
  }
  Sg
}

##' @param x A numerical object of distances.
##' @export
##' @useDynLib geoBayes spcorr
##' @name spcovariance
spcovariance.numeric <- function (x, corrfcn, ssq, 
                                  phi, omg, kappa, ...)
{
  ## Correlation function
  icf <- .geoBayes_correlation(corrfcn)
  if (!.geoBayes_corrfcn$needkappa[icf]) kappa <- 0

  xin <- as.vector(x)
  lenx <- length(xin)
  spcorr <- .Fortran("spcorr", as.double(xin), as.double(phi), as.double(kappa),
                     as.integer(lenx), icf, PACKAGE = "geoBayes")[[1]]
  out <- spcorr*ssq
  ii <- xin == 0
  if (omg > 0 && any(ii)) out[ii] <- out[ii] + ssq*omg
  x[] <- out
  x
}


##' @export
##' @useDynLib geoBayes spcorr
##' @name spcovariance
spcovariance.dist <- function (x, corrfcn, ssq, 
                               phi, omg, kappa, ...)
{
  ## Correlation function
  icf <- .geoBayes_correlation(corrfcn)
  if (!.geoBayes_corrfcn$needkappa[icf]) kappa <- 0

  nloc <- attr(x,"Size")
  labels <- attr(x, "Labels")
  x <- as.vector(x)
  spcorr <- spcovariance.numeric(x, corrfcn, ssq, phi, omg=0, kappa)
  Sg <- diag(ssq*(1 + omg), nloc, nloc)
  Sg[upper.tri(Sg)] <- spcorr
  Sg[lower.tri(Sg)] <- t(Sg)[lower.tri(Sg)]
  if (is.null(labels)) {
    dimnames(Sg) <- list(seq_len(nloc), seq_len(nloc))
  } else {
      dimnames(Sg) <- list(labels, labels)
  }
  Sg
}
