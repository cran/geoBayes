##' Spatial joint log likelihood
##'
##' Computes the joint log likelihood log f(y,z|parameters) or log
##' f(y,mu|parameters) for each (y,z) or (y,mu) in runs and for parameters
##'   in pargrid up to a constant which does not depend on the parameters.
##' @title Spatial log likelihood
##' @param pargrid A data frame with components "linkp", "phi", "omg",
##' "kappa". Each row gives a combination of the parameters to compute
##' the log-likelihood.
##' @param runs A list with outputs from the function
##' \code{\link{mcsglmm}} or \code{\link{mcstrga}}.
##' @param transf Whether to use a transformed sample for the
##'   computations. If \code{"no"} or \code{FALSE}, it doesn't. If
##'   \code{"mu"} or \code{TRUE}, it uses the samples for the mean. If
##'   \code{"wo"} it uses an alternative transformation. The latter
##'   can be used only for the families indicated by
##'   \code{.geoBayes_models$haswo}.
##' @return A matrix with number of rows the total number of samples
##'   in runs and number of columns the number of rows in pargrid. The
##'   [i,j] element of the matrix is the value of the loglikelihood at
##'   the ith sample when all samples in runs are put together evaluated
##'   at the jth parameter value.
##' @useDynLib geoBayes llikfcn_no llikfcn_mu llikfcn_wo llikfcn_tr
##' @export
sploglik <- function(pargrid, runs, transf = c("no", "mu", "wo"))
{
  ## Check input
  if (any(sapply(runs, class) != "geomcmc")) {
    stop ("Input runs is not a list with elements of class geomcmc")
  }
  nruns <- length(runs)
  if (nruns == 0) stop ("No runs specified")
  Nout <- sapply(runs, "[[", "Nout")
  Ntot <- sum(Nout)

  ## Extract model
  modelvars <- c("response", "weights", "modelmatrix", "family",
                 "betm0", "betQ0", "ssqdf", "ssqsc",
                 "dispersion", "tsqdf", "tsqsc", "locations",
                 "longlat", "corrfcn")
  models <- lapply(runs, "[", modelvars)
  model <- models[[1]]
  if (nruns > 1 && !all(sapply(models[2:nruns], identical, model))) {
    stop("MCMC chains don't all correspond to the same model")
  }
  y <- model$response
  n <- length(y)
  l <- model$weights
  F <- model$modelmatrix
  p <- NCOL(F)
  family <- model$family
  ## ifam <- .geoBayes_family(family)
  betm0 <- model$betm0
  betQ0 <- model$betQ0
  ssqdf <- model$ssqdf
  ssqsc <- model$ssqsc
  dispersion <- model$dispersion
  tsqdf <- model$tsqdf
  tsqsc <- model$tsqsc
  corrfcn <- model$corrfcn
  icf <- .geoBayes_correlation(corrfcn)
  loc <- model$locations
  dm <- sp::spDists(loc, longlat = model$longlat)

  ## Extract the grid for the new parameters
  pargrid <- check_pargrid(pargrid, family, corrfcn)
  phi <- pargrid$phi
  omg <- pargrid$omg
  kappa <- pargrid$kappa
  nu <- pargrid$nu
  kg <- NROW(pargrid)

  ## Transformed sample
  getsample <- transfsample(runs, model, transf)
  sample <- matrix(unlist(getsample$sample), n)
  itr <- getsample$itr
  transf <- getsample$transf
  real_transf <- getsample$real_transf
  ifam <- getsample$ifam
  froutine <- paste0("llikfcn_", real_transf)

  if (ifam == 0) {
    tsq <- tsqsc
  } else {
    tsq <- dispersion
    tsqdf <- 0
  }

  lglk <- matrix(0, Ntot, kg)
  fcall <- .Fortran(froutine,
                    lglk = lglk,
                    as.double(phi), as.double(omg), as.double(nu),
                    as.double(kappa),
                    as.double(sample), as.integer(Ntot), as.double(y),
                    as.double(l), as.double(F), as.double(dm),
                    as.double(betm0), as.double(betQ0), as.double(ssqdf),
                    as.double(ssqsc), as.double(tsqdf), as.double(tsq),
                    as.integer(icf), as.integer(n), as.integer(p),
                    as.integer(kg), as.integer(ifam), as.integer(itr),
                    PACKAGE = "geoBayes")
  lglk <- fcall$lglk
  lglk
}
