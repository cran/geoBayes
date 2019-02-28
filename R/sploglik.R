##' Spatial joint log likelihood
##'
##' Computes the joint log likelihood log f(y,T(z)|parameters) where
##' T(z) is the transformation, for each (y,z) in runs and for parameters
##' in pargrid up to a constant which does not depend on the
##' parameters. The parameters beta and sigma^2 are integrated out.
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
  Nout <- sapply(runs, function(x) x$MCMC$Nout)
  Ntot <- sum(Nout)

  ## Extract data and model
  nm_DATA <- c("response", "weights", "modelmatrix", "locations",
               "longlat")
  nm_MODEL <- c("family", "corrfcn", "betm0", "betQ0", "ssqdf", "ssqsc",
                "tsqdf", "tsqsc", "dispersion")
  DATA <- runs[[1]]$DATA[nm_DATA]
  MODEL <- runs[[1]]$MODEL[nm_MODEL]
  if (nruns > 1) {
    for (i in 2:nruns) {
      if (!identical(runs[[i]]$DATA[nm_DATA], DATA)) {
        stop("MCMC chains don't all correspond to the same data.")
      }
      if (!identical(runs[[i]]$MODEL[nm_MODEL], MODEL)) {
        stop("MCMC chains don't all correspond to the same model.")
      }
    }
  }
  y <- DATA$response
  n <- as.integer(length(y))
  l <- DATA$weights
  F <- DATA$modelmatrix
  p <- NCOL(F)
  loc <- DATA$locations
  dm <- sp::spDists(loc, longlat = DATA$longlat)
  family <- MODEL$family
  ## ifam <- .geoBayes_family(family)
  corrfcn <- MODEL$corrfcn
  icf <- .geoBayes_correlation(corrfcn)
  betm0 <- MODEL$betm0
  betQ0 <- MODEL$betQ0
  ssqdf <- MODEL$ssqdf
  ssqsc <- MODEL$ssqsc
  tsqdf <- MODEL$tsqdf
  tsqsc <- MODEL$tsqsc
  dispersion <- MODEL$dispersion

  ## Extract the grid for the new parameters
  pargrid <- check_pargrid(pargrid, family, corrfcn)
  phi <- pargrid$phi
  omg <- pargrid$omg
  kappa <- pargrid$kappa
  nu <- pargrid$nu
  kg <- NROW(pargrid)

  ## Transformed sample
  getsample <-
    transfsample(runs,
                 list(response = y, family = family), transf)
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

##' Spatial joint log likelihood
##'
##' Computes the joint log likelihood log f(y,T(z)|parameters) where
##' T(z) is the transformation, for each (y,z) in runs and for parameters
##' in runs up to a constant which does not depend on the parameters.
##' The parameters beta and sigma^2 are integrated out. 
##' @title Spatial log likelihood
##' @param runs A list with outputs from the function
##' \code{\link{mcsglmm}} or \code{\link{mcstrga}}.
##' @param transf Whether to use a transformed sample for the
##'   computations. If \code{"no"} or \code{FALSE}, it doesn't. If
##'   \code{"mu"} or \code{TRUE}, it uses the samples for the mean. If
##'   \code{"wo"} it uses an alternative transformation. The latter
##'   can be used only for the families indicated by
##'   \code{.geoBayes_models$haswo}. The input can also be a vector
##'   (of the same length as \code{runs} to allow for different
##'   transformation to be used when evaluating each likelihood.
##' @return A matrix with number of rows the total number of samples
##'   in runs and number of columns the length of \code{runs}. The
##'   [i,j] element of the matrix is the value of the loglikelihood at
##'   the ith sample when all samples in \code{runs} are put together evaluated
##'   at the jth parameter value.
##' @useDynLib geoBayes llikfcn_no llikfcn_mu llikfcn_wo llikfcn_tr
##' @export
sploglik_cross <- function(runs, transf = c("no", "mu", "wo"))
{
  ## Check input
  if (any(sapply(runs, class) != "geomcmc")) {
    stop ("Input runs is not a list with elements of class geomcmc")
  }
  nruns <- length(runs)
  if (nruns == 0) stop ("No runs specified")
  Nout <- sapply(runs, function(x) x$MCMC$Nout)
  Ntot <- sum(Nout)

  ## Extract data and model
  nm_DATA <- c("response", "weights", "modelmatrix", "locations",
               "longlat")
  nm_MODEL <- c("betm0", "betQ0", "ssqdf", "ssqsc",
                "tsqdf", "tsqsc", "dispersion")
  DATA <- runs[[1]]$DATA[nm_DATA]
  MODEL <- runs[[1]]$MODEL[nm_MODEL]
  if (nruns > 1) {
    for (i in 2:nruns) {
      if (!identical(runs[[i]]$DATA[nm_DATA], DATA)) {
        stop("MCMC chains don't all correspond to the same data.")
      }
      if (!identical(runs[[i]]$MODEL[nm_MODEL], MODEL)) {
        stop("MCMC chains don't all correspond to the same model.")
      }
    }
  }
  y <- DATA$response
  n <- as.integer(length(y))
  l <- DATA$weights
  F <- DATA$modelmatrix
  p <- NCOL(F)
  loc <- DATA$locations
  dm <- sp::spDists(loc, longlat = DATA$longlat)
  corrfcn <- MODEL$corrfcn
  betm0 <- MODEL$betm0
  betQ0 <- MODEL$betQ0
  ssqdf <- MODEL$ssqdf
  ssqsc <- MODEL$ssqsc
  tsqdf <- MODEL$tsqdf
  tsqsc <- MODEL$tsqsc
  dispersion <- MODEL$dispersion

  ifam_list <- sapply(runs, function(r) .geoBayes_family(r$MODEL$family))
  if (any(ifam_list == 0) && any(ifam_list > 0))
    stop ("Cannot mix transformed Gaussian with other families.")
  icf_list <- sapply(runs, function(r) .geoBayes_correlation(r$MODEL$corrfcn))
  nu_list <- sapply(runs, function(r) r$FIXED$linkp_num)
  kappa_list <- sapply(runs, function(r) r$FIXED$kappa)
  phi_list <- sapply(runs, function(r) r$FIXED$phi)
  omg_list <- sapply(runs, function(r) r$FIXED$omg)
  if (length(phi_list) == 0) {
    phi_mcmc <- TRUE
    phi_sample <- do.call("c", lapply(runs, function(r) r$MCMC$phi))
    phi_pars <- sapply(runs, function(r) r$MODEL$phipars)
  } else if (length(phi_list) == nruns) {
    phi_mcmc <- FALSE
  } else {
    stop ("Models are not consistent w.r.t fixed/random phi.")
  }
  if (length(omg_list) == 0) {
    omg_mcmc <- TRUE
    omg_sample <- do.call("c", lapply(runs, function(r) r$MCMC$omg))
    omg_pars <- sapply(runs, function(r) r$MODEL$omgpars)
  } else if (length(omg_list) == nruns) {
    omg_mcmc <- FALSE
  } else {
    stop ("Models are not consistent w.r.t fixed/random omg.")
  }

  ## Transformed sample
  sample <- array(, c(n, Ntot, nruns))
  itr <- array(, c(n, nruns))
  transf <- rep(transf, length = nruns)
  for (i in 1:nruns) {
    getsample <-
      transfsample(runs, list(response = y, family = ifam_list[i]), transf[[i]],
                   verb = TRUE)
    sample[, , i] <- do.call("c", getsample$sample)
    itr[, i] <- getsample$itr
  }

  if (ifam_list[1] == 0) {
    tsq <- tsqsc
  } else {
    tsq <- dispersion
    tsqdf <- 0
  }

  lglk <- matrix(0, Ntot, nruns)

  if (phi_mcmc && omg_mcmc) {
    fcall <- .Fortran("llikfcnmc_11",
                      lglk = lglk,
                      as.double(nu_list),
                      as.double(kappa_list),
                      as.double(sample),
                      as.double(phi_sample),
                      as.double(omg_sample),
                      as.double(phi_pars),
                      as.double(omg_pars),
                      as.integer(Ntot), as.double(y),
                      as.double(l), as.double(F), as.double(dm),
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(tsqdf), as.double(tsq),
                      as.integer(icf_list), as.integer(n), as.integer(p),
                      as.integer(nruns), as.integer(ifam_list), as.integer(itr),
                      PACKAGE = "geoBayes")
  } else if (!phi_mcmc && omg_mcmc) {
    fcall <- .Fortran("llikfcnmc_01",
                      lglk = lglk,
                      as.double(nu_list),
                      as.double(phi_list),
                      as.double(kappa_list),
                      as.double(sample),
                      as.double(omg_sample),
                      as.double(omg_pars),
                      as.integer(Ntot), as.double(y),
                      as.double(l), as.double(F), as.double(dm),
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(tsqdf), as.double(tsq),
                      as.integer(icf_list), as.integer(n), as.integer(p),
                      as.integer(nruns), as.integer(ifam_list), as.integer(itr),
                      PACKAGE = "geoBayes")
  } else if (phi_mcmc && !omg_mcmc) {
    fcall <- .Fortran("llikfcnmc_10",
                      lglk = lglk,
                      as.double(nu_list),
                      as.double(omg_list),
                      as.double(kappa_list),
                      as.double(sample),
                      as.double(phi_sample),
                      as.double(phi_pars),
                      as.integer(Ntot), as.double(y),
                      as.double(l), as.double(F), as.double(dm),
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(tsqdf), as.double(tsq),
                      as.integer(icf_list), as.integer(n), as.integer(p),
                      as.integer(nruns), as.integer(ifam_list), as.integer(itr),
                      PACKAGE = "geoBayes")
  } else if (!phi_mcmc && !omg_mcmc) {
    fcall <- .Fortran("llikfcnmc_00",
                      lglk = lglk,
                      as.double(nu_list),
                      as.double(phi_list),
                      as.double(omg_list),
                      as.double(kappa_list),
                      as.double(sample),
                      as.integer(Ntot), as.double(y),
                      as.double(l), as.double(F), as.double(dm),
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(tsqdf), as.double(tsq),
                      as.integer(icf_list), as.integer(n), as.integer(p),
                      as.integer(nruns), as.integer(ifam_list), as.integer(itr),
                      PACKAGE = "geoBayes")
  }
  lglk <- fcall$lglk
  lglk
}
