##' Function to compute the Bayes factors from MCMC samples.
##'
##' Computes the Bayes factors using \code{method} with respect to
##' \code{reference}.
##' @title Computation of Bayes factors at the skeleton points
##' @param runs A list with outputs from the function
##'   \code{\link{mcsglmm}} or \code{\link{mcstrga}}.
##' @param bfsize1 A scalar or vector of the same length as
##'   \code{runs} with all integer values or all values in (0, 1]. How
##'   many samples (or what proportion of the sample) to use for
##'   estimating the Bayes factors at the first stage. The remaining
##'   sample will be used for estimating the Bayes factors in the
##'   second stage. Setting it to 1 will perform only the first stage.
##' @param method Which method to use to calculate the Bayes factors:
##' Reverse logistic or Meng-Wong.
##' @param reference Which model goes in the denominator.
##' @param transf Whether to use a transformed sample for the
##'   computations. If \code{"no"} or \code{FALSE}, it doesn't. If
##'   \code{"mu"} or \code{TRUE}, it uses the samples for the mean. If
##'   \code{"wo"} it uses an alternative transformation. The latter
##'   can be used only for the families indicated by
##'   \code{.geoBayes_models$haswo}.
##' @return A list with components
##' \itemize{
##' \item \code{logbf} A vector containing logarithm of the Bayes factors.
##' \item \code{logLik1} \code{logLik2} Matrices with the values of
##' the log-likelihood computed from the samples for each model at the
##' first and second stages.
##' \item \code{isweights} A vector with the importance sampling
##' weights for computing the Bayes factors at new points that will be
##' used at the second stage. Used internally in
##' \code{\link{bf2new}} and \code{\link{bf2optim}}.
##' \item \code{controlvar} A matrix with the control variates
##' computed at the samples that will be used in the second stage.
##' \item \code{sample2} The MCMC sample for mu or z that will be
##' used in the second stage. Used internally in
##' \code{\link{bf2new}} and \code{\link{bf2optim}}.
##' \item \code{N1}, \code{N2} Vectors containing the sample sizes
##' used in the first and second stages.
##' \item \code{distmat} Matrix of distances between locations.
##' \item \code{betm0}, \code{betQ0}, \code{ssqdf}, \code{ssqsc},
##' \code{tsqdf}, \code{tsqsc}, \code{dispersion}, \code{response},
##' \code{weights}, \code{modelmatrix}, \code{locations},
##' \code{family}, \code{corrfcn}, \code{transf} Model parameters used
##' internally in.
##' \code{\link{bf2new}} and \code{\link{bf2optim}}.
##' \item \code{pnts} A list containing the skeleton points. Used
##' internally in \code{\link{bf2new}} and \code{\link{bf2optim}}.
##' }
##' @references Geyer, C. J. (1994). Estimating normalizing constants
##' and reweighting mixtures. Technical report, University of
##' Minnesota.
##'
##' Meng, X. L., & Wong, W. H. (1996). Simulating ratios of
##' normalizing constants via a simple identity: A theoretical
##' exploration. \emph{Statistica Sinica}, 6, 831-860.
##'
##' Roy, V., Evangelou, E., and Zhu, Z. (2015). Efficient estimation
##' and prediction for the Bayesian spatial generalized linear mixed
##' model with flexible link functions. \emph{Biometrics}, 72(1), 289-298.
##' @examples \dontrun{
##' data(rhizoctonia)
##' ### Define the model
##' corrf <- "spherical"
##' kappa <- 0
##' ssqdf <- 1
##' ssqsc <- 1
##' betm0 <- 0
##' betQ0 <- .01
##' family <- "binomial.probit"
##' ### Skeleton points
##' philist <- c(100, 140, 180)
##' omglist <- c(.5, 1)
##' parlist <- expand.grid(linkp=0, phi=philist, omg=omglist, kappa = kappa)
##' ### MCMC sizes
##' Nout <- 100
##' Nthin <- 1
##' Nbi <- 0
##' ### Take MCMC samples
##' runs <- list()
##' for (i in 1:NROW(parlist)) {
##'   runs[[i]] <- mcsglmm(Infected ~ 1, family, rhizoctonia, weights = Total,
##'                        atsample = ~ Xcoord + Ycoord,
##'                        Nout = Nout, Nthin = Nthin, Nbi = Nbi,
##'                        betm0 = betm0, betQ0 = betQ0,
##'                        ssqdf = ssqdf, ssqsc = ssqsc,
##'                        phi = parlist$phi[i], omg = parlist$omg[i],
##'                        linkp = parlist$linkp[i], kappa = parlist$kappa[i],
##'                        corrfcn = corrf,
##'                        corrtuning=list(phi = 0, omg = 0, kappa = 0))
##' }
##' bf <- bf1skel(runs)
##' bf$logbf
##' }
##' @importFrom sp spDists
##' @useDynLib geoBayes bfsp_no bfsp_mu bfsp_wo bfsp_tr
##' @export
bf1skel <- function(runs, bfsize1 = 0.80, method = c("RL", "MW"),
                    reference = 1, transf = c("no", "mu", "wo"))
{
  method <- match.arg(method)
  imeth <- match(method, eval(formals()$method))
  if (!all(sapply(runs, inherits, what = "geomcmc"))) {
    stop ("Input runs is not a list with elements of class geomcmc.")
  }
  nruns <- length(runs)
  if (nruns == 0) stop ("No runs specified")
  reference <- as.integer(reference)
  if (isTRUE(reference < 1L | reference > nruns)) {
    stop("Argument reference does not correspond to a run in runs.")
  }
  Nout <- sapply(runs, function(x) x$MCMC$Nout)
  Nout1 <- getsize(bfsize1, Nout, "*")
  Ntot1 <- sum(Nout1)
  Nout2 <- Nout - Nout1
  Ntot2 <- sum(Nout2)

  ## Check if fixed phi and omg
  if (!all(sapply(runs, function(x) length(x$FIXED$phi) == 1))) {
    stop("Each input runs must have a fixed value phi.")
  }
  if (!all(sapply(runs, function(x) length(x$FIXED$omg) == 1))) {
    stop("Each input runs must have a fixed value omg.")
  }

  ## Extract data and model
  nm_DATA <- c("response", "weights", "modelmatrix", "offset", "locations",
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
  offset <- DATA$offset
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

  ## Choose sample
  getsample <-
    transfsample(runs,
                 list(response = y, family = family), transf)
  sample <- matrix(unlist(getsample$sample), n)
  itr <- getsample$itr
  transf <- getsample$transf
  real_transf <- getsample$real_transf
  ifam <- getsample$ifam

  ## Skeleton points
  phi_pnts <- as.double(sapply(runs, function(r) r$FIXED$phi))
  omg_pnts <- as.double(sapply(runs, function(r) r$FIXED$omg))
  nu_pnts <- as.double(sapply(runs, function(r) r$FIXED$linkp_num))
  if (.geoBayes_corrfcn$needkappa[icf]) {
    kappa_pnts <- sapply(runs, function(r) r$FIXED$kappa)
    kappa_pnts <- .geoBayes_getkappa(kappa_pnts, icf)
  } else {
    kappa_pnts <- rep(0, nruns)
  }

  bfroutine <- paste0("bfsp_", real_transf)

  if (nruns == 1) {
    MCMC <- runs[[1]]$MCMC
    out <- list(logbf = 1, logLik1 = MCMC$logLik[1:Ntot1],
                logLik2 = MCMC$logLik[-(1:Ntot1)],
                isweights = rep.int(0, Ntot2),
                controlvar = matrix(1, Ntot2, 1),
                z = sample[[1]][, -(1:Ntot1), drop = FALSE],
                N1 = Nout1, N2 = Nout2,
                betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf,
                ssqsc = ssqsc, tsqdf = tsqdf, tsqsc = tsqsc,
                dispersion = dispersion, response = y,
                weights = l, modelmatrix = F, offset = offset, 
                locations = loc, longlat = DATA$longlat,
                distmat = dm,
                family = family,
                referencebf = 0, corrfcn = corrfcn, transf = transf,
                real_transf = real_transf, itr = itr,
                pnts = list(nu = nu_pnts, phi = phi_pnts, omg = omg_pnts,
                  kappa = kappa_pnts))
    return(out)
  }

  ## Split the sample
  sel1 <- rep(rep(c(TRUE, FALSE), nruns), rbind(Nout1, Nout2))
  z1 <- sample[, sel1, drop = FALSE]
  z2 <- sample[, !sel1, drop = FALSE]

  logbf <- numeric(nruns)
  lglk1 <- matrix(0., Ntot1, nruns)
  lglk2 <- matrix(0., Ntot2, nruns)
  zcv <- matrix(0., Ntot2, nruns)
  weights <- numeric(Ntot2)
  if (ifam == 0) {
    tsq <- tsqsc
  } else {
    tsq <- dispersion
  }
  RUN <- .Fortran(bfroutine,
                  weights = weights,
                  zcv = zcv,
                  logbf = logbf,
                  lglk1 = lglk1,
                  lglk2 = lglk2,
                  as.double(phi_pnts), as.double(omg_pnts),
                  as.double(nu_pnts), as.double(z1),
                  as.integer(Nout1), as.integer(Ntot1),
                  as.double(z2), as.integer(Nout2), as.integer(Ntot2),
                  as.double(y), as.double(l), as.double(F), as.double(offset), 
                  as.double(dm), as.double(betm0), as.double(betQ0),
                  as.double(ssqdf), as.double(ssqsc), max(tsqdf, 0),
                  as.double(tsq), as.double(kappa_pnts), as.integer(icf),
                  as.integer(n), as.integer(p), as.integer(nruns),
                  as.integer(ifam), as.integer(imeth), as.integer(itr),
                  PACKAGE = "geoBayes")
  refbf <- RUN$logbf[reference]
  logbf <- RUN$logbf - refbf
  if (Ntot2 > 0) {
    weights <- RUN$weights
    lglk2 <- RUN$lglk2
    zcv <- RUN$zcv
  } else {
    weights <- lglk2 <- zcv <- NULL
  }
  out <- list(logbf = logbf, logLik1 = RUN$lglk1, logLik2 = lglk2,
              isweights = weights, controlvar = zcv, sample2 = z2,
              N1 = Nout1, N2 = Nout2,
              betm0 = betm0,
              betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc, tsqdf = tsqdf,
              tsqsc = tsqsc, dispersion = dispersion, response = y, weights = l,
              modelmatrix = F, offset = offset, 
              locations = loc, distmat = dm, family = family,
              corrfcn = corrfcn, transf = transf,
              real_transf = real_transf, itr = itr,
              pnts = list(nu = nu_pnts, phi = phi_pnts, omg = omg_pnts,
                kappa = kappa_pnts))
  out
}
