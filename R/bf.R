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
##' linkp <- "probit"
##' ### Skeleton points
##' philist <- c(100, 140, 180)
##' omglist <- c(.5, 1)
##' parlist <- expand.grid(phi=philist, linkp=linkp, omg=omglist, kappa = kappa)
##' ### MCMC sizes
##' Nout <- 100
##' Nthin <- 1
##' Nbi <- 0
##' ### Take MCMC samples
##' runs <- list()
##' for (i in 1:NROW(parlist)) {
##'   runs[[i]] <- mcsglmm(Infected ~ 1, 'binomial', rhizoctonia, weights = Total,
##'                        atsample = ~ Xcoord + Ycoord,
##'                        Nout = Nout, Nthin = Nthin, Nbi = Nbi,
##'                        betm0 = betm0, betQ0 = betQ0,
##'                        ssqdf = ssqdf, ssqsc = ssqsc,
##'                        phistart = parlist$phi[i], omgstart = parlist$omg[i],
##'                        linkp = parlist$linkp[i], kappa = parlist$kappa[i],
##'                        corrfcn = corrf, phisc = 0, omgsc = 0)
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
  classes <- sapply(runs, class)
  if (any(classes != "geomcmc")) {
    stop ("Input runs is not a list with elements of class geomcmc")
  }
  nruns <- length(runs)
  if (nruns == 0) stop ("No runs specified")
  reference <- as.integer(reference)
  if (isTRUE(reference < 1L | reference > nruns)) {
    stop("Argument reference does not correspond to a run in runs")
  }
  Nout <- sapply(runs, "[[", "Nout")
  Nout1 <- getsize(bfsize1, Nout, "*")

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
  n <- as.integer(length(y))
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

  ## Choose sample
  getsample <- transfsample(runs, model, transf)
  sample <- matrix(unlist(getsample$sample), n)
  itr <- getsample$itr
  transf <- getsample$transf
  real_transf <- getsample$real_transf
  ifam <- getsample$ifam

  bfroutine <- paste0("bfsp_", real_transf)

  Ntot1 <- sum(Nout1)
  Nout2 <- Nout - Nout1
  Ntot2 <- sum(Nout2)
  if (nruns == 1) {
    runs <- runs[[1]]
    out <- list(logbf = 1, logLik1 = runs$logLik[1:Ntot1],
                logLik2 = runs$logLik[-(1:Ntot1)],
                isweights = rep.int(0, Ntot2),
                controlvar = matrix(1, Ntot2, 1),
                z = sample[[1]][, -(1:Ntot1), drop = FALSE],
                N1 = Nout1, N2 = Nout2,
                betm0 = runs$betm0, betQ0 = runs$betQ0, ssqdf = runs$ssqdf,
                ssqsc = runs$ssqsc, tsqdf = runs$tsqdf, tsqsc = runs$tsqsc,
                dispersion = runs$dispersion, response = runs$response,
                weights = runs$weights, modelmatrix = runs$modelmatrix,
                locations = runs$locations, longlat = runs$longlat,
                distmat = sp::spDists(runs$locations),
                family = runs$family,
                referencebf = 0, corrfcn = runs$corrfcn, transf = transf,
                real_transf = real_transf, itr = itr,
                pnts = list(nu = runs$nu, phi = runs$phi, omg = runs$omg,
                  kappa = runs$kappa))
    return(out)
  }
  icf <- .geoBayes_correlation(corrfcn)
  loc <- model$locations
  dm <- sp::spDists(loc, longlat = model$longlat)
  fixphi <- sapply(runs, function(r) attr(r[["phi"]], "fixed"))
  if (sum(fixphi) != 0 & sum(fixphi) != nruns) {
    stop ("The parameter phi is not consistently fixed or estimated")
  }
  fixphi <- fixphi[1]
  if (!fixphi) {
    stop("The case where phi is not fixed is not yet implemented")
  }
  fixomg <- sapply(runs, function(r) attr(r[["omg"]], "fixed"))
  if (sum(fixomg) != 0 & sum(fixomg) != nruns) {
    stop ("The parameter omg is not consistently fixed or estimated")
  }
  fixomg <- fixomg[1]
  if (!fixomg) {
    stop("The case where omg is not fixed is not yet implemented")
  }
  fixnu <- sapply(runs, function(r) attr(r[["nu"]], "fixed"))
  if (sum(fixnu) != 0 & sum(fixnu) != nruns) {
    stop ("The parameter nu is not consistently fixed or estimated")
  }
  fixnu <- fixnu[1]
  if (!fixnu) {
    stop("The parameter nu must be fixed in the MCMC")
  }
  if (fixphi) {
    phi_pnts <- sapply(runs, function(r) r[["phi"]][1])
  }
  if (fixomg) {
    omg_pnts <- sapply(runs, function(r) r[["omg"]][1])
  }
  if (fixnu) {
    nu_pnts <- sapply(runs, function(r) r[["nu"]][1])
  }
  if (.geoBayes_corrfcn$needkappa[icf]) {
    kappa_pnts <- sapply(runs, function(r) r[["kappa"]][1])
  } else {
    kappa_pnts <- rep(0, nruns)
  }
  kappa_pnts <- as.double(kappa_pnts)
  if (corrfcn %in% c("matern", "powerexponential") && any(kappa_pnts < 0)) {
    stop ("Argument kappa_pnts cannot be negative")
  }
  if (corrfcn == "powerexponential" && any(kappa_pnts > 2)) {
    stop ("Argument kappa_pnts cannot be more than 2")
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
                  as.double(y), as.double(l), as.double(F),
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
              modelmatrix = F, locations = loc, distmat = dm, family = family,
              corrfcn = corrfcn, transf = transf,
              real_transf = real_transf, itr = itr,
              pnts = list(nu = nu_pnts, phi = phi_pnts, omg = omg_pnts,
                kappa = kappa_pnts))
  out
}
