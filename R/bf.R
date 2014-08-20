##' Function to compute the Bayes factors from MCMC samples.
##'
##' Computes the bayes factors using \code{method} with respect to
##' \code{reference}. 
##' @title Computation of Bayes factors at the skeleton points
##' @param runs A list with outputs from the function
##' \code{\link{mcsglmm}} or \code{\link{mcstrga}}.
##' @param bfsize1 A scalar or vector of the same length as \code{runs}
##' with all integer values or all values in (0, 1]. How many samples
##' (or what proportion of the sample) to use for estimating the Bayes
##' factors at the first stage. The remaining sample will be used for
##' estimating the Bayes factors in the second stage. Setting it to 1
##' will perform only the first stage.
##' @param method Which method to use to calculate the Bayes factors:
##' Reverse logistic or Meng-Wong.
##' @param reference Which model goes in the denominator.
##' @return A list with components
##' \itemize{
##' \item \code{logbf} A vector containing logarithm of the Bayes factors.
##' \item \code{logLik1} \code{logLik2} Matrices with the values of
##' the log-likelihood computed from the samples for each model at the
##' first and second stages.
##' }
##' @references Geyer, C. J. (1994). Estimating Normalizing Constants
##' and Reweighting Mixtures. Technical report, University of
##' Minnesota.
##'  
##' Meng, X. L., & Wong, W. H. (1996). Simulating Ratios of
##' Normalizing Constants via a Simple Identity: A Theoretical
##' Exploration. \emph{Statistica Sinica}, 6, 831-860. 
##' @importFrom sp spDists
##' @export 
bf1skel <- function(runs, bfsize1 = 0.80, method = c("RL", "MW"),
                    reference = 1){
  method <- match.arg(method)
  imeth <- match(method, eval(formals()$method))
  nruns <- length(runs)
  if (nruns == 0) stop ("No runs specified")
  reference <- as.integer(reference)
  if (isTRUE(reference < 1L | reference > nruns)) {
    stop("Argument reference does not correspond to a run in runs")
  }
  classes <- sapply(runs, class)
  if (any(classes != "geomcmc")) {
    stop ("Input object is not of class geomcmc")
  }
  Nout <- sapply(runs, "[[", "Nout")
  bfsize1 <- as.double(bfsize1)
  if (length(bfsize1) > nruns) {
    warning ("The number of elements in bfsize1 exceeds the number of runs;
the extra elements will be discarded")
    bfsize1 <- bfsize1[1:nruns]
  }
  if (any(bfsize1 <= 0)) {
    stop ("Argument bfsize1 must be positive")
  } else if (all(bfsize1 <= 1)) {
    Nout1 <- as.integer(bfsize1*Nout)
    if (any(Nout1 == 0)) stop ("Calculated 0 sizes; give a larger bfsize1")
  } else if (all(bfsize1 >= 1)) {
    Nout1 <- as.integer(bfsize1)
    if (any(Nout1 > Nout)) stop ("The bfsize1 exceeds the number of samples")
  } else {
    stop ("Argument bfsize1 is a mix of proportions and sizes")
  }
  Ntot1 <- sum(Nout1)
  Nout2 <- Nout - Nout1
  Ntot2 <- sum(Nout2)
  if (nruns == 1) {
    runs <- runs[[1]]
    out <- list(logbf = 1, logLik1 = runs$logLik[1:Ntot1],
                logLik2 = runs$logLik[-(1:Ntot1)],
                isweights = rep.int(0, Ntot2),
                controlvar = matrix(1, Ntot2, 1), z = runs$z[, -(1:Ntot1)],
                betm0 = runs$betm0, betQ0 = runs$betQ0, ssqdf = runs$ssqdf,
                ssqsc = runs$ssqsc, tsqdf = runs$tsqdf, tsqsc = runs$tsqsc,
                dispersion = runs$dispersion, response = runs$response,
                weights = runs$weights, modelmatrix = runs$modelmatrix,
                locations = runs$locations, family = runs$family,
                referencebf = 0, corrfcn = runs$corrfcn,
                pnts = list(nu = runs$nu, phi = runs$phi, omg = runs$omg,
                  kappa = runs$kappa))
    return(out)
  }
  modelvars <- c("response", "weights", "modelmatrix", "family",
                 "betm0", "betQ0", "ssqdf", "ssqsc",
                 "dispersion", "tsqdf", "tsqsc", "locations", "corrfcn")
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
  ifam <- match(family, c("gaussian", "binomial", "poisson", "Gamma"), 0L)
  betm0 <- model$betm0
  betQ0 <- model$betQ0
  ssqdf <- model$ssqdf
  ssqsc <- model$ssqsc
  dispersion <- model$dispersion
  tsqdf <- model$tsqdf
  tsqsc <- model$tsqsc
  corrfcn <- model$corrfcn
  icf <- match(corrfcn, c("matern", "spherical", "power"))
  loc <- model$locations
  dm <- sp::spDists(loc)
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
    if (family == "binomial" && any(nu_pnts <= 0) &&
        length(unique(nu_pnts)) > 1) {
      stop ("The link functions don't have the same functional form")
    }
  }
  kappa_pnts <- sapply(runs, function(r) r[["kappa"]][1])
  kappa_pnts <- as.double(kappa_pnts)
  if (any(kappa_pnts < 0) & corrfcn %in% c("matern", "power")) {
    stop ("Argument kappa_pnts cannot be negative")
  }
  if (any(kappa_pnts > 2) & corrfcn == "power") {
    stop ("Argument kappa_pnts cannot be more than 2")
  }

  zsample <- lapply(runs, function(r) r[["z"]][r[["whichobs"]], ])
  z1 <- matrix(unlist(mapply(function(z, n) z[, seq_len(n)],
                             zsample, Nout1)), n, Ntot1)
  z2 <- matrix(unlist(mapply(function(z, n) z[, n < seq_len(NCOL(z))],
                             zsample, Nout1)), n, Ntot2)
  logbf <- numeric(nruns)
  lglk1 <- matrix(0., Ntot1, nruns)
  lglk2 <- matrix(0., Ntot2, nruns)
  zcv <- matrix(0., Ntot2, nruns)
  weights <- numeric(Ntot2)
  tsq <- if (ifam == 0) tsqsc else dispersion 
  RUN <- .Fortran("bfsp", weights = weights, zcv = zcv, logbf = logbf,
                  lglk1 = lglk1, lglk2 = lglk2, phi_pnts, omg_pnts,
                  nu_pnts, z1, Nout1, Ntot1, z2, Nout2, Ntot2,
                  y, l, F, dm, betm0, betQ0,
                  ssqdf, ssqsc, max(tsqdf, 0), tsq, kappa_pnts, icf, 
                  n, p, nruns,
                  ifam, imeth)
  refbf <- RUN$logbf[reference]
  logbf <- RUN$logbf - refbf
  weights <- if (Ntot2 > 0) RUN$weights
  lglk2 <- if (Ntot2 > 0) RUN$lglk2
  zcv <- if (Ntot2 > 0) RUN$zcv
  out <- list(logbf = logbf, logLik1 = RUN$lglk1, logLik2 = lglk2,
              isweights = weights, controlvar = zcv, z = z2, betm0 = betm0,
              betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc, tsqdf = tsqdf,
              tsqsc = tsqsc, dispersion = dispersion, response = y, weights = l,
              modelmatrix = F, locations = loc, family = family,
              referencebf = refbf, corrfcn = corrfcn, 
              pnts = list(nu = nu_pnts, phi = phi_pnts, omg = omg_pnts,
                kappa = kappa_pnts))
  out
}
