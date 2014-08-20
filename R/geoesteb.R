##' Empirical Bayes estimation for SGLMM
##'
##' Runs the MCMC sampling, computes the importance weights, and
##' estimates the parameters.
##' @title Empirical Bayes estimation for SGLMM
##' @param formula A representation of the model in the form
##' \code{response ~ terms}. The response must be set to \code{NA}'s
##' at the prediction locations. At the observed locations the
##' response is assumed to be a total of replicated measurements. The
##' number of replications is inputted using the argument
##' \code{weights}.
##' @param family The distribution of the data.
##' @param data An optional data frame containing the variables in the
##' model.
##' @param weights An optional vector of weights. Number of replicated
##' samples for Gaussian and gamma, number of trials for binomial,
##' time length for Poisson.
##' @param subset An optional vector specifying a subset of
##' observations to be used in the fitting process.
##' @param atsample A formula in the form \code{~ x1 + x2 + ... + xd}
##' with the coordinates of the sampled locations.
##' @param parameters A named list with the components "linkp", "phi",
##' "omg", "kappa" with all components being vectors of the same
##' length k. Then, k different MCMC samples will be taken from the
##' models with parameters fixed at those values. For a square grid
##' the function \code{\link[base]{expand.grid}} can be useful.
##' @param estimate A named list with the components "linkp", "phi",
##' "omg", "kappa". Each component must be numeric with length 1, 2,
##' or 3 with elements in increasing order but for the binomial family
##' linkp is also allowed to be the character "logit" and "probit". If
##' its length is 1, then the corresponding parameter is considered to
##' be fixed at that value. If 2, then the two numbers denote the
##' lower and upper bounds for the optimisation of that parameter
##' (infinities are allowed). If 3, these correspond to lower bound,
##' starting value, upper bound for the estimation of that parameter.
##' @param corrfcn Spatial correlation function.
##' @param Nout A scalar or vector of size k. Number of MCMC samples
##' to take for each run of the MCMC algorithm for the estimation of
##' the Bayes factors. See argument
##' \code{parameters}.
##' @param Nthin A scalar or vector of size k. The thinning of the
##' MCMC algorithm for the estimation of the Bayes factors.
##' @param Nbi A scalar or vector of size k. The burn-in of the MCMC
##' algorithm for the estimation of the Bayes factors.
##' @param Npro A scalar. The number of Gibbs samples to use for
##' estimation of the parameters at the last stage and for prediction
##' at the unsampled locations.
##' @param Nprt The thinning of the Gibbs algorithm.
##' @param Nprb The burn-in of the Gibbs algorithm.
##' @param betm0 Prior mean for beta (a vector or scalar).
##' @param betQ0 Prior standardised precision (inverse variance)
##' matrix.
##' @param ssqdf Degrees of freedom for the scaled inverse chi-square
##' prior for the partial sill parameter.
##' @param ssqsc Scale for the scaled inverse chi-square prior for the
##' partial sill parameter.
##' @param zstart Starting value for the MCMC for the GRF. Defaults to
##' 0.
##' @param dispersion The fixed dispersion parameter.
##' @param bfsize1 A scalar or vector of the same length as \code{...}
##' with all integer values or all values in (0, 1]. How many samples
##' (or what proportion of the sample) to use for estimating the Bayes
##' factors at the first stage. The remaining sample will be used for
##' estimating the Bayes factors in the second stage. Setting it to 1
##' will perform only the first stage.
##' @param reference Which model goes in the denominator of the Bayes
##' factors.
##' @param bfmethod Which method to use to calculate the Bayes
##' factors: Reverse logistic or Meng-Wong.
##' @param useCV Whether to use control variates for finer
##' corrections.
##' @param longlat How to compute the distance between locations. If
##' \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##' distance. See \code{\link[sp]{spDists}}.
##' @param control A list of control parameters for the optimisation.
##' See \code{\link[stats]{optim}}.
##' @param verbose Whether to print messages when completing each
##' stage on screen.
##' @return A list with components
##' \itemize{
##' \item \code{parest} The parameter estimates
##' \item \code{skeleton} The skeleton points used augmented with the
##' logarithm of the Bayes factors at those points. 
##' \item \code{optim} The output from the \code{\link[stats]{optim}}
##' function. 
##' \item \code{mcmcsample} The MCMC samples for the remaining
##' parameters and the random field.
##' \item \code{sys_time} The time taken to complete the MCMC
##' sampling, calculation of the importance weights, the
##' optimization and the final MCMC sampling. 
##' }
##' @importFrom sp spDists
##' @export 
ebsglmm <- function (formula,
                     family = c("gaussian", "binomial", "poisson", "Gamma"),
                     data, weights, subset, atsample, parameters, estimate,
                     corrfcn = c("matern", "spherical", "power"), 
                     Nout, Nthin = 1, Nbi = 0, Npro, Nprt = 1, Nprb = 0, 
                     betm0, betQ0, ssqdf, ssqsc,
                     zstart, dispersion = 1,
                     bfsize1 = 0.8, reference = 1, bfmethod = c("RL", "MW"), 
                     useCV = TRUE, longlat = FALSE, 
                     control = list(), verbose = TRUE) {

  ## Family
  family <- match.arg(family)
  ifam <- match(family, eval(formals()$family))

  ## Correlation function
  corrfcn <- match.arg(corrfcn)
  icf <- match(corrfcn, eval(formals()$corrfcn))

  ## Design matrix and data
  if (missing(data)) data <- environment(formula)
  mfc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"),
             names(mfc), 0L)
  mfc <- mfc[c(1L, m)]
  mfc$drop.unused.levels <- TRUE
  mfc$na.action <- "na.pass"
  mfc[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfc, parent.frame())
  mt <- attr(mf, "terms")
  FF <- model.matrix(mt,mf)
  if (!all(is.finite(F))) stop ("Non-finite values in the design matrix")
  p <- NCOL(F)
  yy <- model.response(mf)
  if (!is.vector(yy)) {
    stop ("The response must be a vector")
  }
  ll <- model.weights(mf)

  ## All locations
  locvars <- all.vars(atsample)
  formula1 <- as.formula(paste('~', paste(c(locvars, all.vars(formula)),
                                          collapse = ' + ')))
  mfc1 <- mfc
  mfc1$formula <- formula1
  mf1 <- eval(mfc1, parent.frame())
  m <- match(locvars, names(mf1))
  loc <- as.matrix(mf1[, m])
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations")

  ## Split sample, prediction
  ii <- is.finite(yy)
  y <- yy[ii]
  n <- sum(ii)
  l <- ll[ii]
  l <- if (is.null(l)) rep.int(1.0, n) else as.double(l)
  if (any(!is.finite(l))) stop ("Non-finite values in the weights")
  if (any(l <= 0)) stop ("Non-positive weights not allowed")
  if (family == "binomial") {
    l <- l - y # Number of failures
  }
  F <- FF[ii, ]
  dm <- sp::spDists(loc[ii, ], longlat = longlat)
  n0 <- sum(!ii)
  if (n0 > 0) {
    F0 <- FF[!ii, ]
    dmdm0 <- sp::spDists(loc[ii, ], loc[!ii, ], longlat = longlat)
  } else {
    F0 <- numeric(p)
    dmdm0 <- numeric(n)
  }

  ## Priors
  if (all(is.finite(betQ0[upper.tri(betQ0)]))) {
    if (length(betQ0) == 1 && betQ0[1] == 0) {
      ## Uniform prior
      betQ0 <- matrix(0, p, p)
      betm0 <- rep(0, p)
    } else if (length(betQ0) == p*p) {
      betQ0 <- matrix(as.double(betQ0), p, p)
      betQ0[lower.tri(betQ0)] <- 0
      betQ0eig <- eigen(t(betQ0), 1, 1)$values
      if (any (betQ0eig < sqrt(.Machine$double.eps))) {
        stop ('betQ0 not > 0 within tolerance')
      } else {
        betm0 <- rep_len(as.double(betm0), p)
        modeldf <- as.double(n + ssqdf)
      }
    } else stop ('Bad betQ0')
  } else stop ('Non-finite betQ0')
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## Read and check parameters
  if(!is.list(parameters)) stop ("Argument parameters must be a list")
  parnm <- c("linkp", "phi", "omg", "kappa")
  if (!all(parnm %in% names(parameters))) {
    stop (paste("Argument parameters must have the names",
                paste(parnm, collapse=" ")))
  } else {
    parameters <- parameters[parnm]
  }
  phi <- parameters[["phi"]]
  nruns <- length(phi)
  omg <- parameters[["omg"]]
  kappa <- parameters[["kappa"]]
  linkp <- parameters[["linkp"]]
  if (length(omg) != nruns | length(kappa) != nruns | length(linkp) != nruns) {
    stop ("Elements in parameters don't have the same length")
  }
  if (any (phi < 0)) stop ("Element phi in parameters must be non-negative")
  if (any (omg < 0)) stop ("Element omg in parameters must be non-negative")

  kappa <- as.double(kappa)
  if (any(kappa < 0) & corrfcn %in% c("matern", "power")) {
    stop ("Argument kappa cannot be negative")
  }
  if (any(kappa > 2) & corrfcn == "power") {
    stop ("Argument kappa cannot be more than 2")
  }
  if (corrfcn == "spherical" & NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }
  
  if (is.character(linkp)) {
    if (family == "binomial") {
      if (all(linkp == "logit")) {
        nu <- rep.int(-1, nruns)
      } else if (all(linkp == "probit")) {
        nu <- rep.int(0, nruns)
      } else stop ("Cannot recognise character link for binomial")
    } else stop ("Character link is only allowed for binomial")
  } else if (!is.numeric(linkp)) {
    stop ("Element linkp in parameters must be numeric, or in the case of
the binomial can also be the character \"logit\" or \"probit\"")
  } else {
    nu <- as.double(linkp)
    if (family == "binomial" && any(nu <= 0)) {
      stop ("The robit link parameter must be positive")
    }
  }

  ## Method for computing the Bayes factors
  bfmethod <- match.arg(bfmethod)
  imeth <- match(bfmethod, eval(formals()$bfmethod))

  ## Read and check estimate argument
  if (!is.list(estimate) | length(estimate) != 4) {
    stop ("Argument estimate must be a list with 4 components")
  }
  if (!all(parnm %in% names(estimate))) {
    stop (paste("Named components in the argument estimate",
                "must be", paste(parnm, collapse = " ")))
  } else {
    estimate <- estimate[parnm]
  }
  lower <- upper <- pstart <- rep.int(NA, 4)
  estim <- rep.int(FALSE, 4)
  linkpe <- estimate[["linkp"]]
  if (is.character(linkpe)) {
    if (family == "binomial") {
      if (linkpe == "logit") {
        if (linkp == "logit") {
          pstart[1] <- -1
          estim[1] <- FALSE
        } else {
          stop ("Link in estimate inconsistent with link in parameters")
        }
      } else if (linkpe == "probit") {
        if (linkp == "probit") {
          pstart[1] <- 0
          estim[1] <- FALSE
        } else {
          stop ("Link in estimate inconsistent with link in parameters")
        }
      } else {
        stop ("Character link for the binomial family can be either
\"logit\" or \"probit\" in argument estimate")
      }
    } else stop ("Character link in argument estimate is only used for
the binomial family")
  } else if (is.numeric(linkpe)) {
    if (length(linkpe) < 1 | length(linkpe) > 3) {
      stop ("The length for a numeric component in estimate must be 1, 2, or 3")
    } else if (length(linkpe) == 1) {
      pstart[1] <- linkpe
      estim[1] <- FALSE
    } else if (length(linkpe) == 2) {
      pstart[1] <- NA
      estim[1] <- TRUE
      lower[1] <- linkpe[1]
      upper[1] <- linkpe[2]
      if (lower[1] >= upper[1]) {
        stop ("The lower bound must be less than the upper bound for linkp
in estimate")
      }
    } else {
      pstart[1] <- linkpe[2]
      estim[1] <- TRUE
      lower[1] <- linkpe[1]
      upper[1] <- linkpe[3]
      if (lower[1] > estim[1] | estim[1] > upper[1] | lower[1] == upper[1]) {
        stop ("The elements in the component linkp in estimate must be ordered")
      }
    }
  } else {
    stop ("The element linkp in estimate must be either numeric or character")
  }
  for (i in 2:4) {
    ppp <- estimate[[parnm[i]]]
    if (!is.numeric(ppp)) {
      stop(paste("The element", parnm[i], "in estimate must be numeric"))
    }
    lppp <- length(ppp)
    if (lppp < 1 | lppp > 3) {
      stop (paste("The element", parnm[i], "in estimate must have 1, 2, or 3
components"))
    }
    if (lppp == 1) {
      pstart[i] <- ppp
      estim[i] <- FALSE
    } else if (lppp == 2) {
      pstart[i] <- NA
      estim[i] <- TRUE
      lower[i] <- ppp[1]
      upper[i] <- ppp[2]
      if (lower[i] >= upper[i]) {
        stop (paste("The lower bound must be less than the upper bound for",
                     parnm[i], "in estimate"))
      }
    } else {
      pstart[i] <- ppp[2]
      estim[i] <- TRUE
      lower[i] <- ppp[1]
      upper[i] <- ppp[3]
      if (lower[i] > estim[i] | estim[i] > upper[i] | lower[i] == upper[i]) {
        stop (paste("The elements in the component", parnm[i],
                     "in estimate must be ordered"))
      }
    }
  }

  ## MCMC samples
  Nout <- rep_len(as.integer(Nout), nruns)
  Nbi <- rep_len(as.integer(Nbi), nruns)
  Nthin <- rep_len(as.integer(Nthin), nruns)
  Ntot <- sum(Nout)
  z <- matrix(0, n, Ntot)
  lglk <- numeric(Ntot)
  if (!missing(zstart)) z[, 1 + cumsum(Nout[-nruns])] <- zstart
  tsqsc <- tsqdf <- 0

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

  tsq <- if (ifam > 0) dispersion else tsqsc

  ## RUN MCMC
  tm1 <- system.time({
    RUN1 <- .Fortran("samplemulti", lglk = lglk, z = z, phi, omg, y, l,
                     F, betm0, betQ0, ssqdf, ssqsc, kappa, icf, nu, tsqdf, tsq,
                     dm, Ntot, Nout, Nbi, Nthin, n, p, nruns, ifam)
  })

  if (isTRUE(verbose)) {
    message ("Completed MCMC sampling: ", round(tm1[1]), " sec")
  }

  ## Prepare data for first stage
  z <- RUN1$z
  lz <- unlist(lapply(1:nruns, function(i)
                      c(rep.int(TRUE, Nout1[i]), rep.int(FALSE, Nout2[i]))))
  z1 <- z[, lz]
  z2 <- z[, !lz]

  logbf <- numeric(nruns)
  lglk1 <- matrix(0., Ntot1, nruns)
  lglk2 <- matrix(0., Ntot2, nruns)
  isweights <- numeric(Ntot2)
  zcv <- matrix(0., Ntot2, nruns)

  ## RUN estimation of Bayes factors
  tm2 <- system.time({
    RUN2 <- .Fortran("bfsp", isweights = isweights, zcv = zcv, logbf = logbf,
                     lglk1 = lglk1, lglk2 = lglk2, phi, omg,
                     nu, z1, Nout1, Ntot1, z2, Nout2, Ntot2,
                     y, l, F, dm, betm0, betQ0,
                     ssqdf, ssqsc, tsqdf, tsq, kappa, icf, n, p, nruns,
                     ifam, imeth)
  })

  if (isTRUE(verbose)) {
    message ("Calculated Bayes factors and importance weights at skeleton
points: ", round(tm2[1]), " sec")
  }

  refbf <- RUN2$logbf[reference]
  logbf <- RUN2$logbf - refbf
  isweights <- RUN2$isweights
  zcv <- RUN2$zcv

  ## Optimisation
  imaxlogbf <- which.max(logbf)
  wbf <- exp(logbf - logbf[imaxlogbf] - log(sum(exp(logbf - logbf[imaxlogbf]))))
  pstart.d <- colSums(cbind(nu, phi, omg, kappa)*wbf)
  i <- is.na(pstart) & estim
  pstart[i] <- pmax(pmin(upper[i], pstart.d[i]), lower[i])
  ## Function to optimise
  fn <- if (useCV) {
    function (par) {
      parin <- pstart
      parin[estim] <- par
      RUN <- .Fortran('calcb_cv', 0.0, parin[2], parin[1], parin[3], parin[4],
                      icf, 1L, 1L, Ntot2, z2, isweights, zcv, n, p, nruns, betm0,
                      betQ0, ssqdf, ssqsc, tsqdf, tsq, y, l, F, dm,
                      ifam)
      -RUN[[1]][1]
    }
  } else {
    function (par) {
      parin <- pstart
      parin[estim] <- par
      RUN <- .Fortran('calcb_st', 0.0, parin[2], parin[1], parin[3], parin[4],
                      icf, 1L, 1L, Ntot2, z2, isweights, n, p, betm0,
                      betQ0, ssqdf, ssqsc, tsqdf, tsq, y, l, F, dm,
                      ifam)
      -RUN[[1]][1]
    }
  }
  method <- if (sum(estim) == 1) "Brent" else "L-BFGS-B"
  tm3 <- system.time({
    op <- stats::optim(pstart[estim], fn, method = method,
                       lower = lower[estim], upper = upper[estim],
                       control = control)
  })

  if (isTRUE(verbose)) {
    message ("Finished optimization: ", round(tm3[1]), " sec")
  }

  parest <- pstart
  parest[estim] <- op$par
  names(parest) <- parnm

  ## Perform prediction
  Npro <- as.integer(Npro)
  if (Npro > 0) {
    Nprt <- as.integer(Nprt)
    Nprb <- as.integer(Nprb)
    lglks <- numeric(Npro)
    zs <- matrix(0, n, Npro)
    z0s <- matrix(0, n0, Npro)
    beta <- matrix(0, p, Npro)
    ssq <- numeric(Npro)
    acc <- 0L
    phis <- rep.int(parest[2], Npro)
    omgs <- rep.int(parest[3], Npro)
    tmppars <- rep.int(0, 4)
    tm4 <- system.time({
      RUN4 <- .Fortran("mcspsample", ll = lglks, z = zs, z0 = z0s,
                       beta = beta, ssq = ssq, phis, omgs, acc = acc,
                       y, l, F, F0, betm0, betQ0, ssqdf, ssqsc,
                       tmppars, 0, tmppars, 0, parest[4], icf, parest[1], 
                       tsq, dm, dmdm0, Npro, Nprb, Nprt, n, n0, p, ifam)
    })
    if (isTRUE(verbose)) {
      message ("Performed Gibbs sampling: ", round(tm4[1]), " sec")
    }
    ll <- RUN4$ll
    zz0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN4$z
    zz0[!ii, ] <- RUN4$z0
    beta <- RUN4$beta
    ssq <- RUN4$ssq
    acc_ratio <- RUN4$acc/Npro
    sample <- list(z = zz0, beta = beta, ssq = ssq, acc_ratio = acc_ratio,
                   whichobs = ii)
  } else {
    sample <- NULL
    tm4 <- NULL
  }
  
  times <- rbind(sampling = tm1, importance = tm2, optimization = tm3,
                 MCMC = tm4)
  out <- list(parest = parest, skeleton = cbind(parameters, logbf = logbf),
              optim = op, mcmcsample = sample, sys_time = times)
  out
}


##' Empirical Bayes estimation for the spatial transformed Gaussian model.
##'
##' Runs the MCMC sampling, computes the importance weights, and
##' estimates the parameters.
##' @title Empirical Bayes estimation for the TGRFM
##' @param formula A representation of the model in the form
##' \code{response ~ terms}. The response must be set to \code{NA}'s
##' at the prediction locations (the function \code{\link{stackdata}}
##' can be used to create such data frame). At the observed locations the
##' response is assumed to be a total of replicated measurements. The
##' number of replications is inputted using the argument
##' \code{weights}.
##' @param data An optional data frame containing the variables in the
##' model.
##' @param weights An optional vector of weights. Number of replicated
##' samples for Gaussian and gamma, number of trials for binomial,
##' time length for Poisson.
##' @param subset An optional vector specifying a subset of
##' observations to be used in the fitting process.
##' @param atsample A formula in the form \code{~ x1 + x2 + ... + xd}
##' with the coordinates of the sampled locations.
##' @param parameters A named list with the components "linkp", "phi",
##' "omg", "kappa" with all components being vectors of the same
##' length k. Then, k different MCMC samples will be taken from the
##' models with parameters fixed at those values. For a square grid
##' the function \code{\link[base]{expand.grid}} can be useful.
##' @param estimate A named list with the components "linkp", "phi",
##' "omg", "kappa". Each component must be numeric with length 1, 2,
##' or 3 with elements in increasing order but for the binomial family
##' linkp is also allowed to be the character "logit" and "probit". If
##' its length is 1, then the corresponding parameter is considered to
##' be fixed at that value. If 2, then the two numbers denote the
##' lower and upper bounds for the optimisation of that parameter
##' (infinities are allowed). If 3, these correspond to lower bound,
##' starting value, upper bound for the estimation of that parameter.
##' @param corrfcn Spatial correlation function.
##' @param Nout A scalar or vector of size k. Number of MCMC samples
##' to take for each run of the MCMC algorithm for the estimation of
##' the Bayes factors. See argument
##' \code{parameters}.
##' @param Nthin A scalar or vector of size k. The thinning of the
##' MCMC algorithm for the estimation of the Bayes factors.
##' @param Nbi A scalar or vector of size k. The burn-in of the MCMC
##' algorithm for the estimation of the Bayes factors.
##' @param Npro A scalar. The number of Gibbs samples to use for
##' estimation of the parameters at the last stage and for prediction
##' at the unsampled locations.
##' @param Nprt The thinning of the Gibbs algorithm.
##' @param Nprb The burn-in of the Gibbs algorithm.
##' @param betm0 Prior mean for beta (a vector or scalar).
##' @param betQ0 Prior standardised precision (inverse variance)
##' matrix.
##' @param ssqdf Degrees of freedom for the scaled inverse chi-square
##' prior for the partial sill parameter.
##' @param ssqsc Scale for the scaled inverse chi-square prior for the
##' partial sill parameter.
##' @param tsqdf Degrees of freedom for the scaled inverse chi-square
##' prior for the measurement error parameter.
##' @param tsqsc Scale for the scaled inverse chi-square prior for the
##' measurement error parameter.
##' @param zstart Starting value for the MCMC for the GRF. Defaults to
##' 0.
##' @param dispersion The fixed dispersion parameter.
##' @param bfsize1 A scalar or vector of the same length as \code{...}
##' with all integer values or all values in (0, 1]. How many samples
##' (or what proportion of the sample) to use for estimating the Bayes
##' factors at the first stage. The remaining sample will be used for
##' estimating the Bayes factors in the second stage. Setting it to 1
##' will perform only the first stage.
##' @param reference Which model goes in the denominator of the Bayes
##' factors.
##' @param bfmethod Which method to use to calculate the Bayes
##' factors: Reverse logistic or Meng-Wong.
##' @param useCV Whether to use control variates for finer
##' corrections.
##' @param longlat How to compute the distance between locations. If
##' \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##' distance. See \code{\link[sp]{spDists}}.
##' @param control A list of control parameters for the optimisation.
##' See \code{\link[stats]{optim}}.
##' @param verbose Whether to print messages when completing each
##' stage on screen.
##' @return A list with components
##' \itemize{
##' \item \code{parest} The parameter estimates
##' \item \code{skeleton} The skeleton points used augmented with the
##' logarithm of the Bayes factors at those points. 
##' \item \code{optim} The output from the \code{\link[stats]{optim}}
##' function. 
##' \item \code{mcmcsample} The MCMC samples for the remaining
##' parameters and the random field.
##' \item \code{sys_time} The time taken to complete the MCMC
##' sampling, calculation of the importance weights, the
##' optimization and the final MCMC sampling.
##' }
##' @importFrom sp spDists
##' @export 
ebstrga <- function (formula,
                     data, weights, subset, atsample, parameters, estimate,
                     corrfcn = c("matern", "spherical", "power"), 
                     Nout, Nthin = 1, Nbi = 0, Npro, Nprt = 1, Nprb = 0, 
                     betm0, betQ0, ssqdf, ssqsc,
                     tsqdf, tsqsc, zstart, dispersion = 1,
                     bfsize1 = 0.8, reference = 1, bfmethod = c("RL", "MW"), 
                     useCV = TRUE, longlat = FALSE, 
                     control = list(), verbose = TRUE) {

  ## Family
  family <- "transformed-gaussian"
  ifam <- 0L

  ## Correlation function
  corrfcn <- match.arg(corrfcn)
  icf <- match(corrfcn, eval(formals()$corrfcn))

  ## Design matrix and data
  if (missing(data)) data <- environment(formula)
  mfc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"),
             names(mfc), 0L)
  mfc <- mfc[c(1L, m)]
  mfc$drop.unused.levels <- TRUE
  mfc$na.action <- "na.pass"
  mfc[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfc, parent.frame())
  mt <- attr(mf, "terms")
  FF <- model.matrix(mt,mf)
  if (!all(is.finite(F))) stop ("Non-finite values in the design matrix")
  p <- NCOL(F)
  yy <- model.response(mf)
  if (!is.vector(yy)) {
    stop ("The response must be a vector")
  }
  ll <- model.weights(mf)

  ## All locations
  locvars <- all.vars(atsample)
  formula1 <- as.formula(paste('~', paste(c(locvars, all.vars(formula)),
                                          collapse = ' + ')))
  mfc1 <- mfc
  mfc1$formula <- formula1
  mf1 <- eval(mfc1, parent.frame())
  m <- match(locvars, names(mf1))
  loc <- as.matrix(mf1[, m])
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations")

  ## Split sample, prediction
  ii <- is.finite(yy)
  y <- yy[ii]
  n <- sum(ii)
  l <- ll[ii]
  l <- if (is.null(l)) rep.int(1.0, n) else as.double(l)
  if (any(!is.finite(l))) stop ("Non-finite values in the weights")
  if (any(l <= 0)) stop ("Non-positive weights not allowed")
  if (family == "binomial") {
    l <- l - y # Number of failures
  }
  F <- FF[ii, ]
  dm <- sp::spDists(loc[ii, ], longlat = longlat)
  n0 <- sum(!ii)
  if (n0 > 0) {
    F0 <- FF[!ii, ]
    dmdm0 <- sp::spDists(loc[ii, ], loc[!ii, ], longlat = longlat)
  } else {
    F0 <- numeric(p)
    dmdm0 <- numeric(n)
  }

  ## Priors
  if (all(is.finite(betQ0[upper.tri(betQ0)]))) {
    if (length(betQ0) == 1 && betQ0[1] == 0) {
      ## Uniform prior
      betQ0 <- matrix(0, p, p)
      betm0 <- rep(0, p)
    } else if (length(betQ0) == p*p) {
      betQ0 <- matrix(as.double(betQ0), p, p)
      betQ0[lower.tri(betQ0)] <- 0
      betQ0eig <- eigen(t(betQ0), 1, 1)$values
      if (any (betQ0eig < sqrt(.Machine$double.eps))) {
        stop ('betQ0 not > 0 within tolerance')
      } else {
        betm0 <- rep_len(as.double(betm0), p)
        modeldf <- as.double(n + ssqdf)
      }
    } else stop ('Bad betQ0')
  } else stop ('Non-finite betQ0')
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## Read and check parameters
  if(!is.list(parameters)) stop ("Argument parameters must be a list")
  parnm <- c("linkp", "phi", "omg", "kappa")
  if (!all(parnm %in% names(parameters))) {
    stop (paste("Argument parameters must have the names",
                paste(parnm, collapse=" ")))
  } else {
    parameters <- parameters[parnm]
  }
  phi <- parameters[["phi"]]
  nruns <- length(phi)
  omg <- parameters[["omg"]]
  kappa <- parameters[["kappa"]]
  linkp <- parameters[["linkp"]]
  if (length(omg) != nruns | length(kappa) != nruns | length(linkp) != nruns) {
    stop ("Elements in parameters don't have the same length")
  }
  if (any (phi < 0)) stop ("Element phi in parameters must be non-negative")
  if (any (omg < 0)) stop ("Element omg in parameters must be non-negative")
  
  kappa <- as.double(kappa)
  if (any(kappa < 0) & corrfcn %in% c("matern", "power")) {
    stop ("Argument kappa cannot be negative")
  }
  if (any(kappa > 2) & corrfcn == "power") {
    stop ("Argument kappa cannot be more than 2")
  }
  if (corrfcn == "spherical" & NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }

  nu <- as.double(linkp)

  ## Method for computing the Bayes factors
  bfmethod <- match.arg(bfmethod)
  imeth <- match(bfmethod, eval(formals()$bfmethod))

  ## Read and check estimate argument
  if (!is.list(estimate) | length(estimate) != 4) {
    stop ("Argument estimate must be a list with 4 components")
  }
  if (!all(parnm %in% names(estimate))) {
    stop (paste("Named components in the argument estimate",
                "must be", paste(parnm, collapse = " ")))
  } else {
    estimate <- estimate[parnm]
  }
  lower <- upper <- pstart <- rep.int(NA, 4)
  estim <- rep.int(FALSE, 4)
  linkpe <- estimate[["linkp"]]
  if (is.character(linkpe)) {
    if (family == "binomial") {
      if (linkpe == "logit") {
        if (linkp == "logit") {
          pstart[1] <- -1
          estim[1] <- FALSE
        } else {
          stop ("Link in estimate inconsistent with link in parameters")
        }
      } else if (linkpe == "probit") {
        if (linkp == "probit") {
          pstart[1] <- 0
          estim[1] <- FALSE
        } else {
          stop ("Link in estimate inconsistent with link in parameters")
        }
      } else {
        stop ("Character link for the binomial family can be either
\"logit\" or \"probit\" in argument estimate")
      }
    } else stop ("Character link in argument estimate is only used for
the binomial family")
  } else if (is.numeric(linkpe)) {
    if (length(linkpe) < 1 | length(linkpe) > 3) {
      stop ("The length for a numeric component in estimate must be 1, 2, or 3")
    } else if (length(linkpe) == 1) {
      pstart[1] <- linkpe
      estim[1] <- FALSE
    } else if (length(linkpe) == 2) {
      pstart[1] <- NA
      estim[1] <- TRUE
      lower[1] <- linkpe[1]
      upper[1] <- linkpe[2]
      if (lower[1] >= upper[1]) {
        stop ("The lower bound must be less than the upper bound for linkp
in estimate")
      }
    } else {
      pstart[1] <- linkpe[2]
      estim[1] <- TRUE
      lower[1] <- linkpe[1]
      upper[1] <- linkpe[3]
      if (lower[1] > estim[1] | estim[1] > upper[1] | lower[1] == upper[1]) {
        stop ("The elements in the component linkp in estimate must be ordered")
      }
    }
  } else {
    stop ("The element linkp in estimate must be either numeric or character")
  }
  for (i in 2:4) {
    ppp <- estimate[[parnm[i]]]
    if (!is.numeric(ppp)) {
      stop(paste("The element", parnm[i], "in estimate must be numeric"))
    }
    lppp <- length(ppp)
    if (lppp < 1 | lppp > 3) {
      stop (paste("The element", parnm[i], "in estimate must have 1, 2, or 3
components"))
    }
    if (lppp == 1) {
      pstart[i] <- ppp
      estim[i] <- FALSE
    } else if (lppp == 2) {
      pstart[i] <- NA
      estim[i] <- TRUE
      lower[i] <- ppp[1]
      upper[i] <- ppp[2]
      if (lower[i] >= upper[i]) {
        stop (paste("The lower bound must be less than the upper bound for",
                     parnm[i], "in estimate"))
      }
    } else {
      pstart[i] <- ppp[2]
      estim[i] <- TRUE
      lower[i] <- ppp[1]
      upper[i] <- ppp[3]
      if (lower[i] > estim[i] | estim[i] > upper[i] | lower[i] == upper[i]) {
        stop (paste("The elements in the component", parnm[i],
                     "in estimate must be ordered"))
      }
    }
  }

  ## MCMC samples
  Nout <- rep_len(as.integer(Nout), nruns)
  Nbi <- rep_len(as.integer(Nbi), nruns)
  Nthin <- rep_len(as.integer(Nthin), nruns)
  Ntot <- sum(Nout)
  z <- matrix(0, n, Ntot)
  lglk <- numeric(Ntot)
  if (!missing(zstart)) z[, 1 + cumsum(Nout[-nruns])] <- zstart
  dispersion <- NA

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

  tsq <- if (ifam > 0) dispersion else tsqsc

  ## RUN MCMC
  tm1 <- system.time({
    RUN1 <- .Fortran("samplemulti", lglk = lglk, z = z, phi, omg, y, l,
                     F, betm0, betQ0, ssqdf, ssqsc, kappa, icf, nu, tsqdf, tsq,
                     dm, Ntot, Nout, Nbi, Nthin, n, p, nruns, ifam)
  })

  if (isTRUE(verbose)) {
    message ("Completed MCMC sampling: ", round(tm1[1]), " sec")
  }

  ## Prepare data for first stage
  z <- RUN1$z
  lz <- unlist(lapply(1:nruns, function(i)
                      c(rep.int(TRUE, Nout1[i]), rep.int(FALSE, Nout2[i]))))
  z1 <- z[, lz]
  z2 <- z[, !lz]

  logbf <- numeric(nruns)
  lglk1 <- matrix(0., Ntot1, nruns)
  lglk2 <- matrix(0., Ntot2, nruns)
  isweights <- numeric(Ntot2)
  zcv <- matrix(0., Ntot2, nruns)

  ## RUN estimation of Bayes factors
  tm2 <- system.time({
    RUN2 <- .Fortran("bfsp", isweights = isweights, zcv = zcv, logbf = logbf,
                     lglk1 = lglk1, lglk2 = lglk2, phi, omg,
                     nu, z1, Nout1, Ntot1, z2, Nout2, Ntot2,
                     y, l, F, dm, betm0, betQ0,
                     ssqdf, ssqsc, tsqdf, tsq, kappa, icf, n, p, nruns,
                     ifam, imeth)
  })

  if (isTRUE(verbose)) {
    message ("Calculated Bayes factors and importance weights at skeleton
points: ", round(tm2[1]), " sec")
  }

  refbf <- RUN2$logbf[reference]
  logbf <- RUN2$logbf - refbf
  isweights <- RUN2$isweights
  zcv <- RUN2$zcv

  ## Optimisation
  imaxlogbf <- which.max(logbf)
  pmaxlogbf <- c(nu[imaxlogbf], phi[imaxlogbf], omg[imaxlogbf],
                 kappa[imaxlogbf])
  i <- is.na(pstart) & estim
  pstart[i] <- pmax(pmin(upper[i], pmaxlogbf[i]), lower[i])
  ## Function to optimise
  fn <- if (useCV) {
    function (par) {
      parin <- pstart
      parin[estim] <- par
      RUN <- .Fortran('calcb_cv', 0.0, parin[2], parin[1], parin[3], parin[4],
                      icf, 1L, 1L, Ntot2, z2, isweights, zcv, n, p, nruns, betm0,
                      betQ0, ssqdf, ssqsc, tsqdf, tsq, y, l, F, dm,
                      ifam)
      -RUN[[1]][1]
    }
  } else {
    function (par) {
      parin <- pstart
      parin[estim] <- par
      RUN <- .Fortran('calcb_st', 0.0, parin[2], parin[1], parin[3], parin[4],
                      icf, 1L, 1L, Ntot2, z2, isweights, n, p, betm0,
                      betQ0, ssqdf, ssqsc, tsqdf, tsq, y, l, F, dm,
                      ifam)
      -RUN[[1]][1]
    }
  }
  method <- if (sum(estim) == 1) "Brent" else "L-BFGS-B"
  tm3 <- system.time({
    op <- stats::optim(pstart[estim], fn, method = method,
                       lower = lower[estim], upper = upper[estim],
                       control = control)
  })

  if (isTRUE(verbose)) {
    message ("Finished optimization: ", round(tm3[1]), " sec")
  }

  ## Output
  parest <- pstart
  parest[estim] <- op$par
  names(parest) <- parnm

  ## Perform prediction
  Npro <- as.integer(Npro)
  if (Npro > 0) {
    Nprt <- as.integer(Nprt)
    Nprb <- as.integer(Nprb)
    lglks <- numeric(Npro)
    zs <- matrix(0, n, Npro)
    z0s <- matrix(0, n0, Npro)
    beta <- matrix(0, p, Npro)
    ssq <- numeric(Npro)
    tsqs <- numeric(Npro)
    acc <- 0L
    phis <- rep.int(parest[2], Npro)
    omgs <- rep.int(parest[3], Npro)
    tmppars <- rep.int(0, 4)
    tm4 <- system.time({
      RUN4 <- .Fortran("trgasample", ll = lglks, z = zs, z0 = z0s,
                       beta = beta, ssq = ssq, tsq = tsqs,
                       phis, omgs, acc = acc,
                       y, l, F, F0, betm0, betQ0, ssqdf, ssqsc,
                       tsqdf, tsqsc, 
                       tmppars, 0, tmppars, 0, parest[4], icf, parest[1],
                       dm, dmdm0, Npro, Nprb, Nprt, n, n0, p, ifam)
    })
    if (isTRUE(verbose)) {
      message ("Performed Gibbs sampling: ", round(tm4[1]), " sec")
    }
    ll <- RUN4$ll
    zz0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN4$z
    zz0[!ii, ] <- RUN4$z0
    beta <- RUN4$beta
    ssq <- RUN4$ssq
    acc_ratio <- RUN4$acc/Npro
    sample <- list(z = zz0, beta = beta, ssq = ssq, acc_ratio = acc_ratio,
                   whichobs = ii)
  } else {
    sample <- NULL
    tm4 <- NULL
  }
  
  times <- rbind(sampling = tm1, importance = tm2, optimization = tm3,
                 MCMC = tm4)
  out <- list(parest = parest, skeleton = cbind(parameters, logbf = logbf),
              optim = op, mcmcsample = sample, sys_time = times)
  out
}
