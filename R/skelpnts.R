### skelpnts.R ---
##
## Created: Wed, 1 Jun, 2016 14:05 (BST)
## Last-Updated: Sat, 4 Jun, 2016 15:29 (BST)
##
######################################################################
##
### Commentary: This functions are used for deriving appropriate
### skeleton points.
##
######################################################################


##' Log-likelihood approximation.
##'
##' Computes and approximation to the log-likelihood for the given
##' parameters using integrated nested Laplace approximations.
##' @title Log-likelihood approximation
##' @param par_vals A data frame with the components "linkp", "phi",
##'   "omg", "kappa". The approximation will be computed at each row
##'   of the data frame.
##' @param formula A representation of the model in the form
##'   \code{response ~ terms}.
##' @param family The distribution of the response. Can be one of the
##'   options in \code{\link{.geoBayes_models}} or
##'   \code{"transformed.gaussian"}.
##' @param data An optional data frame containing the variables in the
##'   model.
##' @param weights An optional vector of weights. Number of replicated
##'   samples for Gaussian and gamma, number of trials for binomial,
##'   time length for Poisson.
##' @param subset An optional vector specifying a subset of
##'   observations to be used in the fitting process.
##' @param offset See \code{\link[stats]{lm}}.
##' @param atsample A formula in the form \code{~ x1 + x2 + ... + xd}
##'   with the coordinates of the sampled locations.
##' @param corrfcn Spatial correlation function. Can be one of the
##'   choices in \code{\link{.geoBayes_corrfcn}}.
##' @param np The number of integration points for the spatial
##'   variance parameter sigma^2. The total number of points will be
##'   \code{2*np + 1}.
##' @param betm0 Prior mean for beta (a vector or scalar).
##' @param betQ0 Prior standardised precision (inverse variance)
##'   matrix. Can be a scalar, vector or matrix. The first two imply a
##'   diagonal with those elements. Set this to 0 to indicate a flat
##'   improper prior.
##' @param ssqdf Degrees of freedom for the scaled inverse chi-square
##'   prior for the partial sill parameter.
##' @param ssqsc Scale for the scaled inverse chi-square prior for the
##'   partial sill parameter.
##' @param tsqdf Degrees of freedom for the scaled inverse chi-square
##'   prior for the measurement error parameter.
##' @param tsqsc Scale for the scaled inverse chi-square prior for the
##'   measurement error parameter.
##' @param dispersion The fixed dispersion parameter.
##' @param longlat How to compute the distance between locations. If
##'   \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##'   distance. See \code{\link[sp]{spDists}}.
##' @return A list with components
##' \itemize{
##' \item \code{par_vals} A data frame of the parameter values.
##' \item \code{aloglik} The approximate log-likelihood at thos
##'   parameter values.
##' }
##' @useDynLib geoBayes llikparscalc
##' @importFrom stats model.offset
##' @export
##' @references Evangelou, E., & Roy, V. (2019). Estimation and prediction for
##'   spatial generalized linear mixed models with parametric links
##'   via reparameterized importance sampling. Spatial Statistics, 29,
##'   289-315.
alik_inla <- function (par_vals, formula,
                       family = "gaussian",
                       data, weights, subset, offset, atsample,
                       corrfcn = "matern",
                       np, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsqsc,
                       dispersion = 1, longlat = FALSE)
{
  ## Family
  ifam <- .geoBayes_family(family)
  if (ifam) {
    family <- .geoBayes_models$family[ifam]
  } else {
    stop ("This family has not been implemented.")
  }

  ## Logical input
  longlat <- as.logical(longlat)

  ## Correlation function and parameters
  icf <- .geoBayes_correlation(corrfcn)
  corrfcn <- .geoBayes_corrfcn$corrfcn[icf]
  if (!.geoBayes_corrfcn$needkappa[icf]) par_vals$kappa <- 0
  kappa <- .geoBayes_getkappa(par_vals$kappa, icf)
  phi <- as.double(par_vals$phi)
  if (any(phi < 0)) stop ("Input phi must >= 0")
  omg <- as.double(par_vals$omg)
  if (any(omg < 0)) stop ("Input omg must >= 0")

  ## Check the link parameter
  if (!.geoBayes_models$needlinkp[ifam]) par_vals$linkp <- 0
  nu <- .geoBayes_getlinkp(par_vals$linkp, ifam)

  ## All parameters
  allpars <- data.frame(nu = nu, phi = phi, omg = omg, kappa = kappa)
  nu <- as.double(allpars[["nu"]])
  phi <- as.double(allpars[["phi"]])
  omg <- as.double(allpars[["omg"]])
  kappa <- as.double(allpars[["kappa"]])
  npars <- nrow(allpars)

  ## Check dispersion input
  tsq <- as.double(dispersion)
  if (tsq <= 0) stop ("Invalid argument dispersion")
  tsqdf <- 0

  ## Design matrix and data
  if (missing(data)) data <- environment(formula)
  if (length(formula) != 3) stop ("The formula input is incomplete.")
  if ("|" == all.names(formula[[2]], TRUE, 1)) formula[[2]] <- formula[[2]][[2]]
  mfc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "offset"),
             names(mfc), 0L)
  mfc <- mfc[c(1L, m)]
  mfc$formula <- formula
  mfc$drop.unused.levels <- TRUE
  mfc$na.action <- "na.pass"
  mfc[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfc, parent.frame())
  mt <- attr(mf, "terms")
  FF <- model.matrix(mt,mf)
  if (!all(is.finite(FF))) stop ("Non-finite values in the design matrix")
  p <- NCOL(FF)
  yy <- unclass(model.response(mf))
  if (!is.vector(yy)) {
    stop ("The response must be a vector")
  }
  yy <- as.double(yy)
  ll <- as.double(model.weights(mf))
  oofset <- as.vector(model.offset(mf))
  if (!is.null(oofset)) {
    if (length(oofset) != NROW(yy)) {
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                    length(oofset), NROW(yy)), domain = NA)
    } else {
      oofset <- as.double(oofset)
    }
  } else {
    oofset <- double(NROW(yy))
  }

  ## All locations
  atsample <- update(atsample, NULL ~ .)
  locvars <- all.vars(atsample)
  formula1 <- as.formula(paste('~', paste(c(locvars, all.vars(formula)),
                                          collapse = ' + ')))
  mfc1 <- mfc
  mfc1$formula <- formula1
  mf1 <- eval(mfc1, parent.frame())
  m <- match(locvars, names(mf1))
  loc <- as.matrix(mf1[, m])
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations")

  ## Check corrfcn with loc
  if (corrfcn == "spherical" & NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }

  ## Split sample, prediction
  ii <- is.finite(yy)
  y <- yy[ii]
  n <- sum(ii)
  l <- ll[ii]
  l <- if (is.null(l)) rep.int(1.0, n) else as.double(l)
  if (any(!is.finite(l))) stop ("Non-finite values in the weights")
  if (any(l <= 0)) stop ("Non-positive weights not allowed")
  if (grepl("^binomial(\\..+)?$", family)) {
    l <- l - y # Number of failures
  }
  F <- FF[ii, , drop = FALSE]
  dm <- sp::spDists(loc[ii, , drop = FALSE], longlat = longlat)
  offset <- oofset[ii]

  ## Prior for ssq
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## Prior for beta
  betaprior <- getbetaprior(betm0, betQ0, p)
  betm0 <- betaprior$betm0
  betQ0 <- betaprior$betQ0

  np <- as.integer(np)
  if (np < 1) stop ("Input np must be at least 1.")
  fval <- double(npars)
  ssq <- 0                    # Not used

  f90 <- .Fortran("llikparscalc", fval, nu, phi, omg, kappa, npars,
                  as.double(y), as.double(l),
                  F, as.double(offset), betm0, betQ0, ssqdf, ssqsc,
                  dm, tsq, tsqdf, as.integer(n), as.integer(p),
                  as.integer(np), ssq, ifam, icf, PACKAGE = "geoBayes")
  list(par_vals = allpars, aloglik = f90[[1]])
}


##' Approximate log-likelihood maximisation
##'
##' Uses the "L-BFGS-B" method of the function
##' \code{\link[stats]{optim}} to maximise the log-likelihood for the
##' parameters \code{linkp}, \code{phi}, \code{omg}, \code{kappa}.
##' @title Log-likelihood maximisation
##' @param paroptim A named list with the components "linkp", "phi",
##'   "omg", "kappa". Each component must be numeric with length 1, 2,
##'   or 3 with elements in increasing order. If the compontent's
##'   length is 1, then the corresponding parameter is considered to
##'   be fixed at that value. If 2, then the two numbers denote the
##'   lower and upper bounds for the optimisation of that parameter
##'   (infinities are allowed). If 3, these correspond to lower bound,
##'   starting value, upper bound for the estimation of that
##'   parameter.
##' @param formula A representation of the model in the form
##'   \code{response ~ terms}.
##' @param family The distribution of the response.
##' @param data An optional data frame containing the variables in the
##'   model.
##' @param weights An optional vector of weights. Number of replicated
##'   samples for Gaussian and gamma, number of trials for binomial,
##'   time length for Poisson.
##' @param subset An optional vector specifying a subset of
##'   observations to be used in the fitting process.
##' @param offset See \code{\link[stats]{lm}}.
##' @param atsample A formula in the form \code{~ x1 + x2 + ... + xd}
##'   with the coordinates of the sampled locations.
##' @param corrfcn Spatial correlation function. See
##'   \code{\link{geoBayes_correlation}} for details.
##' @param np The number of integration points for the spatial
##'   variance parameter sigma^2. The total number of points will be
##'   \code{2*np + 1}.
##' @param betm0 Prior mean for beta (a vector or scalar).
##' @param betQ0 Prior standardised precision (inverse variance)
##'   matrix. Can be a scalar, vector or matrix. The first two imply a
##'   diagonal with those elements. Set this to 0 to indicate a flat
##'   improper prior.
##' @param ssqdf Degrees of freedom for the scaled inverse chi-square
##'   prior for the partial sill parameter.
##' @param ssqsc Scale for the scaled inverse chi-square prior for the
##'   partial sill parameter.
##' @param dispersion The fixed dispersion parameter.
##' @param longlat How to compute the distance between locations. If
##'   \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##'   distance. See \code{\link[sp]{spDists}}.
##' @param control A list of control parameters for the optimisation.
##'   See \code{\link[stats]{optim}}.
##' @return The output from the function \code{\link[stats]{optim}}.
##'   The \code{"value"} element is the log-likelihood, not the
##'   negative log-likelihood.
##' @importFrom stats optim model.frame optimHess model.offset
##' @importFrom optimx optimr
##' @useDynLib geoBayes llikparsval
##' @export
##' @references Evangelou, E., & Roy, V. (2019). Estimation and prediction for
##'   spatial generalized linear mixed models with parametric links
##'   via reparameterized importance sampling. Spatial Statistics, 29,
##'   289-315.
alik_optim <- function (paroptim,
                      formula,
                      family = "gaussian",
                      data, weights, subset, offset, atsample,
                      corrfcn = "matern",
                      np, betm0, betQ0, ssqdf, ssqsc, dispersion = 1,
                      longlat = FALSE, control = list())
{
  INPUT <- environment()
  log.phi <- FALSE
  log.omg <- FALSE
  ## Family
  ifam <- .geoBayes_family(family)
  if (ifam) {
    family <- .geoBayes_models$family[ifam]
  } else {
    stop ("This family has not been implemented.")
  }

  ## Logical input
  longlat <- as.logical(longlat)

  ## Correlation function
  icf <- .geoBayes_correlation(corrfcn)
  corrfcn <- .geoBayes_corrfcn$corrfcn[icf]

  ## Check dispersion input
  tsq <- as.double(dispersion)
  if (tsq <= 0) stop ("Invalid argument dispersion")
  tsqdf <- 0

  ## Design matrix and data
  if (missing(data)) data <- environment(formula)
  if (length(formula) != 3) stop ("The formula input is incomplete.")
  if ("|" == all.names(formula[[2]], TRUE, 1)) formula[[2]] <- formula[[2]][[2]]
  mfc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "offset"),
             names(mfc), 0L)
  mfc <- mfc[c(1L, m)]
  mfc$formula <- formula
  mfc$drop.unused.levels <- TRUE
  mfc$na.action <- "na.pass"
  mfc[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfc, parent.frame())
  mt <- attr(mf, "terms")
  FF <- model.matrix(mt,mf)
  if (!all(is.finite(FF))) stop ("Non-finite values in the design matrix")
  p <- NCOL(FF)
  yy <- unclass(model.response(mf))
  if (!is.vector(yy)) {
    stop ("The response must be a vector")
  }
  yy <- as.double(yy)
  ll <- model.weights(mf)
  oofset <- as.vector(model.offset(mf))
  if (!is.null(oofset)) {
    if (length(oofset) != NROW(yy)) {
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                    length(oofset), NROW(yy)), domain = NA)
    } else {
      oofset <- as.double(oofset)
    }
  } else {
    oofset <- double(NROW(yy))
  }

  ## All locations
  atsample <- update(atsample, NULL ~ .)
  locvars <- all.vars(atsample)
  formula1 <- as.formula(paste('~', paste(c(locvars, all.vars(formula)),
                                          collapse = ' + ')))
  mfc1 <- mfc
  mfc1$formula <- formula1
  mf1 <- eval(mfc1, parent.frame())
  m <- match(locvars, names(mf1))
  loc <- as.matrix(mf1[, m])
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations")

  ## Check corrfcn with loc
  if (corrfcn == "spherical" && NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }

  ## Split sample, prediction
  ii <- is.finite(yy)
  y <- yy[ii]
  n <- sum(ii)
  l <- ll[ii]
  l <- if (is.null(l)) rep.int(1.0, n) else as.double(l)
  if (any(!is.finite(l))) stop ("Non-finite values in the weights")
  if (any(l <= 0)) stop ("Non-positive weights not allowed")
  if (grepl("^binomial(\\..+)?$", family)) {
    l <- l - y # Number of failures
  }
  F <- FF[ii, , drop = FALSE]
  offset <- oofset[ii]
  dm <- sp::spDists(loc[ii, , drop = FALSE], longlat = longlat)

  ## Prior for ssq
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## Prior for beta
  betaprior <- getbetaprior(betm0, betQ0, p)
  betm0 <- betaprior$betm0
  betQ0 <- betaprior$betQ0

  np <- as.integer(np)
  if (np < 1) stop ("Input np must be at least 1.")
  ssq <- 0 # Not used

  ## Read and check paroptim
  paroptim <- getparoptim(paroptim, ifam, icf)
  pstart <- as.double(paroptim$pstart)
  lower <- paroptim$lower
  upper <- paroptim$upper
  estim <- paroptim$estim
  log.phi <- log.phi && estim[2]
  log.omg <- log.omg && estim[3]
  if (log.phi) {             # Estimating phi so transform
    pstart[2] <- log(pstart[2])
    lower[2] <- log(lower[2])
    upper[2] <- log(upper[2])
  }
  if (log.omg) {             # Estimating omg so transform
    pstart[3] <- log(pstart[3])
    lower[3] <- log(lower[3])
    upper[3] <- log(upper[3])
  }

  fn <- function (par)
  {
    parin <- pstart
    parin[estim] <- par
    nu <- parin[1]
    phi <- parin[2]; if(log.phi) phi <- exp(phi)
    omg <- parin[3]; if(log.omg) omg <- exp(omg)
    kappa <- parin[4]
    fval <- 0
    gval <- rep.int(0, 4)
    ideriv <- rep.int(0L, 4)
    f90 <- try(
      .Fortran("llikparsval", fval, gval, ideriv, nu, phi, omg, kappa,
               as.double(y), as.double(l),
               as.double(F), as.double(offset),
               as.double(betm0), as.double(betQ0),
               as.double(ssqdf), as.double(ssqsc),
               as.double(dm), as.double(tsq), as.double(tsqdf),
               as.integer(n), as.integer(p),
               as.integer(np), as.double(ssq), ifam, icf,
               PACKAGE = "geoBayes"),
      silent = TRUE)
    if (inherits(f90, "try-error")) return (NA)
    -f90[[1]]
  }
  gr <- function (par)
  {
    d1 <- c(1, 1, 1, 1)
    parin <- pstart
    parin[estim] <- par
    nu <- parin[1]
    phi <- parin[2]; if(log.phi) phi <- d1[2] <- exp(phi)
    omg <- parin[3]; if(log.omg) omg <- d1[3] <- exp(omg)
    kappa <- parin[4]
    fval <- 0
    gval <- rep.int(0, 4)
    ideriv <- as.integer(estim)
    f90 <- try(
      .Fortran("llikparsval", fval, gval, ideriv, nu, phi, omg, kappa,
               as.double(y), as.double(l),
               F, as.double(offset), betm0, betQ0, ssqdf, ssqsc,
               dm, tsq, tsqdf, as.integer(n), as.integer(p),
               as.integer(np), ssq, ifam, icf, PACKAGE = "geoBayes"),
      silent = TRUE)
    if (inherits(f90, "try-error")) return (rep(NA, sum(estim)))
    -f90[[2]][estim]*d1[estim]
  }
  method <- "L-BFGS-B" # if (sum(estim) == 1) "Brent" else "L-BFGS-B"
  method <- "nlminb"
  if (isTRUE(control$fnscale < 0)) control$fnscale <- -control$fnscale
###   op <- stats::optim(pstart[estim], fn, gr, method = method,
###                      lower = lower[estim], upper = upper[estim],
###                      hessian = TRUE, control = control)
  if (isTRUE(control$maximize)) control$maximize <- FALSE
  op <- optimx::optimr(pstart[estim], fn, gr, method = method,
                       lower = lower[estim], upper = upper[estim],
                       hessian = FALSE, control = control)
  op$gradient <- gr(op$par)
  op$hessian <- stats::optimHess(op$par, fn, gr, control = control)
  parout <- pstart
  parout[estim] <- op$par
  d1 <- c(1, 1, 1, 1)
  if (log.phi) d1[2] <- parout[2] <- exp(parout[2])
  if (log.omg) d1[3] <- parout[3] <- exp(parout[3])
  names(parout) <- names(d1) <- c("linkp", "phi", "omg", "kappa")
  op$par <- parout
  op$value <- -op$value
  op$hessian <- op$hessian*tcrossprod(d1[estim])
  op$INPUT <- INPUT
###   opscale <- if(!is.null(control$parscale)) control$parscale else 1
###   control$parscale <- NULL
###   op <- stats::nlminb(pstart[estim], fn, gr, scale = opscale,
###                       control = control,
###                       lower = lower[estim], upper = upper[estim])
  op
}


##' Calculate the likelihood approximation at different parameter
##' values. This function is useful for choosing the skeleton set.
##'
##' The input \code{par_vals} is meant to contain vector of parameter
##' values for each parameter. For each element in \code{par_vals},
##' the other parameters are set equal to the maximisers given in
##' \code{likopt} and the approximate likelihood is computed. The
##' cuttoff is calculated using linear interpolation provided by
##' \code{\link[stats]{approx}}.
##' @title Approximate log-likelihood calculation
##' @param likopt Output from the function \code{\link{alik_optim}}.
##' @param par_vals A named list with some of the components "linkp",
##'   "phi", "omg", "kappa".
##' @param likthreshold A threshold value proportion to calculate the
##'   cutoff. The cutoff will be calculated as that proportion
##'   relative to the maximum value of the log-likelihood.
##' @return A list with the log-likelihood approximation and cutoff values.
##' @export
##' @references Evangelou, E., & Roy, V. (2019). Estimation and prediction for
##'   spatial generalized linear mixed models with parametric links
##'   via reparameterized importance sampling. Spatial Statistics, 29,
##'   289-315.
alik_cutoff <- function (likopt, par_vals, likthreshold) {
  ## For each parameter value in par_vals, fix the other parameters at
  ## likopt and compute the approximate log-likelihood.
  if (likthreshold >= 1) stop ("Input likthreshold must < 1.")
  if (likthreshold < 0) stop ("Input likthreshold must >= 0.")
  parnm <- c("linkp", "phi", "omg", "kappa")
  par_vals <- par_vals[parnm]
  names(par_vals) <- parnm
  family <- get("family", likopt$INPUT)
  ifam <- .geoBayes_family(family)
  needlnk <- .geoBayes_models$needlinkp[ifam]
  corrfcn <- get("corrfcn", likopt$INPUT)
  icf <- .geoBayes_correlation(corrfcn)
  needkap <- .geoBayes_corrfcn$needkappa[icf]
  calcpar <- c(needlnk, TRUE, TRUE, needkap) & !sapply(par_vals, is.null)
  parmx <- as.list(likopt$par)
  likmx <- likopt$value
  lthr <- log(likthreshold)
  alk <- as.list(rep(likmx, 4)); names(alk) <- parnm
  cff <- parmx
  input <- substitute(
    alik_inla(par_vals = NULL, formula = formula, family = family,
              data = data, weights = weights, subset = subset,
              atsample = atsample, corrfcn = corrfcn, np = np,
              betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf,
              ssqsc = ssqsc, tsqdf = tsqdf, tsqsc = tsqsc,
              dispersion = dispersion, longlat = longlat), env = likopt$INPUT)
  for (i in 1:4) {
    if (calcpar[i]) {
      tmp <- parmx
      tmp[[i]] <- par_vals[[i]]
      input[[2]] <- as.data.frame(tmp)
      alk[[i]] <- eval(input, likopt$INPUT)$aloglik
    }
  }
  if (is.finite(lthr)) {
    for (i in 1:4) {
      if (calcpar[i]) {
        if (is.finite(lthr)) {
          cff[[i]] <- calc_cutoff(par_vals[[i]], alk[[i]] - likmx, lthr)
        } else {
          cff[[i]] <- c(min(par_vals[[i]]), max(par_vals[[i]]))
        }
      }
    }
  }
  list(par_vals = par_vals, aloglik = alk, cutoff = cff,
       threshold = likmx + lthr, likmax = likmx)
}

##' @importFrom stats approx
calc_cutoff <- function (x, y, y0) {
  i <- which.max(y)
  xmx <- x[i]
  i1 <- x <= xmx
  i2 <- x >= xmx
  x1 <- if (sum(i1) > 1) stats::approx(y[i1], x[i1], y0, rule = 2)$y else x[i1]
  x2 <- if (sum(i2) > 1) stats::approx(y[i2], x[i2], y0, rule = 2)$y else x[i2]
  c(x1, x2)
}

##' Plot likelihood approximation.
##'
##' The plot can be used to visualise the Laplace approximation to the
##' likelihood provided by the function \code{\link{alik_cutoff}}.
##' @title Plot likelihood approximation
##' @param alikobj Output from \code{\link{alik_cutoff}}.
##' @return Draws a plot.
##' @importFrom graphics segments
##' @rdname alik_cutoff
##' @export
alik_plot <- function (alikobj) {
  pvl <- alikobj$par_vals
  alk <- alikobj$aloglik
  cff <- alikobj$cutoff
  lmx <- alikobj$likmax
  thr <- exp(alikobj$threshold - lmx)
  lplt <- sapply(alk, function(x) length(x) > 1)
  npar <- sum(lplt)
  oldpar <- par(mfrow = c(1, npar)); on.exit(par(oldpar))
  parlb <- expression(linkp = nu, phi = phi, omg = omega, kappa = kappa)
  for (i in 1:4) {
    if (lplt[i]) {
      plot(pvl[[i]], exp(alk[[i]]-lmx), type = 'l', xlab = parlb[i], ylab = "")
      segments(c(cff[[i]][1], cff[[i]][1], cff[[i]][2]), c(0, thr, thr),
               c(cff[[i]][1], cff[[i]][2], cff[[i]][2]), c(thr, thr, 0),
               lty = 2)
    }
  }
}


## .. content for \description{} (no empty lines) .. TODO
##
## .. content for \details{} .. TODO
## @title Approximate log-likelihood calculation
## @param likopt Output from the function \code{\link{alik_optim}}.
## @param par_grid A \code{data.frame} with the components "linkp",
##   "phi", "omg", "kappa". The likelihood will be evalutated on each
##   row of this data frame.
## @param likthreshold A threshold value proportion to discard
##   points. The cutoff will be calculated as that proportion
##   relative to the maximum value of the log-likelihood.
## @param formula A representation of the model in the form
##   \code{response ~ terms}.
## @param family The distribution of the response. Can be one of the
##   options in \code{\link{.geoBayes_models}} or
##   \code{"transformed.gaussian"}.
## @param data An optional data frame containing the variables in the
##   model.
## @param weights An optional vector of weights. Number of replicated
##   samples for Gaussian and gamma, number of trials for binomial,
##   time length for Poisson.
## @param subset An optional vector specifying a subset of
##   observations to be used in the fitting process.
## @param atsample A formula in the form \code{~ x1 + x2 + ... + xd}
##   with the coordinates of the sampled locations.
## @param corrfcn Spatial correlation function. Can be one of the
##   choices in \code{\link{.geoBayes_corrfcn}}.
## @param np The number of integration points for the spatial
##   variance parameter sigma^2. The total number of points will be
##   \code{2*np + 1}.
## @param betm0 Prior mean for beta (a vector or scalar).
## @param betQ0 Prior standardised precision (inverse variance)
##   matrix. Can be a scalar, vector or matrix. The first two imply a
##   diagonal with those elements. Set this to 0 to indicate a flat
##   improper prior.
## @param ssqdf Degrees of freedom for the scaled inverse chi-square
##   prior for the partial sill parameter.
## @param ssqsc Scale for the scaled inverse chi-square prior for the
##   partial sill parameter.
## @param tsqdf Degrees of freedom for the scaled inverse chi-square
##   prior for the measurement error parameter.
## @param tsqsc Scale for the scaled inverse chi-square prior for the
##   measurement error parameter.
## @param dispersion The fixed dispersion parameter.
## @param longlat How to compute the distance between locations. If
##   \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##   distance. See \code{\link[sp]{spDists}}.
## @return A list with the log-likelihood approximation and cutoff values.
## calc_alik2 <- function (likopt, par_grid, likthreshold, formula,
##                         family = "gaussian",
##                         data, weights, subset, atsample,
##                         corrfcn = "matern",
##                         np, betm0, betQ0, ssqdf, ssqsc, tsqdf, tsqsc,
##                         dispersion = 1, longlat = FALSE) {
##   ## For each parameter value in par_grid, fix the other parameters at
##   ## likopt and compute the approximate log-likelihood.
##   if (likthreshold >= 1) stop ("Input likthreshold must < 1.")
##   par_grid <- as.data.frame(par_grid, stringsAsFactors = FALSE)
##   parnm <- c("linkp", "phi", "omg", "kappa")
##   ifam <- .geoBayes_family(family)
##   needlnk <- .geoBayes_models$needlinkp[ifam]
##   if (!needlnk) par_grid$linkp <- 0
##   icf <- .geoBayes_correlation(corrfcn)
##   needkap <- .geoBayes_corrfcn$needkappa[icf]
##   if (!needkap) par_grid$kappa <- 0
##   par_grid <- par_grid[, parnm]
##   parmx <- likopt$par
##   likmx <- likopt$value
##   lthr <- log(likthreshold)
##   input <- substitute(
##     alik_inla(linkp = pargrid[, 1], phi = pargrid[, 2], omg = pargrid[, 3],
##              kappa = pargrid[, 4], formula = formula, family = family,
##              data = data, weights = weights, subset = subset,
##              atsample = atsample, corrfcn = corrfcn, np = np,
##              betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf,
##              ssqsc = ssqsc, tsqdf = tsqdf, tsqsc = tsqsc,
##              dispersion = dispersion, longlat = longlat))
##   alk <- eval(input)
##
## }
