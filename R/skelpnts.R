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


##' Log-likelihood approximation
##'
##' Computes and approximation to the log-likelihood for the given
##' parameters.
##' @title Log-likelihood approximation
##' @param linkp Parameter of the link function. For binomial, a
##'   positive number for the degrees of freedom of the robit family
##'   or "logit" or "probit". For the other families any number for
##'   the exponent of the Box-Cox transformation. Input can be a
##'   scalar or a vector.
##' @param phi Spatial range parameter. Input can be a scalar or a
##'   vector.
##' @param omg Relative nugget parameter. Input can be a scalar or a
##'   vector.
##' @param kappa Spatial smoothness parameter. Input can be a scalar
##'   or a vector.
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
##' @return A vector of the same length as the parameters containing
##'   the log-likelihood values.
##' @useDynLib geoBayes llikparscalc
##' @export
likaprxn <- function (linkp, phi, omg, kappa, formula,
                      family = "gaussian",
                      data, weights, subset, atsample,
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
  kappa <- .geoBayes_getkappa(kappa, icf)
  phi <- as.double(phi)
  if (any(phi < 0)) stop ("Input phi must >= 0")
  omg <- as.double(omg)
  if (any(omg < 0)) stop ("Input omg must >= 0")

  ## Check the link parameter
  nu <- .geoBayes_getlinkp(linkp, ifam)

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
  m <- match(c("formula", "data", "subset", "weights"),
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

  ## Prior for ssq
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## Prior for beta
  if (all(is.finite(betQ0[upper.tri(betQ0)]))) {
    if (length(betQ0) == 1 && betQ0[1] == 0) {
      ## Uniform prior
      betQ0 <- matrix(0, p, p)
      betm0 <- rep(0, p)
    } else if (length(betQ0) == 1 || length(betQ0) == p) {
      if (any(betQ0 <= 0)) stop ('betQ0 not > 0')
      betQ0 <- diag(betQ0, p, p)
      betm0 <- rep(as.double(betm0), length.out = p)
      modeldf <- as.double(n + ssqdf)
    } else if (length(betQ0) == p*p) {
      betQ0 <- matrix(as.double(betQ0), p, p)
      betQ0[lower.tri(betQ0)] <- 0
      betQ0eig <- eigen(t(betQ0), 1, 1)$values
      if (any (betQ0eig < sqrt(.Machine$double.eps))) {
        stop ('betQ0 not > 0 within tolerance')
      }
      betm0 <- rep(as.double(betm0), length.out = p)
      modeldf <- as.double(n + ssqdf)
    } else stop ('Bad betQ0')
  } else stop ('Non-finite betQ0')

  np <- as.integer(np)
  if (np < 1) stop ("Input np must be at least 1.")
  fval <- double(npars)
  ssq <- 0                    # Not used

  f90 <- .Fortran("llikparscalc", fval, nu, phi, omg, kappa, npars,
                  as.double(y), as.double(l),
                  F, betm0, betQ0, ssqdf, ssqsc,
                  dm, tsq, tsqdf, as.integer(n), as.integer(p),
                  as.integer(np), ssq, ifam, icf, PACKAGE = "geoBayes")
  f90[[1]]
}


##' Log-likelihood maximisation
##'
##' Uses the "L-BFGS-B" method of the function
##' \code{\link[stats]{optim}} to maximise the log-likelihood for the
##' parameters \code{linkp}, \code{phi}, \code{omg}, \code{kappa}.
##' @title Log-likelihood maximisation
##' @param paroptim A named list with the components "linkp", "phi",
##'   "omg", "kappa". Each component must be numeric with length 1, 2,
##'   or 3 with elements in increasing order but for the binomial
##'   family linkp is also allowed to be the character "logit" and
##'   "probit". If the compontent's length is 1, then the
##'   corresponding parameter is considered to be fixed at that value.
##'   If 2, then the two numbers denote the lower and upper bounds for
##'   the optimisation of that parameter (infinities are allowed). If
##'   3, these correspond to lower bound, starting value, upper bound
##'   for the estimation of that parameter.
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
##' @importFrom stats optim model.frame
##' @useDynLib geoBayes llikparsval
##' @export
likoptim <- function (paroptim,
                      formula,
                      family = "gaussian",
                      data, weights, subset, atsample,
                      corrfcn = "matern",
                      np, betm0, betQ0, ssqdf, ssqsc, dispersion = 1,
                      longlat = FALSE, control = list())
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
  m <- match(c("formula", "data", "subset", "weights"),
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
  dm <- sp::spDists(loc[ii, , drop = FALSE], longlat = longlat)

  ## Prior for ssq
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## Prior for beta
  if (all(is.finite(betQ0[upper.tri(betQ0)]))) {
    if (length(betQ0) == 1 && betQ0[1] == 0) {
      ## Uniform prior
      betQ0 <- matrix(0, p, p)
      betm0 <- rep(0, p)
    } else if (length(betQ0) == 1 || length(betQ0) == p) {
      if (any(betQ0 <= 0)) stop ('betQ0 not > 0')
      betQ0 <- diag(betQ0, p, p)
      betm0 <- rep(as.double(betm0), length.out = p)
      modeldf <- as.double(n + ssqdf)
    } else if (length(betQ0) == p*p) {
      betQ0 <- matrix(as.double(betQ0), p, p)
      betQ0[lower.tri(betQ0)] <- 0
      betQ0eig <- eigen(t(betQ0), 1, 1)$values
      if (any (betQ0eig < sqrt(.Machine$double.eps))) {
        stop ('betQ0 not > 0 within tolerance')
      }
      betm0 <- rep(as.double(betm0), length.out = p)
      modeldf <- as.double(n + ssqdf)
    } else stop ('Bad betQ0')
  } else stop ('Non-finite betQ0')

  np <- as.integer(np)
  if (np < 1) stop ("Input np must be at least 1.")
  ssq <- 0 # Not used

  ## Read and check paroptim
  paroptim <- getparoptim(paroptim, ifam, icf)
  pstart <- as.double(paroptim$pstart)
  lower <- paroptim$lower
  upper <- paroptim$upper
  estim <- paroptim$estim

  fn <- function (par)
  {
    parin <- pstart
    parin[estim] <- par
    nu <- parin[1]
    phi <- parin[2]
    omg <- parin[3]
    kappa <- parin[4]
    fval <- 0
    f90 <- .Fortran("llikparsval", fval, nu, phi, omg, kappa,
                    as.double(y), as.double(l),
                    F, betm0, betQ0, ssqdf, ssqsc,
                    dm, tsq, tsqdf, as.integer(n), as.integer(p),
                    as.integer(np), ssq, ifam, icf, PACKAGE = "geoBayes")
    -f90[[1]]
  }
  gr <- NULL
  method <- if (sum(estim) == 1) "Brent" else "L-BFGS-B"
  if (isTRUE(control$fnscale < 0)) control$fnscale <- -control$fnscale
  op <- stats::optim(pstart[estim], fn, gr, method = method,
                     lower = lower[estim], upper = upper[estim],
                     hessian = TRUE, control = control)
  parout <- pstart
  parout[estim] <- op$par
  names(parout) <- c("linkp", "phi", "omg", "kappa")
  op$par <- parout
  op$value <- -op$value
  op
}
