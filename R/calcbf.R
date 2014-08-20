##' Compute the Bayes factors.
##'
##' Computes the Bayes factors using the importance weights.
##' @title Compute the Bayes factors at other points
##' @param bfspobj Output from the function \code{\link{bf1skel}} which
##' contains the Bayes factors and importance sampling weights.
##' @param linkp,phi,omg,kappa Scalar or vector or \code{NULL}. If
##' scalar or vector, the Bayes factors are calculated at those
##' values. If \code{NULL} then the unique values from the MCMC
##' chains that were inputted in \code{\link{bf1skel}} will be used.
##' @param useCV Whether to use control variates for finer
##' corrections. 
##' @return An array of size length(linkp) x length(phi) x length(omg)
##' x length(kappa) containing the Bayes factors for each combination
##' of the parameters.
##' @importFrom sp spDists
##' @references Doss, H. (2010). Estimation of large families of Bayes
##' factors from Markov chain output. \emph{Statistica Sinica}, 20(2),
##' 537. 
##' @export 
bf2new <- function (bfspobj, linkp, phi, omg, kappa, useCV = TRUE) {
  ## Extract model variables
  y <- bfspobj$response
  n <- length(y)
  l <- bfspobj$weights
  F <- bfspobj$modelmatrix
  p <- NCOL(F)
  family <- bfspobj$family
  ifam <- match(family, c("gaussian", "binomial", "poisson", "Gamma"), 0L)
  corrfcn <- bfspobj$corrfcn
  icf <- match(corrfcn, c("matern", "spherical", "power"))
  betm0 <- bfspobj$betm0
  betQ0 <- bfspobj$betQ0
  ssqdf <- bfspobj$ssqdf
  ssqsc <- bfspobj$ssqsc
  dispersion <- bfspobj$dispersion
  tsqdf <- bfspobj$tsqdf
  tsqsc <- bfspobj$tsqsc
  tsq <- if (ifam == 0) tsqsc else dispersion
  loc <- bfspobj$locations
  dm <- sp::spDists(loc)
  z <- bfspobj$z
  Ntot <- NCOL(z)
  isweights <- bfspobj$isweights
  zcv <- bfspobj$controlvar
  kg <- NCOL(zcv)

  ## Examine link function parameter
  if (is.null(linkp)) {
    linkp <- nu <- unique(bfspobj$pnts$nu)
  } else if (family == "binomial") {
    if (is.character(linkp)) {
      if (length(linkp) != 1) {
        stop ("When using character linkp, it must not be a vector")
      }
      if (linkp == "logit") {
        nu <- unique(bfspobj$pnts$nu)
        if (!isTRUE(nu == -1)) {
          stop ("The logit link is not consistent with the model in bfspobj")
        }
      } else if (linkp == "probit") {
        nu <- unique(bfspobj$pnts$nu)
        if (!isTRUE(nu == 0)) {
          stop ("The logit link is not consistent with the model in bfspobj")
        }
      } else {
        stop ("Unrecognised linkp")
      }
    } else if (is.numeric(linkp)) {
      nu <- as.double(linkp)
      if (any(nu <= 0)) {
        stop ("Link parameter must be postive for the robit link")
      }
    } else {
      stop ("Unrecognised linkp")
    }
  } else if (is.numeric(linkp)) {
    nu <- as.double(linkp)
  } else {
    stop ("Unrecognised linkp")
  }
  n_nu <- length(nu)

  ## Examine covariance parameters
  if (is.null(phi)) {
    phi <- unique(bfspobj$pnts$phi)
  } else if (!is.numeric(phi)) {
    stop ("Argument phi must be numeric or NULL")
  }
  phi <- as.double(phi)
  if (any(phi < 0)) stop ("Argument phi must be non-negative")
  if (is.null(omg)) {
    omg <- unique(bfspobj$pnts$omg)
  } else if (!is.numeric(omg)) {
    stop ("Argument omg must be numeric or NULL")
  }
  n_phi <- length(phi)
  omg <- as.double(omg)
  if (any(omg < 0)) stop ("Argument omg must be non-negative")
  n_omg <- length(omg)
  if (is.null(kappa)) {
    kappa <- unique(bfspobj$pnts$kappa)
  } else if (!is.numeric(kappa)) {
    stop ("Argument kappa must be numeric or NULL")
  }
  kappa <- as.double(kappa)
  if (any(kappa < 0) & corrfcn %in% c("matern", "power")) {
    stop ("Argument kappa cannot be negative")
  }
  if (any(kappa > 2) & corrfcn == "power") {
    stop ("Argument kappa cannot be more than 2")
  }
  n_kappa <- length(kappa)
  covpars <- expand.grid(phi = phi, omg = omg, kappa = kappa)
  n_cov <- NROW(covpars)

  bfact <- numeric(n_nu*n_cov)

  if (useCV) {
    RUN <- .Fortran('calcb_cv', bfact, covpars$phi, nu, covpars$omg,
                    covpars$kappa, icf, 
                    n_cov, n_nu, Ntot, z, isweights, zcv, n, p, kg, betm0,
                    betQ0, ssqdf, ssqsc, max(tsqdf, 0), tsq, y, l, F, dm, ifam)
  } else {
    RUN <- .Fortran('calcb_st', bfact, covpars$phi, nu, covpars$omg,
                    covpars$kappa, icf, 
                    n_cov, n_nu, Ntot, z, isweights, n, p, betm0,
                    betQ0, ssqdf, ssqsc, max(tsqdf, 0), tsq, y, l, F, dm, ifam)
  }
  logbf <- array(RUN[[1]], c(n_nu, n_phi, n_omg, n_kappa))
  logbf <- logbf - bfspobj$referencebf
  maxid <- arrayInd(which.max(logbf), c(n_nu, n_phi, n_omg, n_kappa))
  out <- list(logbf = logbf, linkp = linkp, phi = phi,
              omg = omg, corrfcn = corrfcn, kappa = kappa, indmax = maxid)
  class(out) <- c("bfsp", "list")
  out
}


##' Creates a plot of the estimated Bayes factors from the function
##' \code{\link{bf2new}}. 
##' 
##' Depending on whether \code{ipar} is of length 1 or 2, this
##' function creates a line or a contour plot of the estimated Bayes
##' factors. 
##' @title Plot the estimated Bayes factors
##' @param x Output from the function \code{\link{bf2new}}.
##' @param ipar One or two unique integers in 1:4 indicating with
##' respect to which parameters to plot.
##' @param xlab,ylab,main Graphical parameters
##' @param ... Other graphical parameters to be passed on either
##' \code{plot} or \code{contour}.
##' @method contour bfsp
##' @return Produces the plot mentioned in the Details.
##' @export 
contour.bfsp <- function (x, ipar, xlab, ylab, main, ...) {
  ipar <- unique(as.integer(ipar))
  if (any(ipar < 1) | any(ipar > 4)) {
    stop("Argument ipar must be a vector of integers in 1:4")
  }
  npar <- length(ipar)
  if (length(ipar) < 1) stop ("Missing argument ipar")
  if (length(ipar) > 2) stop ("Cannot plot more than two parameters at a time")
  pars <- c("linkp", "phi", "omg", "kappa")
  maxid <- x$indmax
  ii <- c(list(), maxid)
  ii[ipar] <- TRUE
  bf <- do.call('[', c(list(x$logbf), ii, list(drop = FALSE)))
  jj <- 1:4; jj[sort(ipar)] <- ipar
  bf <- drop(aperm(bf, jj))
  if (length(ipar) == 1) {
    par1nm <- pars[ipar]
    par1 <- x[[par1nm]]
    if (missing(xlab)) xlab <- par1nm
    if (missing(ylab)) ylab <- "Logarithm of Bayes factor"
    if (missing(main)) main <- ""
    plot(par1, bf, type = "l", xlab = xlab, ylab = ylab, main = main, ...)
  } else {
    par1nm <- pars[ipar[1]]
    par2nm <- pars[ipar[2]]
    par1 <- x[[par1nm]]
    par2 <- x[[par2nm]]
    if (missing(xlab)) xlab <- par1nm
    if (missing(ylab)) ylab <- par2nm
    if (missing(main)) main <- "Logarithm of Bayes factor"
    contour(par1, par2, bf, xlab = xlab, ylab = ylab, main = main, ...)
    points(par1[maxid[ipar[1]]], par2[maxid[ipar[2]]])
  }
  invisible()
}


##' Estimation by empirical Bayes.
##'
##' This function is a wrap around \code{\link{bf2new}} using the
##' "L-BFGS-B" method of the function \code{\link[stats]{optim}} to
##' estimate the parameters.
##' @title Empirical Bayes estimator
##' @param bfspobj Output from the function \code{\link{bf1skel}} which
##' contains the Bayes factors and importance sampling weights.
##' @param estimate A named list with the components "linkp", "phi",
##' "omg", "kappa". Each component must be numeric with length 1, 2,
##' or 3 with elements in increasing order but for the binomial family
##' linkp is also allowed to be the character "logit" and "probit". If
##' its length is 1, then the corresponding parameter is considered to
##' be fixed at that value. If 2, then the two numbers denote the
##' lower and upper bounds for the optimisation of that parameter
##' (infinities are allowed). If 3, these correspond to lower bound,
##' starting value, upper bound for the estimation of that parameter.
##' @param useCV Whether to use control variates for finer
##' corrections.
##' @param control A list of control parameters for the optimisation.
##' See \code{\link[stats]{optim}}.
##' @return The output from the function \code{\link[stats]{optim}}.
##' @importFrom stats optim
##' @export 
bf2optim <- function (bfspobj, estimate, useCV = TRUE,
                      control = list()) {

  ## Extract model variables
  y <- bfspobj$response
  n <- length(y)
  l <- bfspobj$weights
  F <- bfspobj$modelmatrix
  p <- NCOL(F)
  family <- bfspobj$family
  ifam <- match(family, c("gaussian", "binomial", "poisson", "Gamma"), 0L)
  corrfcn <- bfspobj$corrfcn
  icf <- match(corrfcn, c("matern", "spherical", "power"))
  betm0 <- bfspobj$betm0
  betQ0 <- bfspobj$betQ0
  ssqdf <- bfspobj$ssqdf
  ssqsc <- bfspobj$ssqsc
  dispersion <- bfspobj$dispersion
  tsqdf <- bfspobj$tsqdf
  tsqsc <- bfspobj$tsqsc
  tsq <- if (ifam == 0) tsqsc else dispersion
  loc <- bfspobj$locations
  dm <- sp::spDists(loc)
  z <- bfspobj$z
  Ntot <- NCOL(z)
  isweights <- bfspobj$isweights
  zcv <- bfspobj$controlvar
  nruns <- NCOL(zcv)
  linkp <- if (family == "binomial") {
    if (bfspobj$pnts$nu[1] < 0) {
      "logit"
    } else if (bfspobj$pnts$nu[1] < 0) {
      "probit"
    } else "robit"
  } else "boxcox"
  parnm <- c("linkp", "phi", "omg", "kappa")
  
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
          stop ("Link in estimate inconsistent with link in bfspobj")
        }
      } else if (linkpe == "probit") {
        if (linkp == "probit") {
          pstart[1] <- 0
          estim[1] <- FALSE
        } else {
          stop ("Link in estimate inconsistent with link in bfspobj")
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

  ## Function to optimise
  fn <- if (useCV) {
    function (par) {
      parin <- pstart
      parin[estim] <- par
      RUN <- .Fortran('calcb_cv', 0.0, parin[2], parin[1], parin[3], parin[4],
                      icf, 1L, 1L, Ntot, z, isweights, zcv, n, p, nruns, betm0,
                      betQ0, ssqdf, ssqsc, max(tsqdf, 0), tsq, y, l, F, dm,
                      ifam)
      -RUN[[1]][1]
    }
  } else {
    function (par) {
      parin <- pstart
      parin[estim] <- par
      RUN <- .Fortran('calcb_st', 0.0, parin[2], parin[1], parin[3], parin[4],
                      icf, 1L, 1L, Ntot, z, isweights, n, p, betm0,
                      betQ0, ssqdf, ssqsc, max(tsqdf, 0), tsq, y, l, F, dm,
                      ifam)
      -RUN[[1]][1]
    }
  }
  
  method <- if (sum(estim) == 1) "Brent" else "L-BFGS-B"
  logbf <- bfspobj$logbf
  imaxlogbf <- which.max(logbf)
  wbf <- exp(logbf - logbf[imaxlogbf] - log(sum(exp(logbf - logbf[imaxlogbf]))))
  pstart.d <- colSums(data.frame(bfspobj$pnts[c("nu", "phi", "omg", "kappa")])*
                      wbf)
  i <- is.na(pstart) & estim
  pstart[i] <- pmax(pmin(upper[i], pstart.d[i]), lower[i])
  op <- stats::optim(pstart[estim], fn, method = method,
                     lower = lower[estim], upper = upper[estim],
                     control = control)
  parout <- pstart
  parout[estim] <- op$par
  op$par <- parout
  op$value <- -op$value
  op
}
