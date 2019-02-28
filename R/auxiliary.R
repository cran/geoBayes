######################################################################
##
### Commentary: Some auxiliary functions.
##
######################################################################

getbetaprior <- function (betm0, betQ0, p) {
  if (length(betQ0) == 1 && betQ0[1] == 0) {
    ## Uniform prior
    betQ0 <- matrix(0, p, p)
    betm0 <- rep(0, p)
  } else if (length(betQ0) == 1 || length(betQ0) == p) {
    if (any(betQ0 <= 0)) stop ('betQ0 not > 0')
    betQ0 <- diag(as.vector(betQ0), p, p)
    betm0 <- rep(as.double(betm0), length.out = p)
  } else if (length(betQ0) == p*p) {
    betQ0 <- matrix(as.double(betQ0), p, p)
    betQ0[lower.tri(betQ0)] <- 0
    betQ0eig <- eigen(t(betQ0), 1, 1)$values
    if (any (betQ0eig < sqrt(.Machine$double.eps))) {
      stop ('betQ0 not > 0 within tolerance')
    }
    betm0 <- rep(as.double(betm0), length.out = p)
  } else stop ('Bad betQ0')
  if (!all(is.finite(betm0))) {
    stop ('Non-finite betm0')
  }
  if (!all(is.finite(betQ0[upper.tri(betQ0, diag = TRUE)]))) {
    stop ('Non-finite betQ0')
  }
  list(betm0 = as.double(betm0), betQ0 = as.double(betQ0))
}

## Checks if the paroptim argument is valid and returns a list with
## components pstart, lower, upper, estim.
getparoptim <- function (paroptim, ifam, icf) {
  paroptim <- as.list(paroptim)
  parnmall <- c("linkp", "phi", "omg", "kappa")
  if (!is.numeric(ifam)) ifam <- .geoBayes_family(ifam)
  if (!is.numeric(icf)) icf <- .geoBayes_correlation(icf)
  paroptim <- paroptim[parnmall]
  lower <- upper <- pstart <- rep.int(0, 4) # nu, phi, omg, kappa
  estim <- rep.int(FALSE, 4)

  ##
  needlinkp <- .geoBayes_models$needlinkp[ifam]
  needkappa <- .geoBayes_corrfcn$needkappa[icf]
  need <- c(needlinkp, TRUE, TRUE, needkappa)

  getnonneg <- function(x) {
    if (any(x < 0)) stop ("All values must be non-negative.")
    as.double(x)
  }
  getf <- list(linkp = function(x) .geoBayes_getlinkp(x, ifam),
               phi = getnonneg, omg = getnonneg,
               kappa = function(x) .geoBayes_getkappa(x, icf))

  ## Check parameters
  for (i in 1:4) {
    if (!need[i]) next        # Skip those not needed
    ppp <- paroptim[[parnmall[i]]]
    if (!is.numeric(ppp)) {
      stop(paste("The element", parnmall[i], "in paroptim must be numeric"))
    }
    ppp <- getf[[i]](ppp)
    lppp <- length(ppp)
    if (lppp < 1 || lppp > 3) {
      stop (paste("The element", parnmall[i],
                  "in paroptim must have 1, 2, or 3 components"))
    }
    if (lppp == 1) {
      pstart[i] <- ppp
      estim[i] <- FALSE
    } else if (lppp == 2) {
      pstart[i] <- .5*(ppp[1]+ppp[2])
      estim[i] <- TRUE
      lower[i] <- ppp[1]
      upper[i] <- ppp[2]
      if (lower[i] >= upper[i]) {
        stop (paste("The lower bound must be less than the upper bound for",
                     parnmall[i], "in paroptim"))
      }
    } else {
      pstart[i] <- ppp[2]
      estim[i] <- TRUE
      lower[i] <- ppp[1]
      upper[i] <- ppp[3]
      if (lower[i] > pstart[i] || pstart[i] > upper[i] || lower[i] == upper[i]) {
        stop (paste("The elements in the component", parnmall[i],
                     "in paroptim must be ordered"))
      }
    }
  }
  list(pstart = pstart, lower = lower, upper = upper, estim = estim)
}


## Take the sample sizes in N and split them according to fun
## using sz as a parameter.
## allow0 will allow sizes to be 0.
getsize <- function(sz, N, FUN = "*", allow0 = FALSE)
{
  nm <- deparse(substitute(sz))
  FUN <- match.fun(FUN)
  sz <- as.double(sz)
  nruns <- length(N)
  if (length(sz) > nruns) {
    warning(paste("The number of elements in", nm, "exceeds the
number of runs; the extra elements will be discarded."))
  }
  sz <- rep(sz, length.out = nruns)
  if (any(sz <= 0)) {
    stop(paste("Argument", nm, "must be positive"))
  } else if (all(sz < 1)) {
    Nout <- as.integer(FUN(N, sz))
    if (!allow0 && any(Nout == 0)) {
      stop(paste("Calculated 0 sizes; give a larger", nm))
    }
  } else if (all(sz == 1)) {
    Nout <- as.integer(FUN(N, sz))
  } else if (all(sz >= 1)) {
    Nout <- as.integer(sz)
    if (any(Nout > N)) {
      stop(paste(nm, "is too big."))
    }
  } else {
    stop (paste("Argument", nm, "must be all < 1 or all >= 1."))
  }
  Nout
}

## Compute transformed sample.
##
## Computes the transformed sample according to transf.
## @title Compute transformed sample
## @param runs A list of outputs from \code{mcsglmm} or \code{mctrga}.
## @param model A list with elements response and family.
## @param transf The type of transformation to use.
## @return A list with elements sample, transf, itr, ifam.
## @useDynLib geoBayes transformz
transfsample <- function (runs, model, transf = c("no", "mu", "wo"), verb = TRUE)
{
  y <- model$response
  n <- as.integer(length(y))
  family <- model$family
  ifam <- .geoBayes_family(family)
  if (!is.logical(transf)) {
    transf <- transf0 <- match.arg(transf)
  } else if (transf) {
    transf <- transf0 <- "mu"
  } else {
    transf <- transf0 <- "no"
  }
  if (transf == "wo" && !.geoBayes_models$haswo[ifam]) {
    if (verb)
      warning("Can't use workaround for this family. Computing no transformation.")
    transf0 <- "no"
  }
  if (transf == "wo" && family == "poisson.boxcox") {
    ypo <- y > 0
    if (all(ypo)) {
      if (verb)
        message("For poisson.boxcox with all positive responses, the workaround is the same as the mean transformation.")
      transf0 <- "mu"
    } else if (any(ypo)) {
      transf0 <- "tr"
    }
  }
  if (transf0 == "mu") {
    itr <- rep.int(1L, n)
    ftrw <- function (r) {
      o <- r$MCMC[["whichobs"]]
      out <- r$MCMC[["mu"]][o, , drop = FALSE]
      if (is.null(out) && any(o)) {
        out <- r$MCMC[["z"]][o, , drop = FALSE]
        nu <- r$FIXED[["linkp_num"]][1]
        nn <- as.integer(length(out))
        out[] <- .Fortran("transformz", out, nu, nn, ifam,
                          PACKAGE = "geoBayes")[[1]]
      }
      out
    }
    sample <- lapply(runs, ftrw)
  } else if (transf0 == "tr") {
    ifam <- -ifam
    if (family == "poisson.boxcox") { # Poisson y = 0
      itr <- 2L - ypo
      ftrw <- function (r) {
        i2 <- itr == 2
        o <- r$MCMC[["whichobs"]]
        out <- r$MCMC[["mu"]][o, , drop = FALSE]
        z <- r$MCMC[["z"]][o, ][i2, ]
        nu <- r$FIXED[["linkp_num"]][1]
        nn <- as.integer(length(z))
        out[i2, ] <- .Fortran("transformz", z, nu, nn, ifam,
                              PACKAGE = "geoBayes")[[1]]
        out
      }
    } else { # All other cases
      stop ("BUG: Only the poisson.boxcox family is treated differently.")
    }
    sample <- lapply(runs, ftrw)
  } else if (transf0 == "wo") {
    ifam <- -ifam
    itr <- rep.int(2L, n)
    ftrw <- function (r) {
      o <- r$MCMC[["whichobs"]]
      out <- r$MCMC[["z"]][o, , drop = FALSE]
      nu <- r$FIXED[["linkp_num"]][1]
      nn <- as.integer(length(out))
      out[] <- .Fortran("transformz", out, nu, nn, ifam,
                        PACKAGE = "geoBayes")[[1]]
      out
    }
    sample <- lapply(runs, ftrw)
  } else {
    ftrw <- function (r) {
      o <- r$MCMC[["whichobs"]]
      out <- r$MCMC[["z"]][o, , drop = FALSE]
      out
    }
    sample <- lapply(runs, ftrw)
    itr <- rep.int(0L, n)
  }
  list(sample = sample, transf = transf,
       real_transf = transf0, itr = as.integer(itr),
       ifam = ifam)
}



## Check and return the parameter grid
check_pargrid <- function(pargrid, family, corrfcn)
{
  pargrid <- as.data.frame(pargrid, stringsAsFactors = FALSE)
  ##if(!is.list(pargrid)) stop ("Argument pargrid must be a list")
  parnm <- c("linkp", "phi", "omg", "kappa")
  ifam <- .geoBayes_family(family)
  needlinkp <- .geoBayes_models$needlinkp[ifam]
  if (!needlinkp) pargrid$linkp <- 0
  icf <- .geoBayes_correlation(corrfcn)
  needkappa <- .geoBayes_corrfcn$needkappa[icf]
  if (!needkappa) {
    pargrid$kappa <- 0
  } else {
    pargrid$kappa <- .geoBayes_getkappa(pargrid$kappa, icf)
  }
  if (!all(parnm %in% names(pargrid))) {
    stop (paste("Argument", deparse(substitute(pargrid)), "must have the names",
                paste(c(if (needlinkp) "linkp", "phi", "omg",
                        if (needkappa) "kappa"), collapse=", ")))
  } else {
    pargrid <- pargrid[parnm]
  }
  nruns <- nrow(pargrid)
  phi <- pargrid[["phi"]]
  omg <- pargrid[["omg"]]
  kappa <- pargrid[["kappa"]]
  linkp <- pargrid[["linkp"]]
  if (any (phi < 0)) stop ("Element phi given must be non-negative")
  if (any (omg < 0)) stop ("Element omg given must be non-negative")

  ## Check if linkp conforms with family
  nu <- .geoBayes_getlinkp(linkp, family)

  ## Output
  data.frame(linkp = linkp, phi = phi, omg = omg, kappa = kappa, nu = nu)
}

check_gengamma_prior <- function (pars)
  ## Checks if the generalised gamma prior parameters are correct.
{
  pars <- as.double(pars)
  if (4 != length(pars))
    stop ("Generalised gamma parameters must be a vector of length 4.")
  if (!all(is.finite(pars)))
    stop ("Non-finite generalised gamma parameters.")
  if (pars[1] <= 0 || pars[4] < 0 || pars[2]*pars[3] < 0) {
    stop ("Invalid values in pars: pars[1] <= 0 || pars[4] < 0 || pars[2]*pars[3] < 0")
  }
  pars
}

check_unif_prior <- function (pars)
  ## Checks if the uniform prior parameters are correct.
{
  pars <- as.double(pars)
  if (2 != length(pars))
    stop ("Uniform parameters must be a vector of length 2.")
  if (!all(is.finite(pars)))
    stop ("Non-finite uniform prior parameters.")
  if (pars[1] >= pars[2]) {
    stop ("Invalid values in pars: pars[1] >= pars[2]")
  }
  pars
}
