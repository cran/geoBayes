######################################################################
##
### Commentary: Some auxiliary functions.
##
######################################################################


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
      pstart[i] <- NA ## .5*(ppp[1]+ppp[2])
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

##' Compute transformed sample.
##'
##' Computes the transformed sample according to transf.
##' @title Compute transformed sample
##' @param runs A list of lists with elements z, mu, nu, whichobs
##'   which gives the samples and link parameter.
##' @param model A list with elements response and family.
##' @param transf The type of transformation to use.
##' @return A list with elements sample, transf, itr, ifam.
##' @useDynLib geoBayes transformz
transfsample <- function (runs, model, transf = c("no", "mu", "wo"))
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
    warning("Can't use workaround for this family. Computing no transformation.")
    transf0 <- "no"
  }
  if (transf == "wo" && family == "poisson.boxcox") {
    ypo <- y > 0
    if (all(ypo)) {
      message("For poisson.boxcox with all positive responses, the workaround is the same as the mean transformation.")
      transf0 <- "mu"
    } else if (any(ypo)) {
      transf0 <- "tr"
    }
  }
  if (transf0 == "mu") {
    sample <- lapply(runs, function(r) r[["mu"]][r[["whichobs"]], ])
    itr <- rep.int(1L, n)
  } else if (transf0 == "tr") {
    ifam <- -ifam
    if (family == "poisson.boxcox") { # Poisson y = 0
      itr <- 2L - ypo
      ftrw <- function (r) {
        i2 <- itr == 2
        o <- r[["whichobs"]]
        out <- r[["mu"]][o, ]
        z <- r[["z"]][o, ][i2, ]
        nu <- r[["nu"]][1]
        nn <- as.integer(length(z))
        out[i2, ] <- .Fortran("transformz", z, nu, nn, as.integer(ifam),
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
      o <- r[["whichobs"]]
      out <- r[["z"]][o, ]
      nu <- r[["nu"]][1]
      nn <- as.integer(length(out))
      out[] <- .Fortran("transformz", out, nu, nn, as.integer(ifam),
                        PACKAGE = "geoBayes")[[1]]
      out
    }
    sample <- lapply(runs, ftrw)
  } else {
    sample <- lapply(runs, function(r) r[["z"]][r[["whichobs"]], ])
    itr <- rep.int(0L, n)
  }
  list(sample = sample, transf = transf,
       real_transf = transf0, itr = as.integer(itr),
       ifam = as.integer(ifam))
}



## Check and return the parameter grid
check_pargrid <- function(pargrid, family, corrfcn)
{
  pargrid <- data.frame(pargrid, stringsAsFactors = FALSE)
  ##if(!is.list(pargrid)) stop ("Argument pargrid must be a list")
  parnm <- c("linkp", "phi", "omg", "kappa")
  needkappa <- corrfcn %in% c("matern", "powerexponential")
  if (!needkappa) pargrid$kappa <- 0
  if (!all(parnm %in% names(pargrid))) {
    stop (paste("Argument", deparse(substitute(pargrid)), "must have the names",
                paste(c("linkp", "phi", "omg", if (needkappa) "kappa"),
                      collapse=" ")))
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

  ## Check kappa and corrfcn
  kappa <- as.double(kappa)
  if ((corrfcn %in% c("matern", "powerexponential")) && any(kappa < 0)) {
    stop ("Component kappa cannot be negative")
  }
  if ((corrfcn == "powerexponential") && any(kappa > 2)) {
    stop ("Component kappa cannot be more than 2 for the powerexponential correlation")
  }

  ## Check if linkp conforms with family
  nu <- .geoBayes_getlinkp(linkp, family)

  ## Output
  data.frame(linkp = linkp, phi = phi, omg = omg, kappa = kappa, nu = nu)
}
