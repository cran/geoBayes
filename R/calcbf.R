##' Compute the Bayes factors.
##'
##' Computes the Bayes factors using the importance weights at the new
##' points. The new points are taken from the grid derived by
##' expanding the parameter values inputted. The arguments
##' \code{linkp} \code{phi} \code{omg} \code{kappa} correspond to the
##' link function, spatial range, relative nugget, and correlation
##' function parameters respectively.
##' @title Compute the Bayes factors at new points
##' @param bf1obj Output from the function \code{\link{bf1skel}} which
##' contains the Bayes factors and importance sampling weights.
##' @param linkp,phi,omg,kappa Optional scalar or vector or
##' \code{NULL}. If scalar or vector, the Bayes factors are calculated
##' at those values with respect to the reference model used in
##' \code{\link{bf1skel}}. If missing or \code{NULL} then the unique
##' values from the MCMC chains that were inputted in
##' \code{\link{bf1skel}} will be used.
##' @param useCV Whether to use control variates for finer
##' corrections.
##' @return An array of size \code{length(linkp) * length(phi) *
##' length(omg) * length(kappa)} containing the Bayes factors for each
##' combination of the parameters.
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
##' bfall <- bf2new(bf, phi = seq(100, 200, 10), omg = seq(0, 2, .2))
##' plotbf2(bfall, c("phi", "omg"))
##' }
##' @importFrom sp spDists
##' @references Doss, H. (2010). Estimation of large families of Bayes
##' factors from Markov chain output. \emph{Statistica Sinica}, 20(2),
##' 537.
##'
##' Roy, V., Evangelou, E., and Zhu, Z. (2015). Efficient estimation
##' and prediction for the Bayesian spatial generalized linear mixed
##' model with flexible link functions. \emph{Biometrics}, 72(1),
##'   289-298.
##' @useDynLib geoBayes calcb_no_st calcb_mu_st calcb_wo_st calcb_tr_st
##' @useDynLib geoBayes calcb_no_cv calcb_mu_cv calcb_wo_cv calcb_tr_cv
##' @export
bf2new <- function (bf1obj, linkp, phi, omg, kappa, useCV = TRUE)
{
  ## Logical input
  useCV <- as.logical(useCV)

  ## Extract model variables
  transf <- bf1obj$transf
  real_transf <- bf1obj$real_transf
  itr <- bf1obj$itr
  y <- bf1obj$response
  n <- length(y)
  l <- bf1obj$weights
  F <- bf1obj$modelmatrix
  p <- NCOL(F)
  family <- bf1obj$family
  ifam <- .geoBayes_family(family)
  corrfcn <- bf1obj$corrfcn
  icf <- .geoBayes_correlation(corrfcn)
  betm0 <- bf1obj$betm0
  betQ0 <- bf1obj$betQ0
  ssqdf <- bf1obj$ssqdf
  ssqsc <- bf1obj$ssqsc
  dispersion <- bf1obj$dispersion
  tsqdf <- bf1obj$tsqdf
  tsqsc <- bf1obj$tsqsc
  tsq <- if (ifam == 0) tsqsc else dispersion
  loc <- bf1obj$locations
  dm <- sp::spDists(loc)
  sample <- bf1obj$sample2
  Ntot <- NCOL(sample)
  isweights <- bf1obj$isweights
  zcv <- bf1obj$controlvar
  kg <- NCOL(zcv)

  ## Examine link function parameter
  if (missing(linkp) || is.null(linkp)) {
    nu <- unique(bf1obj$pnts$nu)
  } else {
    nu <- unique(.geoBayes_getlinkp(linkp, ifam))
  }
  n_nu <- length(nu)

  ## Examine covariance parameters
  if (missing(phi) || is.null(phi)) {
    phi <- unique(bf1obj$pnts$phi)
  } else if (!is.numeric(phi)) {
    stop ("Argument phi must be numeric or NULL")
  }
  phi <- as.double(phi)
  if (any(phi < 0)) stop ("Argument phi must be non-negative")
  n_phi <- length(phi)
  if (missing(omg) || is.null(omg)) {
    omg <- unique(bf1obj$pnts$omg)
  } else if (!is.numeric(omg)) {
    stop ("Argument omg must be numeric or NULL")
  }
  omg <- as.double(omg)
  if (any(omg < 0)) stop ("Argument omg must be non-negative")
  n_omg <- length(omg)
  if (missing(kappa) || is.null(kappa)) {
    kappa <- unique(bf1obj$pnts$kappa)
  } else {
    kappa <- unique(.geoBayes_getkappa(kappa, icf))
  }
  n_kappa <- length(kappa)
  covpars <- expand.grid(phi = phi, omg = omg, kappa = kappa)
  n_cov <- NROW(covpars)

  bfact <- numeric(n_nu*n_cov)

  ## Check for non-finite values in control variates
  if (useCV && any(!is.finite(zcv))) {
    warning ("Computed non-finite control variates, probably due to numerical
overflow. Control variates corrections will not be used.")
    useCV <- FALSE
  }

  froutine <- paste0("calcb_", real_transf)
  if (real_transf == "wo" || real_transf == "tr") {
    ifam <- -ifam
  }

  if (useCV) {
    froutine <- paste(froutine, '_cv', sep='')
    QRin <- QRintercept(zcv)
    RUN <- .Fortran(froutine, bfact,
                    as.double(covpars$phi), as.double(nu),
                    as.double(covpars$omg),
                    as.double(covpars$kappa), as.integer(icf),
                    as.integer(n_cov), as.integer(n_nu), as.integer(Ntot),
                    as.double(sample), as.double(isweights),
                    as.double(QRin), as.integer(n), as.integer(p),
                    as.double(betm0),
                    as.double(betQ0), as.double(ssqdf), as.double(ssqsc),
                    max(tsqdf, 0), as.double(tsq), as.double(y),
                    as.double(l), as.double(F), as.double(dm),
                    as.integer(ifam), as.integer(itr),
                    PACKAGE = "geoBayes")
  } else {
    froutine <- paste(froutine, '_st', sep='')
    RUN <- .Fortran(froutine, bfact,
                    as.double(covpars$phi), as.double(nu),
                    as.double(covpars$omg),
                    as.double(covpars$kappa), as.integer(icf),
                    n_cov, n_nu, Ntot, as.double(sample), as.double(isweights),
                    as.integer(n), as.integer(p), as.double(betm0),
                    as.double(betQ0), as.double(ssqdf), as.double(ssqsc),
                    max(tsqdf, 0), as.double(tsq), as.double(y), as.double(l),
                    as.double(F), as.double(dm), as.integer(ifam),
                    as.integer(itr),
                    PACKAGE = "geoBayes")
  }

  logbf <- array(RUN[[1]], c(n_nu, n_phi, n_omg, n_kappa))
  logbf <- logbf + bf1obj$logbf[1] # Constant to match bf1obj$logbf
  maxid <- arrayInd(which.max(logbf), c(n_nu, n_phi, n_omg, n_kappa))
  out <- list(logbf = logbf, linkp = nu, phi = phi,
              omg = omg, corrfcn = corrfcn, kappa = kappa, indmax = maxid)
  class(out) <- c("bfsp", "list")
  out
}


##' This function plots the estimated logarithm Bayes factors from the
##' function \code{\link{bf2new}}.
##'
##' Depending on whether \code{pars} has length 1 or 2, this function
##' creates a line or a contour plot of the estimated Bayes factors.
##' If its length is 3 or 4, then it produces multiple profile plots.
##' In this case the variable is fixed at different values and the
##' maximum Bayes factor corresponding to the fixed value is plotted
##' against that value.
##' @title Plot the estimated Bayes factors
##' @param bf2obj Output from the function \code{\link{bf2new}}.
##' @param pars A vector with the names of the parameters to plot.
##' @param profile Whether it should produce a profile plot or a
##' contour plot if the length of pars is 2.
##' @param ... Other input to be passed to either \code{plot} or
##' \code{contour}.
##' @return This function returns nothing.
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
##' bfall <- bf2new(bf, phi = seq(100, 200, 10), omg = seq(0, 2, .2))
##' plotbf2(bfall, c("phi", "omg"))
##' plotbf2(bfall, c("phi", "omg"), profile = TRUE, type = "b", ylab="log(BF)")
##' }
##' @importFrom graphics plot points par title
##' @export
plotbf2 <- function (bf2obj, pars = c("linkp", "phi", "omg", "kappa"),
                     profile = length(pars) > 2, ...) {
  N <- 4
  pars <- unique(pars)
  pars <- match.arg(pars, several.ok = TRUE)
  usepars <- sapply(bf2obj[pars], length) > 1
  pars <- pars[usepars]
  ipar <- match(pars, eval(formals()$pars))
  npars <- length(pars)
  if (npars == 0) stop ("No parameters to plot.")
  dots <- list(...)
  nmdots <- names(dots)

  ## Determine type of plot
  if (profile) {
    ptype <- "profile"
  } else if (npars == 1) {
    ptype <- "line"
  } else if (npars == 2) {
    ptype <- "contour"
  } else {
    warning ("Only a profile plot is available for more than 3 parameters")
    ptype <- "profile"
  }

  grknm <- list(nu = expression(nu), phi = expression(phi),
                omg = expression(omega), kappa = expression(kappa))

  if (ptype == "line") {
    maxid <- bf2obj$indmax
    ii <- as.list(maxid)
    ii[ipar] <- TRUE
    bf <- drop(do.call('[', c(bf2obj["logbf"], ii, list(drop = FALSE))))
    pdata <- data.frame(bf2obj[pars], logbf = bf)
    if ("type" %in% nmdots) {
      plot(pdata, ...)
    } else {
      plot(pdata, type = "l", ...)
    }
  } else if (ptype == "contour") {
    maxid <- bf2obj$indmax
    ii <- as.list(maxid)
    ii[ipar] <- TRUE
    bf <- do.call('[', c(list(bf2obj$logbf), ii, list(drop = FALSE)))
    jj <- seq(N); jj[sort(ipar)] <- ipar
    bf <- drop(aperm(bf, jj))
    CC <- call("contour", bf2obj[[pars[1]]], bf2obj[[pars[2]]], bf,
               quote(...))
    if (!("xlab" %in% nmdots)) CC["xlab"] = pars[1]
    if (!("ylab" %in% nmdots)) CC["ylab"] = pars[2]
    eval(CC)
    points((bf2obj[[pars[1]]])[maxid[ipar[1]]],
           (bf2obj[[pars[2]]])[maxid[ipar[2]]])
  } else if (ptype == "profile") {
    ii <- as.list(rep(TRUE, N))
    if (prod(par("mfrow")) < npars) {
      oldpar <- par(mfrow = c(1, npars))
      on.exit(par(oldpar), add = TRUE)
    }
    for (i in 1:npars) {
      bf <- apply(bf2obj[["logbf"]], ipar[i], max)
      pdata <- data.frame(bf2obj[pars[i]], logbf = bf)
      CC <- call("plot", pdata, quote(...))
      if (!("type" %in% nmdots)) CC["type"] = "l"
      if (!("xlab" %in% nmdots)) CC["xlab"] = ""
      eval(CC)
      if (!("xlab" %in% nmdots)) {
        title(xlab = switch(pars[i],
                linkp = expression(nu),
                phi = expression(phi),
                omg = expression(omega),
                kappa = expression(kappa)))
      }
    }
  }
  invisible()
}


##' Estimation by empirical Bayes.
##'
##' This function is a wrap around \code{\link{bf2new}} using the
##' "L-BFGS-B" method of the function \code{\link[stats]{optim}} to
##' estimate the parameters.
##' @title Empirical Bayes estimator
##' @param bf1obj Output from the function \code{\link{bf1skel}} which
##' contains the Bayes factors and importance sampling weights.
##' @param paroptim A named list with the components "linkp", "phi",
##' "omg", "kappa". Each component must be numeric with length 1, 2,
##' or 3 with elements in increasing order but for the binomial family
##' linkp is also allowed to be the character "logit" and "probit". If
##' the compontent's length is 1, then the corresponding parameter is
##' considered to be fixed at that value. If 2, then the two numbers
##' denote the lower and upper bounds for the optimisation of that
##' parameter (infinities are allowed). If 3, these correspond to
##' lower bound, starting value, upper bound for the estimation of
##' that parameter.
##' @param useCV Whether to use control variates for finer
##' corrections.
##' @param control A list of control parameters for the optimisation.
##' See \code{\link[stats]{optim}}.
##' @return The output from the function \code{\link[stats]{optim}}.
##'   The \code{"value"} element is the log-Bayes factor, not the
##'   negative log-Bayes factor.
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
##' est <- bf2optim(bf, list(linkp = linkp, phi = c(100, 200), omg = c(0, 2)))
##' est
##' }
##' @importFrom stats optim
##' @useDynLib geoBayes calcb_no_st calcb_mu_st calcb_wo_st calcb_tr_st
##' @useDynLib geoBayes calcb_no_cv calcb_mu_cv calcb_wo_cv calcb_tr_cv
##' @useDynLib geoBayes calcbd_no calcbd_mu calcbd_wo calcbd_tr
##' @export
bf2optim <- function (bf1obj, paroptim, useCV = TRUE,
                      control = list()) {

  ## Logical input
  useCV <- isTRUE(useCV)

  ## Extract model variables
  transf <- bf1obj$transf
  real_transf <- bf1obj$real_transf
  itr <- bf1obj$itr
  y <- bf1obj$response
  n <- length(y)
  l <- bf1obj$weights
  F <- bf1obj$modelmatrix
  p <- NCOL(F)
  family <- bf1obj$family
  ifam <- .geoBayes_family(family)
  corrfcn <- bf1obj$corrfcn
  icf <- .geoBayes_correlation(corrfcn)
  betm0 <- bf1obj$betm0
  betQ0 <- bf1obj$betQ0
  ssqdf <- bf1obj$ssqdf
  ssqsc <- bf1obj$ssqsc
  dispersion <- bf1obj$dispersion
  tsqdf <- bf1obj$tsqdf
  tsqsc <- bf1obj$tsqsc
  tsq <- if (ifam == 0) tsqsc else dispersion
  loc <- bf1obj$locations
  dm <- sp::spDists(loc)
  sample <- bf1obj$sample2
  Ntot <- NCOL(sample)
  isweights <- bf1obj$isweights
  zcv <- bf1obj$controlvar
  kg <- NCOL(zcv)


  ## Check paroptim
  ## parnmall <- c("linkp", "phi", "omg", "kappa")
  paroptim <- getparoptim(paroptim, ifam, icf)
  pstart <- as.double(paroptim$pstart)
  lower <- paroptim$lower
  upper <- paroptim$upper
  estim <- paroptim$estim

  ## Check for non-finite values in control variates
  if (useCV && any(!is.finite(zcv))) {
    warning ("Computed non-finite control variates, probably due to numerical
overflow. Control variates corrections will not be used.")
    useCV <- FALSE
  }

  froutine <- paste0("calcb_", real_transf)
  groutine <- paste0("calcbd_", real_transf)
  if (real_transf == "wo" || real_transf == "tr") {
    ifam <- -ifam
  }

  ## Function to optimise
  if (useCV) {
    froutine <- paste(froutine, "_cv", sep = "")
    QRin <- QRintercept(zcv)
    fn <- function (par) {
      parin <- pstart
      parin[estim] <- par
      parin[2:4] <- exp(parin[2:4])
      RUN <- .Fortran(froutine, 0,
                      parin[2], parin[1], parin[3], parin[4],
                      as.integer(icf), 1L, 1L, as.integer(Ntot),
                      as.double(sample), as.double(isweights),
                      as.double(QRin),
                      as.integer(n), as.integer(p),
                      as.double(betm0),
                      as.double(betQ0), as.double(ssqdf), as.double(ssqsc),
                      max(tsqdf, 0), as.double(tsq), as.double(y),
                      as.double(l), as.double(F), as.double(dm),
                      as.integer(ifam), as.integer(itr),
                      PACKAGE = "geoBayes")
      -RUN[[1]][1]
    }
  } else {
    QRin <- numeric(Ntot)
    froutine <- paste(froutine, "_st", sep = "")
    fn <- function (par) {
      parin <- pstart
      parin[estim] <- par
      parin[2:4] <- exp(parin[2:4])
      RUN <- .Fortran(froutine, 0,
                      parin[2], parin[1], parin[3], parin[4],
                      as.integer(icf), 1L, 1L,
                      as.integer(Ntot), as.double(sample),
                      as.double(isweights), as.integer(n), as.integer(p),
                      as.double(betm0),
                      as.double(betQ0), as.double(ssqdf), as.double(ssqsc),
                      max(tsqdf, 0), as.double(tsq), as.double(y), as.double(l),
                      as.double(F), as.double(dm),
                      as.integer(ifam), as.integer(itr),
                      PACKAGE = "geoBayes")
      -RUN[[1]][1]
    }
  }

### Compute derivative
  gr <- function (par) {
    parin <- pstart
    parin[estim] <- par
    parin[2:4] <- exp(parin[2:4])
    RUN <- .Fortran(groutine, 0, 0, 0, 0, 0,
                    parin[2], parin[1], parin[3], parin[4],
                    as.integer(icf), as.integer(Ntot), as.double(sample),
                    as.double(isweights),
                    as.double(QRin),
                    as.integer(n), as.integer(p),
                    as.double(betm0),
                    as.double(betQ0), as.double(ssqdf), as.double(ssqsc),
                    max(tsqdf, 0), as.double(tsq), as.double(y), as.double(l),
                    as.double(F), as.double(dm),
                    as.integer(ifam), as.integer(itr), as.integer(useCV),
                    PACKAGE = "geoBayes")
    gg <- -do.call("c", RUN[2:5])*c(1, parin[2:4])
    gg[estim]
  }


  logbf <- bf1obj$logbf
  imaxlogbf <- which.max(logbf)
  wbf <- exp(logbf - logbf[imaxlogbf] - log(sum(exp(logbf - logbf[imaxlogbf]))))
  pstart.d <- colSums(data.frame(bf1obj$pnts[c("nu", "phi", "omg", "kappa")])*
                      wbf)
  i <- is.na(pstart) & estim
  pstart[i] <- pmax(pmin(upper[i], pstart.d[i]), lower[i])
  scale <- control$parscale; if (is.null(scale)) scale <- 1
  pstart[2:4] <- log(pstart[2:4])
  lower[2:4] <- log(lower[2:4])
  upper[2:4] <- log(upper[2:4])
  method <- "L-BFGS-B"
  if (isTRUE(control$fnscale < 0)) control$fnscale <- -control$fnscale
###   op <- stats::nlminb(pstart[estim], fn, gr, scale = scale,
###                       control = control,
###                       lower = lower[estim], upper = upper[estim])
###   op$value <- -op$objective + bf1obj$logbf[1]
###   op$hessian <- optimHess(op$par, fn, gr, control = control)
###   op <- optimx::optimx(pstart[estim], fn, gr, method = method,
###                        lower = lower[estim], upper = upper[estim],
###                        control = control, hessian = TRUE)
  op <- stats::optim(pstart[estim], fn, gr, method = method,
                     lower = lower[estim], upper = upper[estim],
                     control = control, hessian = TRUE)
  op$value <- -op$value + bf1obj$logbf[1]
  parout <- pstart
  parout[estim] <- op$par
  names(parout) <- c("linkp", "phi", "omg", "kappa")
  op$par <- c(parout[1], exp(parout[2:4]))
  op$hessian <- op$hessian*tcrossprod(c(1, op$par[2:4])[estim])
  op$fixed <- !estim
  op
}


QRintercept <- function(zcv) {
  ## Retrun the column from Q-R that corresponds to the intercept.
  QR <- qr(zcv)
  Rcv <- qr.R(QR)
  Qcv <- qr.Q(QR)
###   tol <- max(abs(diag(Rcv))) * NROW(zcv) * .Machine$double.eps
###   ii <- abs(diag(Rcv)) > tol
###   zcvrnk <- sum(ii)
  zcvii <- QR$pivot == 1
##   zcvQR <- Qcv
##   zcvQR <- .Fortran("dtrsm", "r", "u", "t", "n", as.integer(Ntot),
##                     as.integer(zcvrnk), 1.0, as.double(Rcv),
##                     as.integer(kg), as.double(zcvQR), as.integer(Ntot))[[10]]
  ## zcvQR <- t(backsolve(Rcv,t(Qcv)))
  backsolve(Rcv,t(Qcv))[zcvii, ]
}
