##' Draw MCMC samples from the Spatial GLMM with known link function
##'
##' The four-parameter prior for \code{phi} is defined by
##' \deqn{\propto (\phi - \theta_4)^{\theta_2 -1} \exp\{-(\frac{\phi -
##' \theta_4}{\theta_1})^{\theta_3}\}}{propto (phi -
##' phiprior[4])^(phiprior[2]-1) *
##' exp(-((phi-phiprior[4])/phiprior[1])^phiprior[3])} for \eqn{\phi >
##' \theta_4}{phi > phiprior[4]}. The prior for \code{omg} is similar.
##' The prior parameters correspond to scale, shape, exponent, and
##' location. See \code{arXiv:1005.3274} for details of this
##' distribution.
##'
##' The GEV (Generalised Extreme Value) link is defined by \deqn{\mu =
##' 1 - \exp\{-\max(0, 1 + \nu x)^{\frac{1}{\nu}}\}}{mu = 1 -
##' \exp[-max(0, 1 + nu x)^(1/nu)]} for any real \eqn{\nu}{nu}. At
##' \eqn{\nu = 0}{nu = 0} it reduces to the complementary log-log
##' link.
##' @title MCMC samples from the Spatial GLMM
##' @param formula A representation of the model in the form
##'   \code{response ~ terms}. The response must be set to \code{NA}'s
##'   at the prediction locations (see the examples on how to do this
##'   using the function \code{\link{stackdata}}). At the observed
##'   locations the response is assumed to be a total of replicated
##'   measurements. The number of replications is inputted using the
##'   argument \code{weights}.
##' @param family The distribution of the data. The
##'   \code{"GEVbinomial"} family is the binomial family with link the
##'   GEV link (see Details).
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
##' @param linkp Parameter of the link function. A scalar value.
##' @param phi Optional starting value for the MCMC for the
##'   spatial range parameter \code{phi}. Defaults to the mean of its
##'   prior. If \code{corrtuning[["phi"]]} is 0, then this argument is required and
##'   it corresponds to the fixed value of \code{phi}. This can be a
##'   vector of the same length as Nout.
##' @param omg Optional starting value for the MCMC for the
##'   relative nugget parameter \code{omg}. Defaults to the mean of
##'   its prior. If \code{corrtuning[["omg"]]} is 0, then this argument is required
##'   and it corresponds to the fixed value of \code{omg}. This can be
##'   a vector of the same length as Nout.
##' @param kappa Optional starting value for the MCMC for the
##'   spatial correlation parameter \code{kappa} (Matern smoothness or
##'   exponential power). Defaults to the mean of
##'   its prior. If \code{corrtuning[["kappa"]]} is 0 and it is needed for
##'   the chosen correlation function, then this argument is required
##'   and it corresponds to the fixed value of \code{kappa}. This can be
##'   a vector of the same length as Nout.
##' @param Nout Number of MCMC samples to return. This can be a vector
##'   for running independent chains.
##' @param Nthin The thinning of the MCMC algorithm.
##' @param Nbi The burn-in of the MCMC algorithm.
##' @param betm0 Prior mean for beta (a vector or scalar).
##' @param betQ0 Prior standardised precision (inverse variance)
##'   matrix. Can be a scalar, vector or matrix. The first two imply a
##'   diagonal with those elements. Set this to 0 to indicate a flat
##'   improper prior.
##' @param ssqdf Degrees of freedom for the scaled inverse chi-square
##'   prior for the partial sill parameter.
##' @param ssqsc Scale for the scaled inverse chi-square prior for the
##'   partial sill parameter.
##' @param corrpriors A list with the components \code{phi},
##'   \code{omg} and \code{kappa} as needed. These correspond to the
##'   prior distribution parameters. For \code{phi} and \code{omg} it
##'   must be a vector of length 4. The generalized inverse gamma
##'   prior is assumed and the input corresponds to the parameters
##'   scale, shape, exponent, location in that order (see Details).
##'   For \code{kappa} it must be a vector of length 2. A uniform
##'   prior is assumed and the input corresponds to the lower and
##'   upper bounds in that order.
##' @param corrtuning A vector or list with the components \code{phi},
##'   \code{omg} and \code{kappa} as needed. These correspond to the
##'   random walk parameter for the Metropolis-Hastings step. Smaller values
##'   increase the acceptance ratio. Set this to 0 for fixed
##'   parameter value.
##' @param dispersion The fixed dispersion parameter.
##' @param longlat How to compute the distance between locations. If
##'   \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##'   distance. See \code{\link[sp]{spDists}}.
##' @param test Whether this is a trial run to monitor the acceptance
##'   ratio of the random walk for \code{phi} and \code{omg}. If set
##'   to \code{TRUE}, the acceptance ratio will be printed on the
##'   screen every 100 iterations of the MCMC. Tune the \code{phisc}
##'   and \code{omgsc} parameters in order to achive 20 to 30\%
##'   acceptance. Set this to a positive number to change the default
##'   100. No thinning or burn-in are done when testing.
##' @return A list containing the objects \code{MODEL}, \code{DATA},
##'   \code{FIXED}, \code{MCMC} and \code{call}. The MCMC samples are
##'   stored in the object \code{MCMC} as follows:
##' \itemize{
##'  \item \code{z} A matrix containing the MCMC samples for the
##' spatial random field. Each column is one sample.
##'  \item \code{mu} A matrix containing the MCMC samples for the
##' mean response (a transformation of z). Each column is one sample.
##'  \item \code{beta} A matrix containing the MCMC samples for the
##' regressor coefficients. Each column is one sample.
##'  \item \code{ssq} A vector with the MCMC samples for the partial
## sill parameter.
##'  \item \code{phi} A vector with the MCMC samples for the spatial
##' range parameter, if sampled.
##'  \item \code{omg} A vector with the MCMC samples for the relative
##' nugget parameter, if sampled.
##'  \item \code{logLik} A vector containing the value of the
##' log-likelihood evaluated at each sample.
##'  \item \code{acc_ratio} The acceptance ratio for the joint update
##' of the parameters \code{phi} and \code{omg}, if sampled.
##'  \item \code{sys_time} The total computing time for the MCMC sampling.
##'  \item \code{Nout}, \code{Nbi},  \code{Nthin} As in input. Used
##' internally in other functions.
##' }
##' The other objects contain input variables. The object \code{call}
##'   contains the function call.
##' @examples \dontrun{
##' data(rhizoctonia)
##'
##' ### Create prediction grid
##' predgrid <- mkpredgrid2d(rhizoctonia[c("Xcoord", "Ycoord")],
##'                          par.x = 100, chull = TRUE, exf = 1.2)
##'
##' ### Combine observed and prediction locations
##' rhizdata <- stackdata(rhizoctonia, predgrid$grid)
##' ##'
##' ### Define the model
##' corrf <- "spherical"
##' family <- "binomial.probit"
##' kappa <- 0
##' ssqdf <- 1
##' ssqsc <- 1
##' betm0 <- 0
##' betQ0 <- .01
##' phiprior <- c(100, 1, 1000, 100) # U(100, 200)
##' phisc <- 3
##' omgprior <- c(2, 1, 1, 0)        # Exp(mean = 2)
##' omgsc <- .1
##' ##'
##' ### MCMC sizes
##' Nout <- 100
##' Nthin <- 1
##' Nbi <- 0
##'
##' ### Trial run
##' emt <- mcsglmm(Infected ~ 1, family, rhizdata, weights = Total,
##'                atsample = ~ Xcoord + Ycoord,
##'                Nout = Nout, Nthin = Nthin, Nbi = Nbi,
##'                betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
##'                corrpriors = list(phi = phiprior, omg = omgprior), 
##'                corrfcn = corrf, kappa = kappa,
##'                corrtuning = list(phi = phisc, omg = omgsc, kappa = 0),
##'                dispersion = 1, test = 10)
##'
##' ### Full run
##' emc <- update(emt, test = FALSE)
##'
##' emcmc <- mcmcmake(emc)
##' summary(emcmc[, c("phi", "omg", "beta", "ssq")])
##' plot(emcmc[, c("phi", "omg", "beta", "ssq")])
##' }
##' @importFrom sp spDists
##' @importFrom stats model.matrix model.response model.weights
##'   as.formula update model.offset
##' @useDynLib geoBayes mcspsamtry mcspsample
##' @export
mcsglmm <- function (formula, family = "gaussian",
                     data, weights, subset, offset, 
                     atsample, corrfcn = "matern",
                     linkp, phi, omg, kappa,
                     Nout, Nthin = 1, Nbi = 0, betm0, betQ0, ssqdf, ssqsc,
                     corrpriors, corrtuning,
                     dispersion = 1, longlat = FALSE, test = FALSE) {
  cl <- match.call()
  ## Family
  ifam <- .geoBayes_family(family)
  if (ifam) {
    family <- .geoBayes_models$family[ifam]
  } else {
    stop ("This family has not been implemented.")
  }
  if (.geoBayes_models$needlinkp[ifam]) {
    if (missing(linkp))
      stop ("Missing input linkp.")
  } else {
    linkp <- 0
  }

  ## Correlation function
  icf <- .geoBayes_correlation(corrfcn)
  corrfcn <- .geoBayes_corrfcn$corrfcn[icf]
  needkappa <- .geoBayes_corrfcn$needkappa[icf]

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
  atsample <- update(atsample, NULL ~ . + 0) # No response and no intercept
  mfatc <- mfc
  mfatc$weights = NULL
  mfatc$formula = atsample
  mfat <- eval(mfatc, parent.frame())
  loc <- as.matrix(mfat)
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations")
  if (corrfcn == "spherical" && NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }

  ## Split sample, prediction
  ii <- is.finite(yy)
  y <- yy[ii]
  k <- sum(ii)
  l <- ll[ii]
  l <- if (is.null(l)) rep.int(1.0, k) else as.double(l)
  if (any(!is.finite(l))) stop ("Non-finite values in the weights")
  if (any(l <= 0)) stop ("Non-positive weights not allowed")
  if (grepl("^binomial(\\..+)?$", family)) {
    l <- l - y # Number of failures
  }
  F <- FF[ii, , drop = FALSE]
  offset <- oofset[ii]
  dm <- sp::spDists(loc[ii, , drop = FALSE], longlat = longlat)
  k0 <- sum(!ii)
  if (k0 > 0) {
    F0 <- FF[!ii, , drop = FALSE]
    dmdm0 <- sp::spDists(loc[ii, , drop = FALSE], loc[!ii, , drop = FALSE],
                         longlat = longlat)
    offset0 <- oofset[!ii]
  } else {
    F0 <- dmdm0 <- offset0 <- numeric(0)
    dim(F0) <- c(0, p)
    dim(dmdm0) <- c(k, 0)
  }

  ## Prior for ssq
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## Prior for beta
  betaprior <- getbetaprior(betm0, betQ0, p)
  betm0 <- betaprior$betm0
  betQ0 <- betaprior$betQ0

  ## Other fixed parameters
  dispersion <- as.double(dispersion)
  if (dispersion <= 0) stop ("Invalid argument dispersion")

  nu <- .geoBayes_getlinkp(linkp, ifam)

  ## MCMC samples
  Nout <- as.integer(Nout)
  if (any(Nout < 0)) stop ("Negative MCMC sample size entered.")
  nch <- length(Nout)         # Number of chains
  Nmc <- Nout                 # Size of each chain
  Nout <- sum(Nout)           # Total MCMC size
  Nbi <- as.integer(Nbi)
  Nthin <- as.integer(Nthin)
  lglk <- numeric(Nout)
  z <- matrix(0, k, Nout)
  z0 <- matrix(0, k0, Nout)
  beta <- matrix(0, p, Nout)
  ssq <- numeric(Nout)

  ## Starting values for correlation parameters
  phisc <- corrtuning[["phi"]]
  if (is.null(phisc) || !is.numeric(phisc) || phisc < 0)
    stop ("Invalid tuning parameter for phi.")
  if (phisc > 0) {
    phipars <- check_gengamma_prior(corrpriors[["phi"]])
  } else phipars <- rep.int(0, 4)
  if (missing(phi)) {
    if (phisc == 0) {
      stop ("Argument phi needed for fixed phi")
    } else {
      if(phipars[2] == -1) {
        tmp <- .1/abs(phipars[3])
      } else {
        tmp <- abs((phipars[2]+1)/phipars[3])
      }
      phistart <- phipars[4] + phipars[1]*gamma(tmp)/
        gamma(phipars[2]/phipars[3])
    }
  } else {
    phistart <- as.double(phi)
    if (phisc > 0 && phistart <= phipars[4]) {
      stop ("Starting value for phi not in the support of its prior")
    }
  }
  phi <- numeric(Nout)
  phi[cumsum(c(1, Nmc[-nch]))] <- phistart

  omgsc <- corrtuning[["omg"]]
  if (is.null(omgsc) || !is.numeric(omgsc) || omgsc < 0)
    stop ("Invalid tuning parameter for omg.")
  if (omgsc > 0) {
    omgpars <- check_gengamma_prior(corrpriors[["omg"]])
  } else omgpars <- rep.int(0, 4)
  if (missing(omg)) {
    if (omgsc == 0) {
      stop ("Argument omg needed for fixed omg")
    } else {
      if(omgpars[2] == -1) {
        tmp <- .1/abs(omgpars[3])
      } else {
        tmp <- abs((omgpars[2]+1)/omgpars[3])
      }
      omgstart <- omgpars[4] + omgpars[1]*gamma(tmp)/
        gamma(omgpars[2]/omgpars[3])
    }
  } else {
    omgstart <- as.double(omg)
    if (omgsc > 0 && omgstart <= omgpars[4]) {
      stop ("Starting value for omg not in the support of its prior")
    }
  }
  omg <- numeric(Nout)
  omg[cumsum(c(1, Nmc[-nch]))] <- omgstart

  if (needkappa) {
    kappasc <- corrtuning[["kappa"]]
  } else {
    kappasc <- 0
    kappa <- 0
  }
  if (is.null(kappasc) || !is.numeric(kappasc) || kappasc < 0)
    stop ("Invalid tuning parameter for kappa.")
  if (kappasc > 0) {
    kappapars <- check_unif_prior(corrpriors[["kappa"]])
  } else kappapars <- c(0, 0)
  if (missing(kappa)) {
    if (kappasc == 0) {
      stop ("Argument kappa needed for fixed kappa")
    } else {
      kappastart <- (kappapars[1] + kappapars[2])*.5
    }
  } else {
    kappastart <- as.double(kappa)
  }
  if (kappasc > 0) {
    kappastart <- .geoBayes_getkappa(kappastart, icf)
    kappapars <- .geoBayes_getkappa(kappapars, icf)
    if (kappastart >= kappapars[2] || kappastart <= kappapars[1]) {
      stop ("Starting value for kappa not in the support of its prior")
    }
  }
  kappa <- numeric(Nout)
  kappa[cumsum(c(1, Nmc[-nch]))] <- kappastart

  ## Run code
  if (test > 0) { # Running a test
    if (is.logical(test)) test <- 100
    test <- as.integer(test)
    acc <- 0L
    tm <- system.time({
      RUN <- .Fortran("mcspsamtry", ll = lglk, z = z, phi = phi, omg = omg,
                      kappa = kappa,
                      acc = acc,
                      as.double(y), as.double(l), as.double(F),
                      as.double(offset), 
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(phipars),
                      as.double(omgpars), as.double(kappapars),
                      as.double(phisc), as.double(omgsc), as.double(kappasc),
                      as.integer(icf),
                      as.double(nu), as.double(dispersion), as.double(dm),
                      as.integer(Nout), as.integer(test), as.integer(k),
                      as.integer(p), as.integer(ifam),
                      PACKAGE = "geoBayes")
    })
    ## Store samples
    ll <- RUN$ll
    zz0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN$z
    mm0 <- NULL
    beta <- NULL
    ssq <- NULL
    phi <- RUN$phi
###     attr(phi, 'fixed') <- phisc == 0
    omg <- RUN$omg
###     attr(omg, 'fixed') <- omgsc == 0
###     attr(nu, 'fixed') <- TRUE
    kappa <- RUN$kappa
    acc_ratio <- RUN$acc/Nout
    Nthin <- 1
    Nbi <- 0
###     out <- list(z = zz0, beta = beta, ssq = ssq, phi = phi, omg = omg, nu = nu,
###                 logLik = ll, acc_ratio = acc_ratio, sys_time = tm,
###                 Nout = Nout, Nbi = Nbi, Nthin = Nthin,
###                 response = y, weights = l, modelmatrix = F, family = family,
###                 betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
###                 corrfcn = corrfcn, kappa = kappa,
###                 dispersion = dispersion, locations = loc[ii, , drop = FALSE],
###                 longlat = longlat, whichobs = ii)
  } else {
    acc <- integer(nch)
    tm <- system.time({
      RUN <- .Fortran("mcspsample", ll = lglk, z = z, z0 = z0,
                      mu = z, mu0 = z0,
                      beta = beta, ssq = ssq,
                      phi = phi, omg = omg, kappa = kappa, acc = acc,
                      as.double(y), as.double(l), as.double(F),
                      as.double(offset), as.double(F0), as.double(offset0), 
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(phipars), as.double(omgpars),
                      as.double(kappapars),
                      as.double(phisc), as.double(omgsc), as.double(kappasc),
                      as.integer(icf),
                      as.double(nu), as.double(dispersion), as.double(dm),
                      as.double(dmdm0), as.integer(nch), as.integer(Nmc),
                      as.integer(Nout), as.integer(Nbi),
                      as.integer(Nthin), as.integer(k), as.integer(k0),
                      as.integer(p), as.integer(ifam),
                      PACKAGE = "geoBayes")
    })
    ## Store samples
    ll <- RUN$ll
    zz0 <- mm0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN$z
    zz0[!ii, ] <- RUN$z0
    mm0[ii, ] <- RUN$mu
    mm0[!ii, ] <- RUN$mu0
    beta <- RUN$beta
    ssq <- RUN$ssq
    phi <- RUN$phi
###    attr(phi, 'fixed') <- phisc == 0
    omg <- RUN$omg
###    attr(omg, 'fixed') <- omgsc == 0
###    attr(nu, 'fixed') <- TRUE
    kappa <- RUN$kappa
    acc_ratio <- RUN$acc/(Nmc*Nthin + max(Nthin, Nbi))
###     out <- list(z = zz0, mu = mm0,
###                 beta = beta, ssq = ssq, phi = phi, omg = omg, nu = nu,
###                 logLik = ll, acc_ratio = acc_ratio, sys_time = tm,
###                 Nout = Nout, Nbi = Nbi, Nthin = Nthin,
###                 response = y, weights = l, modelmatrix = F, family = family,
###                 betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
###                 corrfcn = corrfcn, kappa = kappa,
###                 dispersion = dispersion, locations = loc[ii, , drop = FALSE],
###                 longlat = longlat, whichobs = ii)
  }
  MCMC <- FIXED <- MODEL <- DATA <- list()
  MCMC$z <- zz0
  MCMC$mu <- mm0
  MCMC$beta <- beta
  MCMC$ssq <- ssq
  FIXED$linkp <- as.vector(linkp)
  FIXED$linkp_num <- nu
  if (phisc == 0) {
    FIXED$phi <- phi[1]
  } else {
    MCMC$phi <- phi
  }
  if (omgsc == 0) {
    FIXED$omg <- omg[1]
  } else {
    MCMC$omg <- omg
  }
  if (kappasc == 0) {
    FIXED$kappa <- kappa[1]
  } else {
    MCMC$kappa <- kappa
  }
  MCMC$logLik <- ll
  MCMC$acc_ratio <- acc_ratio
  MCMC$sys_time <- tm
  MCMC$Nout <- Nout
  MCMC$Nbi <- Nbi
  MCMC$Nthin <- Nthin
  MCMC$whichobs <- ii
  DATA$response <- y
  DATA$weights <- l
  DATA$modelmatrix <- F
  DATA$offset <- offset
  DATA$locations <- loc[ii, , drop = FALSE]
  DATA$longlat <- longlat
  MODEL$family <- family
  MODEL$corrfcn <- corrfcn
  MODEL$betm0 <- betm0
  MODEL$betQ0 <- betQ0
  MODEL$ssqdf <- ssqdf
  MODEL$ssqsc <- ssqsc
  MODEL$phipars <- phipars
  MODEL$omgpars <- omgpars
  MODEL$dispersion <- dispersion
  out <- list(MODEL = MODEL, DATA = DATA, FIXED = FIXED, MCMC = MCMC, call = cl)
  class(out) <- "geomcmc"
  out
}


##' Draw MCMC samples from the transformed Gaussian model with known
##' link function
##'
##' Simulates from the posterior distribution of this model.
##' @title MCMC samples from the transformed Gaussian model
##' @param formula A representation of the model in the form
##' \code{response ~ terms}. The response must be set to \code{NA}'s
##' at the prediction locations (see the example in
##' \code{\link{mcsglmm}} for how to do this using
##' \code{\link{stackdata}}). At the observed locations the response
##' is assumed to be a total of replicated measurements. The number of
##' replications is inputted using the argument \code{weights}.
##' @param data An optional data frame containing the variables in the
##' model.
##' @param weights An optional vector of weights. Number of replicated
##' samples.
##' @param subset An optional vector specifying a subset of
##' observations to be used in the fitting process.
##' @param offset See \code{\link[stats]{lm}}.
##' @param atsample A formula in the form \code{~ x1 + x2 + ... + xd}
##' with the coordinates of the sampled locations.
##' @param corrfcn Spatial correlation function. See
##' \code{\link{geoBayes_correlation}} for details.
##' @param linkp Parameter of the link function. A scalar value.
##' @param phi Optional starting value for the MCMC for the
##'   spatial range parameter \code{phi}. Defaults to the mean of its
##'   prior. If \code{corrtuning[["phi"]]} is 0, then this argument is required and
##'   it corresponds to the fixed value of \code{phi}. This can be a
##'   vector of the same length as Nout.
##' @param omg Optional starting value for the MCMC for the
##'   relative nugget parameter \code{omg}. Defaults to the mean of
##'   its prior. If \code{corrtuning[["omg"]]} is 0, then this argument is required
##'   and it corresponds to the fixed value of \code{omg}. This can be
##'   a vector of the same length as Nout.
##' @param kappa Optional starting value for the MCMC for the
##'   spatial correlation parameter \code{kappa} (Matern smoothness or
##'   exponential power). Defaults to the mean of
##'   its prior. If \code{corrtuning[["kappa"]]} is 0 and it is needed for
##'   the chosen correlation function, then this argument is required
##'   and it corresponds to the fixed value of \code{kappa}. This can be
##'   a vector of the same length as Nout.
##' @param Nout Number of MCMC samples to return. This can be a vector
##'   for running independent chains.
##' @param Nthin The thinning of the MCMC algorithm.
##' @param Nbi The burn-in of the MCMC algorithm.
##' @param betm0 Prior mean for beta (a vector or scalar).
##' @param betQ0 Prior standardised precision (inverse variance)
##' matrix. Can be a scalar, vector or matrix. The first two imply a
##' diagonal with those elements. Set this to 0 to indicate a flat
##' improper prior.
##' @param ssqdf Degrees of freedom for the scaled inverse chi-square
##' prior for the partial sill parameter.
##' @param ssqsc Scale for the scaled inverse chi-square prior for the
##' partial sill parameter.
##' @param tsqdf Degrees of freedom for the scaled inverse chi-square
##' prior for the measurement error parameter.
##' @param tsqsc Scale for the scaled inverse chi-square prior for the
##' measurement error parameter.
##' @param corrpriors A list with the components \code{phi},
##'   \code{omg} and \code{kappa} as needed. These correspond to the
##'   prior distribution parameters. For \code{phi} and \code{omg} it
##'   must be a vector of length 4. The generalized inverse gamma
##'   prior is assumed and the input corresponds to the parameters
##'   scale, shape, exponent, location in that order (see Details).
##'   For \code{kappa} it must be a vector of length 2. A uniform
##'   prior is assumed and the input corresponds to the lower and
##'   upper bounds in that order.
##' @param corrtuning A vector or list with the components \code{phi},
##'   \code{omg} and \code{kappa} as needed. These correspond to the
##'   random walk parameter for the Metropolis-Hastings step. Smaller values
##'   increase the acceptance ratio. Set this to 0 for fixed
##'   parameter value.
##' @param longlat How to compute the distance between locations. If
##' \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##' distance. See \code{\link[sp]{spDists}}.
##' @param test Whether this is a trial run to monitor the acceptance
##' ratio of the random walk for \code{phi} and \code{omg}. If set to
##' \code{TRUE}, the acceptance ratio will be printed on the screen
##' every 100 iterations of the MCMC. Tune the \code{phisc} and
##' \code{omgsc} parameters in order to achive 20 to 30\% acceptance.
##' Set this to a positive number to change the default 100. No
##' thinning or burn-in are done when testing.
##' @return A list containing the objects \code{MODEL}, \code{DATA},
##'   \code{FIXED}, \code{MCMC} and \code{call}. The MCMC samples are
##'   stored in the object \code{MCMC} as follows:
##' \itemize{
##'  \item \code{z} A matrix containing the MCMC samples for the
##' spatial random field. Each column is one sample.
##'  \item \code{mu} A matrix containing the MCMC samples for the
##' mean response (a transformation of z). Each column is one sample.
##'  \item \code{beta} A matrix containing the MCMC samples for the
##' regressor coefficients. Each column is one sample.
##'  \item \code{ssq} A vector with the MCMC samples for the partial
## sill parameter.
##'  \item \code{tsq} A vector with the MCMC samples for the
##' measurement error variance.
##'  \item \code{phi} A vector with the MCMC samples for the spatial
##' range parameter, if sampled.
##'  \item \code{omg} A vector with the MCMC samples for the relative
##' nugget parameter, if sampled.
##'  \item \code{logLik} A vector containing the value of the
##' log-likelihood evaluated at each sample.
##'  \item \code{acc_ratio} The acceptance ratio for the joint update
##' of the parameters \code{phi} and \code{omg}, if sampled.
##'  \item \code{sys_time} The total computing time for the MCMC sampling.
##'  \item \code{Nout}, \code{Nbi},  \code{Nthin} As in input. Used
##' internally in other functions.
##' }
##' The other objects contain input variables. The object \code{call}
##'   contains the function call.
##' @examples \dontrun{
##' ### Load the data
##' data(rhizoctonia)
##' rhiz <- na.omit(rhizoctonia)
##' rhiz$IR <- rhiz$Infected/rhiz$Total # Incidence rate of the
##'                               # rhizoctonia disease
##'
##' ### Define the model
##' corrf <- "spherical"
##' ssqdf <- 1
##' ssqsc <- 1
##' tsqdf <- 1
##' tsqsc <- 1
##' betm0 <- 0
##' betQ0 <- diag(.01, 2, 2)
##' phiprior <- c(200, 1, 1000, 100) # U(100, 300)
##' phisc <- 1
##' omgprior <- c(3, 1, 1000, 0) # U(0, 3)
##' omgsc <- 1
##' linkp <- 1
##'
##' ## MCMC parameters
##' Nout <- 100
##' Nbi <- 0
##' Nthin <- 1
##'
##' samplt <- mcstrga(Yield ~ IR, data = rhiz,
##'                   atsample = ~ Xcoord + Ycoord, corrf = corrf,
##'                   Nout = Nout, Nthin = Nthin,
##'                   Nbi = Nbi, betm0 = betm0, betQ0 = betQ0,
##'                   ssqdf = ssqdf, ssqsc = ssqsc,
##'                   tsqdf = tsqdf, tsqsc = tsqsc,
##'                   corrprior = list(phi = phiprior, omg = omgprior),
##'                   linkp = linkp,
##'                   corrtuning = list(phi = phisc, omg = omgsc, kappa = 0),
##'                   test=10)
##'
##' sample <- update(samplt, test = FALSE)
##' }
##' @importFrom sp spDists
##' @importFrom stats model.matrix model.response model.weights
##'   as.formula update model.offset
##' @useDynLib geoBayes trgasamtry trgasample
##' @export
mcstrga <- function (formula,
                     data, weights, subset, offset,
                     atsample, corrfcn = "matern",
                     linkp, phi, omg, kappa,
                     Nout, Nthin = 1, Nbi = 0, betm0, betQ0, ssqdf, ssqsc,
                     tsqdf, tsqsc,
                     corrpriors, corrtuning,
                     longlat = FALSE,
                     test = FALSE) {
  cl <- match.call()

  family <- "transformed.gaussian"

  ## Correlation function
  icf <- .geoBayes_correlation(corrfcn)
  corrfcn <- .geoBayes_corrfcn$corrfcn[icf]
  needkappa <- .geoBayes_corrfcn$needkappa[icf]

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
  atsample <- update(atsample, NULL ~ . + 0) # No response and no intercept
  mfatc <- mfc
  mfatc$weights = NULL
  mfatc$formula = atsample
  mfat <- eval(mfatc, parent.frame())
  loc <- as.matrix(mfat)
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations")
  if (corrfcn == "spherical" && NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }

  ## Split sample, prediction
  ii <- is.finite(yy)
  y <- yy[ii]
  k <- sum(ii)
  l <- ll[ii]
  l <- if (is.null(l)) rep.int(1.0, k) else as.double(l)
  if (any(!is.finite(l))) stop ("Non-finite values in the weights")
  if (any(l <= 0)) stop ("Non-positive weights not allowed")
  ybar <- y/l
  F <- FF[ii, , drop = FALSE]
  offset <- oofset[ii]
  dm <- sp::spDists(loc[ii, , drop = FALSE], longlat = longlat)
  k0 <- sum(!ii)
  if (k0 > 0) {
    F0 <- FF[!ii, , drop = FALSE]
    dmdm0 <- sp::spDists(loc[ii, , drop = FALSE], loc[!ii, , drop = FALSE],
                         longlat = longlat)
    offset0 <- oofset[!ii]
  } else {
    F0 <- dmdm0 <- offset0 <- numeric(0)
    dim(F0) <- c(0, p)
    dim(dmdm0) <- c(k, 0)
  }

  ## Prior for ssq
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## Prior for beta
  betaprior <- getbetaprior(betm0, betQ0, p)
  betm0 <- betaprior$betm0
  betQ0 <- betaprior$betQ0

  ## Prior for tsq
  tsqdf <- as.double(tsqdf)
  if (tsqdf <= 0) stop ("Argument tsqdf must > 0")
  tsqsc <- as.double(tsqsc)
  if (tsqsc <= 0) stop ("Argument tsqsc must > 0")

  if (missing(linkp))
    stop ("Missing input linkp.")
  nu <- .geoBayes_getlinkp(linkp, family)

  ## MCMC samples
  Nout <- as.integer(Nout)
  if (any(Nout < 0)) stop ("Negative MCMC sample size entered.")
  nch <- length(Nout)         # Number of chains
  Nmc <- Nout                 # Size of each chain
  Nout <- sum(Nout)           # Total MCMC size
  Nbi <- as.integer(Nbi)
  Nthin <- as.integer(Nthin)
  lglk <- numeric(Nout)
  z <- matrix(0, k, Nout)
  z0 <- matrix(0, k0, Nout)
  beta <- matrix(0, p, Nout)
  ssq <- tsq <- numeric(Nout)

  ## Starting values for correlation parameters
  phisc <- corrtuning[["phi"]]
  if (is.null(phisc) || !is.numeric(phisc) || phisc < 0)
    stop ("Invalid tuning parameter for phi.")
  if (phisc > 0) {
    phipars <- check_gengamma_prior(corrpriors[["phi"]])
  } else phipars <- rep.int(0, 4)
  if (missing(phi)) {
    if (phisc == 0) {
      stop ("Argument phi needed for fixed phi")
    } else {
      if(phipars[2] == -1) {
        tmp <- .1/abs(phipars[3])
      } else {
        tmp <- abs((phipars[2]+1)/phipars[3])
      }
      phistart <- phipars[4] + phipars[1]*gamma(tmp)/
        gamma(phipars[2]/phipars[3])
    }
  } else {
    phistart <- as.double(phi)
    if (phisc > 0 && phistart <= phipars[4]) {
      stop ("Starting value for phi not in the support of its prior")
    }
  }
  phi <- numeric(Nout)
  phi[cumsum(c(1, Nmc[-nch]))] <- phistart

  omgsc <- corrtuning[["omg"]]
  if (is.null(omgsc) || !is.numeric(omgsc) || omgsc < 0)
    stop ("Invalid tuning parameter for omg.")
  if (omgsc > 0) {
    omgpars <- check_gengamma_prior(corrpriors[["omg"]])
  } else omgpars <- rep.int(0, 4)
  if (missing(omg)) {
    if (omgsc == 0) {
      stop ("Argument omg needed for fixed omg")
    } else {
      if(omgpars[2] == -1) {
        tmp <- .1/abs(omgpars[3])
      } else {
        tmp <- abs((omgpars[2]+1)/omgpars[3])
      }
      omgstart <- omgpars[4] + omgpars[1]*gamma(tmp)/
        gamma(omgpars[2]/omgpars[3])
    }
  } else {
    omgstart <- as.double(omg)
    if (omgsc > 0 && omgstart <= omgpars[4]) {
      stop ("Starting value for omg not in the support of its prior")
    }
  }
  omg <- numeric(Nout)
  omg[cumsum(c(1, Nmc[-nch]))] <- omgstart

  if (needkappa) {
    kappasc <- corrtuning[["kappa"]]
  } else {
    kappasc <- 0
    kappa <- 0
  }
  if (is.null(kappasc) || !is.numeric(kappasc) || kappasc < 0)
    stop ("Invalid tuning parameter for kappa.")
  if (kappasc > 0) {
    kappapars <- check_unif_prior(corrpriors[["kappa"]])
  } else kappapars <- c(0, 0)
  if (missing(kappa)) {
    if (kappasc == 0) {
      stop ("Argument kappa needed for fixed kappa")
    } else {
      kappastart <- (kappapars[1] + kappapars[2])*.5
    }
  } else {
    kappastart <- as.double(kappa)
  }
  if (kappasc > 0) {
    kappastart <- .geoBayes_getkappa(kappastart, icf)
    kappapars <- .geoBayes_getkappa(kappapars, icf)
    if (kappastart >= kappapars[2] || kappastart <= kappapars[1]) {
      stop ("Starting value for kappa not in the support of its prior")
    }
  }
  kappa <- numeric(Nout)
  kappa[cumsum(c(1, Nmc[-nch]))] <- kappastart

  ## Run code
  if (test > 0) { # Running a test
    if (is.logical(test)) test <- 100
    test <- as.integer(test)
    acc <- 0L
    tm <- system.time({
      RUN <- .Fortran("trgasamtry", ll = lglk, z = z, phi = phi, omg = omg,
                      kappa = kappa, acc = acc,
                      as.double(ybar), as.double(l), as.double(F),
                      as.double(offset), 
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(tsqdf), as.double(tsqsc),
                      as.double(phipars), as.double(omgpars),
                      as.double(kappapars),
                      as.double(phisc), as.double(omgsc),
                      as.double(kappasc), as.integer(icf),
                      as.double(nu), as.double(dm), as.integer(Nout),
                      as.integer(test), as.integer(k), as.integer(p),
                      PACKAGE = "geoBayes")
    })
    ## Store samples
    ll <- RUN$ll
    zz0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN$z
    mm0 <- NULL
    beta <- NULL
    ssq <- NULL
    phi <- RUN$phi
###     attr(phi, 'fixed') <- phisc == 0
    omg <- RUN$omg
###     attr(omg, 'fixed') <- omgsc == 0
###     attr(nu, 'fixed') <- TRUE
    kappa <- RUN$kappa
    acc_ratio <- RUN$acc/Nout
    Nthin <- 1
    Nbi <- 0
###     out <- list(z = zz0, beta = beta, ssq = ssq, phi = phi, omg = omg, nu = nu,
###                 logLik = ll, acc_ratio = acc_ratio, sys_time = tm,
###                 Nout = Nout, Nbi = Nbi, Nthin = Nthin,
###                 response = y, weights = l, modelmatrix = F, family = family,
###                 betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
###                 corrfcn = corrfcn, kappa = kappa,
###                 tsqdf = tsqdf, tsqsc = tsqsc,
###                 locations = loc[ii, , drop = FALSE],
###                 longlat = longlat, whichobs = ii)
  } else {
    acc <- integer(nch)
    tm <- system.time({
      RUN <- .Fortran("trgasample", ll = lglk, z = z, z0 = z0,
                      mu = z, mu0 = z0, beta = beta, ssq = ssq,
                      tsq = tsq, phi = phi, omg = omg,
                      kappa = kappa, acc = acc,
                      as.double(ybar),
                      as.double(l), as.double(F),
                      as.double(offset), as.double(F0), as.double(offset0), 
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(tsqdf), as.double(tsqsc),
                      as.double(phipars), as.double(omgpars),
                      as.double(kappapars),
                      as.double(phisc), as.double(omgsc), as.double(kappasc),
                      as.integer(icf),
                      as.double(nu), as.double(dm), as.double(dmdm0),
                      as.integer(nch), as.integer(Nmc),
                      as.integer(Nout), as.integer(Nbi), as.integer(Nthin),
                      as.integer(k), as.integer(k0), as.integer(p),
                      PACKAGE = "geoBayes")
    })
    ## Store samples
    ll <- RUN$ll
    zz0 <- mm0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN$z
    zz0[!ii, ] <- RUN$z0
    mm0[ii, ] <- RUN$mu
    mm0[!ii, ] <- RUN$mu0
    beta <- RUN$beta
    ssq <- RUN$ssq
    tsq <- RUN$tsq
    phi <- RUN$phi
###     attr(phi, 'fixed') <- phisc == 0
    omg <- RUN$omg
###     attr(omg, 'fixed') <- omgsc == 0
###     attr(nu, 'fixed') <- TRUE
    kappa <- RUN$kappa
    acc_ratio <- RUN$acc/(Nmc*Nthin + max(Nthin, Nbi))
###     out <- list(z = zz0, mu = mm0, beta = beta, ssq = ssq, tsq = tsq,
###                 phi = phi, omg = omg, nu = nu,
###                 logLik = ll, acc_ratio = acc_ratio, sys_time = tm,
###                 Nout = Nout, Nbi = Nbi, Nthin = Nthin,
###                 response = ybar, weights = l, modelmatrix = F, family = family,
###                 betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
###                 corrfcn = corrfcn, kappa = kappa,
###                 tsqdf = tsqdf, tsqsc = tsqsc,
###                 locations = loc[ii, , drop = FALSE],
###                 longlat = longlat, whichobs = ii)
  }
  MCMC <- FIXED <- MODEL <- DATA <- list()
  MCMC$z <- zz0
  MCMC$mu <- mm0
  MCMC$beta <- beta
  MCMC$ssq <- ssq
  MCMC$tsq <- tsq
  FIXED$linkp <- linkp
  FIXED$linkp_num <- nu
  if (phisc == 0) {
    FIXED$phi <- phi[1]
  } else {
    MCMC$phi <- phi
  }
  if (omgsc == 0) {
    FIXED$omg <- omg[1]
  } else {
    MCMC$omg <- omg
  }
  if (kappasc == 0) {
    FIXED$kappa <- kappa[1]
  } else {
    MCMC$kappa <- kappa
  }
  MCMC$logLik <- ll
  MCMC$acc_ratio <- acc_ratio
  MCMC$sys_time <- tm
  MCMC$Nout <- Nout
  MCMC$Nbi <- Nbi
  MCMC$Nthin <- Nthin
  MCMC$whichobs <- ii
  DATA$response <- ybar
  DATA$weights <- l
  DATA$modelmatrix <- F
  DATA$offset <- offset
  DATA$locations <- loc[ii, , drop = FALSE]
  DATA$longlat <- longlat
  MODEL$family <- family
  MODEL$corrfcn <- corrfcn
  MODEL$betm0 <- betm0
  MODEL$betQ0 <- betQ0
  MODEL$ssqdf <- ssqdf
  MODEL$ssqsc <- ssqsc
  MODEL$tsqdf <- tsqdf
  MODEL$tsqsc <- tsqsc
  MODEL$phipars <- phipars
  MODEL$omgpars <- omgpars
  out <- list(MODEL = MODEL, DATA = DATA, FIXED = FIXED, MCMC = MCMC, call = cl)
  class(out) <- "geomcmc"
  out
}


##' Convert to an \code{\link[coda]{mcmc}} object.
##'
##' This function takes as input the one or more output(s) from
##' function \code{\link{mcsglmm}} or \code{\link{mcstrga}} and
##' returns an \code{\link[coda]{mcmc}} object or an
##' \code{\link[coda]{mcmc.list}} object for coda. The function
##' requires the \code{coda} package to be installed.
##'
##' The spatial random field components are assigned the names
##' \code{z_*} where \code{*} is a number beginning at 1. Similarly,
##' the regressor coefficients are assigned the names \code{beta_*} if
##' not unique, or simply \code{beta} if there is only one regressor.
##' The names \code{ssq}, \code{tsq}, \code{phi}, \code{omg}
##' correspond to the partial sill, measurement error variance,
##' spatial range, and relative nugget parameters respectively.
##' @title Convert to an \code{\link[coda]{mcmc}} object
##' @param ... Output(s) from the functions mentioned in the Details.
##' @return An mcmc object.
##' @examples \dontrun{
##'  ### Load the data
##' data(rhizoctonia)
##' rhiz <- na.omit(rhizoctonia)
##' rhiz$IR <- rhiz$Infected/rhiz$Total # Incidence rate of the
##'                               # rhizoctonia disease
##'  ### Define the model
##' corrf <- "spherical"
##' ssqdf <- 1
##' ssqsc <- 1
##' tsqdf <- 1
##' tsqsc <- 1
##' betm0 <- 0
##' betQ0 <- diag(.01, 2, 2)
##' phiprior <- c(200, 1, 1000, 100) # U(100, 300)
##' phisc <- 1
##' omgprior <- c(3, 1, 1000, 0) # U(0, 3)
##' omgsc <- 1.3
##' linkp <- 1
##' ## MCMC parameters
##' Nout <- 100
##' Nbi <- 0
##' Nthin <- 1
##'  ### Run MCMC
##' sample <- mcstrga(Yield ~ IR, data = rhiz,
##'                   atsample = ~ Xcoord + Ycoord, corrf = corrf,
##'                   Nout = Nout, Nthin = Nthin,
##'                   Nbi = Nbi, betm0 = betm0, betQ0 = betQ0,
##'                   ssqdf = ssqdf, ssqsc = ssqsc,
##'                   tsqdf = tsqdf, tsqsc = tsqsc,
##'                   linkp = linkp,
##'                   corrprior = list(phi = phiprior, omg = omgprior), 
##'                   corrtuning = list(phi = phisc, omg = omgsc, kappa = 0),
##'                   test = FALSE)
##' mcsample <- mcmcmake(sample)
##' plot(mcsample[, c("phi", "omg", "beta_1", "beta_2", "ssq", "tsq")],
##'      density = FALSE)
##' summary(mcsample[, c("phi", "omg", "beta_1", "beta_2", "ssq", "tsq")])
##' }
##' @seealso Functions such as \code{\link[coda]{plot.mcmc}} and
##' \code{\link[coda]{summary.mcmc}} in the \code{coda} package. The
##' function \code{\link[base]{do.call}} can be used to pass arguments
##' stored in a list.
##' @importFrom coda mcmc mcmc.list
##' @export
mcmcmake <- function (...) {
  ### This function takes as input the output from function mcmcrun
  ### and returns an mcmc object or an mcmc.list object for coda. The
  ### function requires the coda package to be installed.
  vnm <- c('z','beta','ssq','tsq','phi','omg', 'kappa')
  input <- list(...)
  nruns <- length(input)
  mcl <- list(); length(mcl) <- nruns
  for (j in seq_len(nruns)) {
    chain <- input[[j]]$MCMC[vnm]
    names(chain) <- vnm
    if (is.matrix(chain$z)) chain$z <- t(chain$z)
    if(isTRUE(ncol(chain$z) > 1)) {
      colnames(chain$z) <- paste0('z_', 1:ncol(chain$z))
    } else {
      chain$z <- as.vector(chain$z)
    }
    if (is.matrix(chain$beta)) chain$beta <- t(chain$beta)
    if(isTRUE(ncol(chain$beta) > 1)) {
      colnames(chain$beta) <- paste0('beta_', 1:ncol(chain$beta))
    } else {
      chain$beta <- as.vector(chain$beta)
    }
    chain <- do.call("cbind", chain)
    thin <- max(1, input[[j]]$MCMC$Nthin)
    start <- max(thin, input[[j]]$MCMC$Nbi) + 1
    mcl[[j]] <- coda::mcmc(chain, start = start, thin = thin)
  }
  if (nruns == 1) {
    return (mcl[[1]])
  } else if (nruns > 1) {
    return (do.call(coda::mcmc.list, mcl))
  } else return (NULL)
}


##' Return subset of MCMC chain.
##'
##'
##' @title Subset MCMC chain
##' @param x Output from the functions \code{\link{mcsglmm}} or
##'   \code{\link{mcstrga}}.
##' @param subset Logical or integer vector.
##' @param ... Further arguments to be passed to or from other methods.
##' @return A similar object as \code{x} with the subsetted chain.
##' @export
subset.geomcmc <- function (x, subset, ...) {
  if (class(x) != "geomcmc") stop ("Wrong class of object x.")
  subset <- as.vector(subset)
  if (is.logical(subset)) subset <- subset & !is.na(subset)
  nm_vec <- c("ssq", "tsq", "phi", "omg", "kappa", "logLik")
  nm_mat <- c("z", "mu", "beta")
  out <- x
  for (nm in nm_vec) {
    out$MCMC[[nm]] <- x$MCMC[[nm]][subset]
  }
  for (nm in nm_mat) {
    out$MCMC[[nm]] <- x$MCMC[[nm]][, subset, drop = FALSE]
  }
  out$MCMC$Nout <- length(out$MCMC[[nm_vec[1]]])
  out
}
