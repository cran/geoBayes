######################################################################
##
### Commentary: Functions to simulate from the model.
##
######################################################################

##' Simulate from a variety of spatial models.
##'
##' The spatial Gaussian random field is simulated using the Cholesky
##' decomposition of the covariance matrix.
##'
##' The sample from a quasi distribution uses a hack which matches the
##' mean and the variance of the distribution. See the source code for
##' details.
##' @title Simulation from a spatial model
##' @param n The number of instances to simulate
##' @param formula A representation of the model in the form
##'   \code{response ~ terms}. The LHS can be omitted. If the LHS
##'   exists, it can be of the form \code{y}, \code{y|z}, or sums of
##'   terms at either side of the \code{|} to specify the names of the
##'   variables to include in the data frame.
##' @param family The distribution of the data to simulate from.
##' @param data An optional data frame containing the variables in the
##'   model.
##' @param weights An optional vector of weights. Number of replicated
##'   samples for Gaussian and gamma, number of trials for binomial,
##'   time length for Poisson.
##' @param subset An optional set of indices. Simulations will be
##'   provided for those locations only.
##' @param atsample A formula of the form \code{~ Xcoord + Ycoord}
##'   specifying the sampled locations.
##' @param beta A vector of the regressor coefficents to use.
##' @param linkp The link function parameter.
##' @param phi The spatial range parameter.
##' @param omg The relative nugget parameter.
##' @param kappa The spatial smoothness parameter.
##' @param ssq The partial sill parameter.
##' @param corrfcn The correlation function to use.
##' @param longlat How to compute the distance between locations. If
##'   \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##'   distance. See \code{\link[sp]{spDists}}.
##' @param dispersion The fixed dispersion parameter. When this is not
##'   1 and the sample is from a binomial or a Poisson distribution,
##'   an approximate sample is returned.
##' @param returnGRF Whether to return the simulate Gaussian random
##'   field as well.
##' @param warndisp Whether to warn when sampling from a quasi
##'   distribution. This is the case for binomial and Poisson when the
##'   dispersion is not one.
##' @return A data frame containing the predictors, sampling
##'   locations, optional weights, and samples.
##' @name rsglmm
##' @examples \dontrun{
##' n <- 100
##' beta <- c(-2, 1)
##' phi <- .2
##' omg <- .3
##' linkp <- 0
##' ssq <- 1
##' l <- rep(10, n)
##' corrf <- "matern"
##' kappa <- .5
##' family <- "poisson"
##' Xcoord <- runif(n)
##' Ycoord <- runif(n)
##' f <- Xcoord + Ycoord
##' formula <- y|z ~ f
##' mydata <- rsglmm(1, formula, family, weights = l,
##'                  atsample = ~ Xcoord + Ycoord, beta = beta, linkp = linkp,
##'                  phi = phi, omg = omg, kappa = kappa, ssq = ssq,
##'                  corrfcn = corrf, returnGRF = TRUE)
##' }
##' @export
##' @importFrom stats rnorm rbinom rpois rgamma rbeta model.frame
##'   update as.formula
##' @importFrom sp spDists
rsglmm <- function(n, formula,
                   family = "gaussian",
                   data, weights, subset, atsample,
                   beta, linkp, phi, omg, kappa, ssq,
                   corrfcn = "matern",
                   longlat = FALSE, dispersion = 1, returnGRF = FALSE,
                   warndisp = TRUE) {

  ## Family
  ifam <- .geoBayes_family(family)
  if (ifam) {
    family <- .geoBayes_models$family[ifam]
  } else {
    stop ("You cannot simulate transformed Gaussian using this function.")
  }

  ## Logical input
  returnGRF <- as.logical(returnGRF)
  if (length(returnGRF) != 1 || is.na(returnGRF))
    stop ("Input returnGRF must be a single logical value.")
  warndisp <- as.logical(warndisp)
  if (length(warndisp) != 1 || is.na(warndisp))
    stop ("Input warndisp must be a single logical value.")

  ## Check dispersion input
  dispersion <- as.double(dispersion)
  if (dispersion <= 0) stop ("Invalid argument dispersion")

  ## Check the link parameter
  nu <- .geoBayes_getlinkp(linkp, ifam)

  ## Response name
  formula <- as.formula(formula)
  if (length(formula) == 2) { # No LHS
    ## RHS <- formula[[2]]
    nmy <- "y"
    nmz <- "z"
  } else { # LHS names provided
    LHS <- formula[[2]]
    ## RHS <- formula[[3]]
    if ("|" %in% all.names(LHS)) { # LHS contains bar
      nmy <- all.vars(LHS[[2]])    # Before the bar
      nmz <- all.vars(LHS[[3]])    # After the bar
    } else {
      nmy <- all.vars(LHS)
      nmz <- "z"
    }
  }
  formula0 <- update(formula, NULL ~ .)

  ## Design matrix and data
  if (missing(data)) data <- environment(formula)
  mfc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"),
             names(mfc), 0L)
  mfc <- mfc[c(1L, m)]
  mfc$formula <- formula0
  mfc$drop.unused.levels <- TRUE
  mfc$na.action <- "na.fail"
  mfc[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfc, parent.frame())
  mt <- attr(mf, "terms")
  FF <- model.matrix(mt,mf)
  nloc <- nrow(FF)

  ## Weights
  weightsc <- mfc$weights
  noweights <- is.null(weightsc)
  if (noweights) {
    ll <- rep.int(1, nloc)
    nmweights <- weights <- NULL
  } else {
    ll <- model.weights(mf)
    if (any(!is.finite(ll))) {
      stop ("Non-finite values in the weights")
    } else if (any(ll <= 0)) {
      stop ("Non-positive weights not allowed")
    }
    weights <- eval(weightsc, envir = data)
    if (length(weightsc) > 1) {
      nmweights <- "(weights)"
    } else {
      nmweights <- as.character(weightsc)
    }
    weights <- list(weights)
    names(weights) <- nmweights
  }

  ## All locations
  atsample <- update(atsample, NULL ~ . + 0) # No response and no intercept
  mfatc <- mfc
  mfatc$weights = NULL
  mfatc$formula = atsample
  mfat <- eval(mfatc, parent.frame())
  loc <- as.matrix(mfat)
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations")

  ## Simulate
  z <- randgrf(n, loc, beta, phi, omg, kappa, ssq, FF, corrfcn, longlat)
  pary <- linkinv(z, linkp, family)
  if (warndisp && (dispersion != 1) &&
      grepl("^(binomial|poisson)(\\..+)?$", family))
    warning ("Sampling from a quasi distribution: results are approximate.")
  distny <- gsub("\\.[[:alnum:]]+", "", family)
  y <- if (n > 0)
         switch(distny,
                gaussian = randgauss(n, nloc, ll, pary, dispersion),
                binomial = randqbinom(n, nloc, ll, pary, dispersion),
                poisson = randqpois(n, nloc, ll, pary, dispersion),
                Gamma = randgamma(n, nloc, ll, pary, dispersion))

  ## Return data.frame
  varsnm <- unique(c(all.vars(formula0), all.vars(atsample)))
  varsfm <- as.formula(paste("~", paste(varsnm, collapse = "+"), "+0"))
  varscall <- list(quote(stats::model.frame), formula = varsfm, data = data)
  outdf <- eval(as.call(varscall), parent.frame())
  varsnm <- names(outdf)
  if (!noweights && !(nmweights %in% varsnm)) {
    outdf <- cbind(outdf, weights)
    varsnm <- c(varsnm, nmweights)
  }
  NN <- nrow(outdf)
  nmall <- make.names(c(varsnm, rep(nmy, length = n),
                        if(returnGRF) rep(nmz, length = n)), unique = TRUE)
  nmyz <- nmall[(ncol(outdf) + 1):(length(nmall))]
  if (missing(subset) || is.null(subset)) subset <- 1:NN
  dfy <- matrix(NA, NN, n)
  dfy[subset, ] <- y
  if (returnGRF) {
    dfz <- matrix(NA, NN, n)
    dfz[subset, ] <- z
  } else dfz <- NULL
  dfyz <- cbind(dfy, dfz)
  colnames(dfyz) <- nmyz
  outdf <- cbind(outdf, dfyz)
  outdf
}

##' @name rsglmm
##' @export
rstrga <- function(n, formula, data, weights, subset, atsample,
                   beta, linkp, phi, omg, kappa, ssq,
                   corrfcn = "matern",
                   longlat = FALSE, dispersion = 1, returnGRF = FALSE) {

  call <- match.call()
  call[[1]] <- quote(rsglmm)
  call$family <- "gaussian"
  call$warndisp <- FALSE
  eval(call, parent.frame())
}


##' @name rsglmm
##' @export
rsgrf <- function (n, formula, data, subset, atsample,
                   beta, phi, omg, kappa, ssq,
                   corrfcn = "matern",
                   longlat = FALSE) {

  ## Response name
  formula <- as.formula(formula)
  if (length(formula) == 2) { # No LHS
    ## RHS <- formula[[2]]
    nmz <- "z"
  } else { # LHS names provided
    LHS <- formula[[2]]
    ## RHS <- formula[[3]]
    if ("|" %in% all.names(LHS)) { # LHS contains bar
      nmz <- all.vars(LHS[[3]])    # After the bar
    } else {
      nmz <- all.vars(LHS)
    }
  }
  formula0 <- update(formula, NULL ~ .)

  ## Design matrix and data
  if (missing(data)) data <- environment(formula)
  mfc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"),
             names(mfc), 0L)
  mfc <- mfc[c(1L, m)]
  mfc$formula <- formula0
  mfc$drop.unused.levels <- TRUE
  mfc$na.action <- "na.fail"
  mfc[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfc, parent.frame())
  mt <- attr(mf, "terms")
  FF <- model.matrix(mt,mf)

  ## All locations
  atsample <- update(atsample, NULL ~ . + 0) # No response and no intercept
  mfatc <- mfc
  mfatc$weights = NULL
  mfatc$formula = atsample
  mfat <- eval(mfatc, parent.frame())
  loc <- as.matrix(mfat)
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations.")

  z <- randgrf(n, loc, beta, phi, omg, kappa, ssq, FF, corrfcn, longlat)

  ## Return data.frame
  varsnm <- unique(c(all.vars(formula0), all.vars(atsample)))
  formulavars <- as.formula(paste("~", paste(varsnm, collapse = "+"), "+0"))
  varscall <- list(quote(stats::model.frame),
                   formula = formulavars, data = data)
  outdf <- eval(as.call(varscall),parent.frame())
  varsnm <- names(outdf)
  NN <- nrow(outdf)
  nmall <- make.names(c(varsnm, rep(nmz, length = n)), unique = TRUE)
  nmz <- nmall[(ncol(outdf) + 1):(length(nmall))]
  if (missing(subset) || is.null(subset)) subset <- 1:NN
  dfz <- matrix(NA, NN, n)
  dfz[subset, ] <- z
  colnames(dfz) <- nmz
  outdf <- cbind(outdf, dfz)
  outdf
}


##' @useDynLib geoBayes spcorr
randgrf <- function (n, loc, beta, phi, omg, kappa, ssq, FF,
                     corrfcn = "matern",
                     longlat = FALSE) {

  ## Check sample size
  n <- as.integer(n)
  if (n < 0) stop ("Negative sample size n entered.")

  ## Logical input
  longlat <- as.logical(longlat)
  if (length(longlat) != 1 || is.na(longlat))
    stop ("Input longlat must be a single logical value.")

  ## Correlation function
  icf <- .geoBayes_correlation(corrfcn)
  corrfcn <- .geoBayes_corrfcn$corrfcn[icf]
  kappa <- .geoBayes_getkappa(kappa, icf)
  phi <- as.double(phi)
  if (phi < 0) stop ("Input phi must >= 0")
  omg <- as.double(omg)
  if (omg < 0) stop ("Input omg must >= 0")

  ## Check beta, ssq
  beta <- as.vector(beta)
  ssq <- as.double(ssq)
  if (ssq <= 0) stop ("Non positive ssq entered.")

  ## Check loc
  loc <- as.matrix(loc)
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations.")
  if (corrfcn == "spherical" && ncol(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }

  ## Check model matrix
  FF <- as.matrix(FF)
  if (!all(is.finite(FF))) stop ("Non-finite values in the design matrix.")
  nloc <- nrow(FF)
  p <- ncol(FF)
  if (length(beta) != p)
    stop ("The length of beta does not match the number of predictors.")

  ## Covariance matrix
  dm <- sp::spDists(loc, longlat = longlat)
  dmup <- dm[upper.tri(dm)]
  ndm <- as.integer(length(dmup))
  spcorr <- .Fortran("spcorr", dmup, as.double(phi), as.double(kappa),
                     ndm, icf, PACKAGE = "geoBayes")[[1]]
  Sg <- diag(ssq*(1 + omg), nloc, nloc)
  Sg[upper.tri(Sg)] <- ssq*spcorr
  Sgch <- chol(Sg)

  ## Simulate
  meanz <- as.vector(FF %*% beta)
  z <- matrix(rnorm(n*nloc), nloc, n)
  z <- meanz + crossprod(Sgch, z)
  z
}

randgauss <- function (n, ns, weight, par, dispersion) {
  mean <- weight*par
  sd <- sqrt(weight*dispersion)
  y <- matrix(rnorm(n*ns, mean, sd), ns, n)
  y
}

randqbinom <- function (n, ns, weight, par, dispersion) {
  if (dispersion == 1) {
    y <- matrix(rbinom(n*ns, weight, par), ns, n)
  } else {
    ## y = dispersion*Bin(m/dispersion,p) +
    ## dispersion*eps*Beta(1+p,2-p)
    mt <- weight/dispersion
    mtm <- floor(mt)
    y <- matrix(rbinom(ns*n, mtm, par), ns, n)*dispersion
    ii <- mt > mtm
    ni <- sum(ii)
    if (ni > 0) {
      eps <- weight[ii] %% dispersion
      bet1 <- 1 + par[ii]
      bet2 <- 2 - par[ii]
      betr <- matrix(rbeta(ni*n, bet1, bet2), ni, n)*eps
      y[ii, ] <- y[ii, ] + betr
    }
  }
  y
}

randqpois <- function (n, ns, weight, par, dispersion) {
  if (dispersion == 1) {
    mean <- weight*par
    y <- matrix(rpois(n*ns, mean), ns, n)
  } else {
    mean <- weight*par/dispersion
    y <- matrix(rpois(n*ns, mean), ns, n)*dispersion
  }
  y
}

randgamma <- function (n, ns, weight, par, dispersion) {
  scale <- par*dispersion
  shape <- weight/dispersion
  y <- matrix(rgamma(n*ns, shape = shape, scale = scale), n, ns)
  y
}
