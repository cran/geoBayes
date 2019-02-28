##' Compute the standard errors for the Bayes factors estimates.
##'
##' Uses the batch means method to compute the standard errors for
##' Bayes factors.
##' @title Batch means, Bayes factors standard errors
##' @param pargrid A data frame with components "linkp", "phi", "omg",
##'   "kappa". Each row gives a combination of the parameters to
##'   compute the new standard errors.
##' @param runs A list with outputs from the function
##'   \code{\link{mcsglmm}} or \code{\link{mcstrga}}.
##' @param bfsize1 A scalar or vector of the same length as
##'   \code{runs} with all integer values or all values in (0, 1]. How
##'   many samples (or what proportion of the sample) to use for
##'   estimating the Bayes factors at the first stage. The remaining
##'   sample will be used for estimating the standard errors in the
##'   second stage. Setting it to 1 will perform only the first stage.
##' @param nbatch1 A scalar or vector of the same length as
##'   \code{runs}. All values must be integers or less than 1. This is
##'   used for calculating how many batches to split each of the
##'   sample in runs for the calculation of the Bayes factors standard
##'   errors for the parameters corresponding to \code{runs}.
##' @param nbatch2 A scalar or vector of the same length as
##'   \code{runs}. All values must be integers or less than 1. This is
##'   used for calculating how many batches to split each of the
##'   sample in runs for the calculation of the Bayes factors standard
##'   errors for the parameters corresponding to \code{pargrid}.
##' @param method Which method to use to calculate the Bayes factors:
##'   Reverse logistic or Meng-Wong.
##' @param bvmethod Which method to use for the calculation of the
##'   batch variance. The standard method splits to disjoint batches.
##'   The second and third method use the spectral variance method
##'   with different lag windows.
##' @param reference Which model goes in the denominator.
##' @param transf Whether to use a transformed sample for the
##'   computations. If \code{"no"} or \code{FALSE}, it doesn't. If
##'   \code{"mu"} or \code{TRUE}, it uses the samples for the mean. If
##'   \code{"wo"} it uses an alternative transformation. The latter
##'   can be used only for the families indicated by
##'   \code{.geoBayes_models$haswo}.
##' @return A list with components
##' \itemize{
##' \item \code{pargrid} The inputted pargrid augmented with the computed standard
##' errors.
##' \item \code{bfEstimate} The estimates of the Bayes factors
##' \item \code{bfSigma} The covariance matrix for the Bayes factors
##' estimates.
##' }
##' @references Roy, V., Tan, A. and Flegal, J. (2018). Estimating
##'   standard errors for importance sampling estimators with multiple
##'   Markov chains, Statistica Sinica, 28 1079-1101.
##'
##' Roy, V., & Evangelou, E. (2018). Selection of proposal
##'   distributions for generalized importance sampling estimators.
##'   arXiv preprint arXiv:1805.00829. 
##' @useDynLib geoBayes bfse_no bfse_mu bfse_wo bfse_tr
##' @export
bmbfse <- function(pargrid, runs, bfsize1 = 0.80, nbatch1 = 0.5, nbatch2 = 0.5,
                   method = c("RL", "MW"),
                   bvmethod = c("Standard", "TukeyHanning", "Bartlett"),
                   reference = 1, transf = c("no", "mu", "wo"))
{
  ## Method
  method <- match.arg(method)
  imeth <- match(method, eval(formals()$method))

  ## Batch variance method
  bvmethod <- match.arg(bvmethod)
  ibvmeth <- match(bvmethod, eval(formals()$bvmethod))

  ## Check input
  if (any(sapply(runs, class) != "geomcmc")) {
    stop ("Input runs is not a list with elements of class geomcmc")
  }
  nruns <- length(runs)
  if (nruns == 0) stop ("No runs specified")

  reference <- as.integer(reference)
  if (isTRUE(reference < 1L | reference > nruns)) {
    stop("Argument reference does not correspond to a run in runs")
  }

  ## Check if fixed phi and omg
  if (!all(sapply(runs, function(x) length(x$FIXED$phi) == 1))) {
    stop("Each input runs must have a fixed value phi.")
  }
  if (!all(sapply(runs, function(x) length(x$FIXED$omg) == 1))) {
    stop("Each input runs must have a fixed value omg.")
  }

  ## Extract data and model
  nm_DATA <- c("response", "weights", "modelmatrix", "locations",
               "longlat")
  nm_MODEL <- c("family", "corrfcn", "betm0", "betQ0", "ssqdf", "ssqsc",
                "tsqdf", "tsqsc", "dispersion")
  DATA <- runs[[1]]$DATA[nm_DATA]
  MODEL <- runs[[1]]$MODEL[nm_MODEL]
  if (nruns > 1) {
    for (i in 2:nruns) {
      if (!identical(runs[[i]]$DATA[nm_DATA], DATA)) {
        stop("MCMC chains don't all correspond to the same data.")
      }
      if (!identical(runs[[i]]$MODEL[nm_MODEL], MODEL)) {
        stop("MCMC chains don't all correspond to the same model.")
      }
    }
  }
  y <- DATA$response
  n <- as.integer(length(y))
  l <- DATA$weights
  F <- DATA$modelmatrix
  p <- NCOL(F)
  loc <- DATA$locations
  dm <- sp::spDists(loc, longlat = DATA$longlat)
  family <- MODEL$family
  ## ifam <- .geoBayes_family(family)
  corrfcn <- MODEL$corrfcn
  icf <- .geoBayes_correlation(corrfcn)
  betm0 <- MODEL$betm0
  betQ0 <- MODEL$betQ0
  ssqdf <- MODEL$ssqdf
  ssqsc <- MODEL$ssqsc
  tsqdf <- MODEL$tsqdf
  tsqsc <- MODEL$tsqsc
  dispersion <- MODEL$dispersion

  ## Skeleton points
  phi_pnts <- as.double(sapply(runs, function(r) r$FIXED$phi))
  omg_pnts <- as.double(sapply(runs, function(r) r$FIXED$omg))
  nu_pnts <- as.double(sapply(runs, function(r) r$FIXED$linkp_num))
  if (.geoBayes_corrfcn$needkappa[icf]) {
    kappa_pnts <- sapply(runs, function(r) r$FIXED$kappa)
    kappa_pnts <- .geoBayes_getkappa(kappa_pnts, icf)
  } else {
    kappa_pnts <- rep(0, nruns)
  }

  ## MCMC Sizes
  Nout <- sapply(runs, function(x) x$MCMC$Nout)
  Nout1 <- getsize(bfsize1, Nout, "*")
  Ntot1 <- sum(Nout1)
  Nout2 <- Nout - Nout1
  Ntot2 <- sum(Nout2)

  ## Sizes for SE
  nb1 <- getsize(nbatch1, Nout1, "^")
  nb2 <- if(Ntot2 > 0) getsize(nbatch2, Nout2, "^") else rep(0L, nruns)

  ## Transformed sample
  getsample <-
    transfsample(runs,
                 list(response = y, family = family), transf)
  sample <- matrix(unlist(getsample$sample), n)
  itr <- getsample$itr
  transf <- getsample$transf
  real_transf <- getsample$real_transf
  ifam <- getsample$ifam

  froutine <- paste0("bfse_", real_transf)
  if (real_transf == "wo" || real_transf == "tr") {
    ifam <- -ifam
  }

  ## New parameter values
  pargrid <- check_pargrid(pargrid, family, corrfcn)
  nnew <- nrow(pargrid)
  phi <- pargrid$phi
  omg <- pargrid$omg
  kappa <- pargrid$kappa
  nu <- pargrid$nu

  ## Split the sample
  sel1 <- rep(rep(c(TRUE, FALSE), nruns), rbind(Nout1, Nout2))
  sample <- matrix(unlist(sample), n)
  z1 <- sample[, sel1, drop = FALSE]
  z2 <- sample[, !sel1, drop = FALSE]

  ## dispersion or tsqsc
  if (ifam == 0) {
    tsq <- tsqsc
  } else {
    tsq <- dispersion
    tsqdf <- 0
  }

  ## Output
  SE <- numeric(nnew)
  bf <- numeric(nruns)
  Sig <- numeric((nruns - 1)^2)
  VT1 <- VT2 <- numeric(nnew)

  RUN <- .Fortran(froutine, bf = bf, Sig = Sig, SE = SE, VT1 = VT1,
                  VT2 = VT2, as.integer(reference),
                  as.double(phi), as.double(omg), as.double(nu), as.double(kappa),
                  as.double(phi_pnts), as.double(omg_pnts),
                  as.double(nu_pnts), as.double(kappa_pnts),
                  as.double(z1), as.integer(Nout1), as.integer(Ntot1),
                  as.double(z2), as.integer(Nout2), as.integer(Ntot2),
                  as.double(y), as.double(l), as.double(F), as.double(dm),
                  as.double(betm0), as.double(betQ0), as.double(ssqdf),
                  as.double(ssqsc), as.double(tsqdf), as.double(tsq),
                  as.integer(icf), as.integer(n), as.integer(p), as.integer(nnew),
                  as.integer(nruns), as.integer(ifam),
                  as.integer(imeth), as.integer(nb1), as.integer(nb2),
                  as.integer(ibvmeth), as.integer(itr), PACKAGE = "geoBayes")

  ## Return
  out <- list()
  pargrid$SE <- RUN$SE/sqrt(Ntot2)
  out$pargrid <- pargrid
  out$bfEstimate <- RUN$bf
  out$bfSigma <- matrix(RUN$Sig/Ntot1, nruns - 1)
  out$VT1 <- RUN$VT1/Ntot2
  out$VT2 <- RUN$VT2/Ntot2
  return (out)
}


groupColMeans <- function (x, i) {
  ## Compute the mean at each column of x by subsetting its rows
  ## according to the factor i. Returns a matrix with as many columns
  ## as x and as many rows as the number of levels in i.
  x <- as.matrix(x)
  rx <- nrow(x)
  cx <- ncol(x)
  m <- tapply(seq_len(rx), i, function(jj) colMeans(x[jj, , drop = FALSE]),
              simplify = FALSE)
  matrix(unlist(m), ncol = cx, byrow = TRUE)
}

bmmcse <- function(x, size)
{
  ## Purpose: Compute batch-means standard error of univariate MC samples
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x   A matrix of MC samples. The SE corresponding to each
  ## column of x will be returned.
  ## size   The size of each batch.
  ## ----------------------------------------------------------------------
  x <- as.matrix(x)
  n <- nrow(x)
  if (missing(size)) {
    size <- "sqroot"
  } else if (length(size) != 1) stop ("Argument size must have length 1.")
  if (is.character(size)) {
    i <- pmatch(size, c("sqroot", "sqrt", "cuberoot", "cubert"), nomatch = 0L,
                duplicates.ok = FALSE)
    if (i == 1 | i == 2) {
      size <- as.integer(sqrt(n))
    } else if (i == 3 | i == 4) {
      size <- as.integer(n^(1/3))
    } else stop ("Character size must be one of sqroot or cuberoot.")
  } else {
    size <- as.numeric(size)
    if (size <= 0) {
      stop ("Argument size must be positive.")
    } else if (size <= 1) {
      size <- as.integer(size*n)
    } else if (size > n) stop ("Argument size is larger than sample size.")
  }
  ## Determine the number of batches and the size of each batch
  nbatches <- as.integer(n/size)
  if (nbatches <= 1L) stop ("Number of batches too small. Reduce size.")
  nrem <- n%%nbatches
  batchl <- rep(c(size + 1L, size), c(nrem, nbatches - nrem))
  batchf <- rep(seq_len(nbatches), batchl)
  ## Compute batch mean and overall mean
  bm <- groupColMeans(x, batchf)
  m <- colSums(bm*batchl)/n
  ## Compute the batch variance
  bmmm <- sweep(bm, 2, m, check.margin = FALSE)
  bv <- colSums(batchl*bmmm^2)/(nbatches - 1)
  sqrt(bv)
}


## .. content for \description{} (no empty lines) ..
##
## .. content for \details{} ..
## @title
## @param sample A list containing the samples from each run. Each
##   element of the list must be a matrix or a vector. If a matrix,
##   each column represents one sample.
## @param theta0 A matrix or vector containing the parameter values
##   corresponding to each run in \code{sample}. Each column of
##   \code{theta0} is a parameter.
## @param theta1 A matrix or vector containing the new parameter
##   values for which to compute the standard errors. Each column of
##   \code{theta1} is a parameter.
## @param llikfun A function that computes the log-likelihood. The
##   function is called as \code{llikfun(x,theta,...)} where x is a
##   sample and theta is a parameter vector. The function must return
##   a scalar for each \code{x,theta} pair.
## @param MoreArgs Further arguments to \code{llikfun}.
## @param bfsize1 A scalar or vector of the same length as
##   \code{sample} with all integer values or all values in (0, 1]. How
##   many samples (or what proportion of the sample) to use for
##   estimating the Bayes factors at the first stage. The remaining
##   sample will be used for estimating the standard errors in the
##   second stage. Setting it to 1 will perform only the first stage.
## @param nbatch1 A scalar or vector of the same length as \code{sample}. All
##   values must be integers or less than 1. This is used for
##   calculating how many batches to split each of the sample in runs
##   for the calculation of the Bayes factors standard errors for the
##   parameters corresponding to \code{theta0}.
## @param nbatch2 A scalar or vector of the same length as \code{sample}. All
##   values must be integers or less than 1. This is used for
##   calculating how many batches to split each of the sample in runs
##   for the calculation of the Bayes factors standard errors for the
##   parameters corresponding to \code{theta1}.
## @param method Which method to use to calculate the Bayes factors:
## Reverse logistic or Meng-Wong.
## @param bvmethod Which method to use for the calculation of the
##   batch variance. The standard method splits to disjoint batches.
##   The second and third method use the spectral variance method
##   with different lag windows.
## @param reference Which model goes in the denominator.
## @return A list with components
## \itemize{
## \item \code{logbf} The estimated Bayes factors at \code{theta0}.
## \item \code{bfSigma} The covariance matrix for the Bayes factors
## estimates in \code{logbf}.
## \item \code{SE} The standard error of the Bayes factors estimates
##   at \code{theta1}
## \item \code{VT1, VT2} The first and second variance terms for
##   \code{SE} above.
## }
bfsecalc <- function(sample, theta0, theta1, llikfun, MoreArgs = NULL,
                     bfsize1 = 0.8, nbatch1 = 0.5, nbatch2 = 0.5,
                     method = c("RL", "MW"),
                     bvmethod = c("Standard", "TukeyHanning", "Bartlett"),
                     reference = 1)
{
  ## Method
  method <- match.arg(method)
  imeth <- match(method, eval(formals()$method))

  ## Batch variance method
  bvmethod <- match.arg(bvmethod)
  ibvmeth <- match(bvmethod, eval(formals()$bvmethod))

  kg <- length(sample)

  reference <- as.integer(reference)
  if (reference < 1 || reference > kg)
    stop ("Argument reference must be in between 1 and number of parameters.")

  nCOL <- function(x) if (length(d <- dim(x)) > 1L) d[2L] else length(x)
  Nout <- sapply(sample, nCOL)
  Ntot <- sum(Nout)
  Nout1 <- getsize(bfsize1, Nout, "*")
  Ntot1 <- sum(Nout1)
  Nout2 <- Nout - Nout1
  Ntot2 <- sum(Nout2)

  ## Size of x
  nROW <- function(x) if (length(d <- dim(x)) > 1L) d[1L] else 1L
  n <- unique(sapply(sample, nROW))
  if (length(n) != 1L) stop ("The dimension of the sampled variable changes.")

  ## Check input theta
  ltht <- nROW(theta0)
  if (nCOL(theta0) != kg)
    stop ("Columns in theta0 must be the length of samples.")
  if (nROW(theta1) != ltht)
    stop ("theta0 and theta1 must have the same number of rows.")
  nnew <- nCOL(theta1)
  dim(theta0) <- c(ltht, kg)
  dim(theta1) <- c(ltht, nnew)

  ## Sizes for SE
  nb1 <- getsize(nbatch1, Nout1, "^")
  nb2 <- getsize(nbatch2, Nout2, "^")

  ## Split the sample
  sel1 <- rep(rep(c(TRUE, FALSE), kg), rbind(Nout1, Nout2))
  sample <- matrix(unlist(sample), n)
  z1 <- sample[, sel1, drop = FALSE]
  z2 <- sample[, !sel1, drop = FALSE]

  ## Compute log-likelihoods
  llik1 <- mapply(llikfun, do.call("c", apply(z1, 2, list)),
                  rep(do.call("c", apply(theta0, 2, list)), each = Ntot1),
                  MoreArgs = MoreArgs)
  llik2 <- mapply(llikfun, do.call("c", apply(z2, 2, list)),
                  rep(do.call("c", apply(theta0, 2, list)), each = Ntot2),
                  MoreArgs = MoreArgs)
  llikn <- mapply(llikfun, do.call("c", apply(z2, 2, list)),
                  rep(do.call("c", apply(theta1, 2, list)), each = Ntot2),
                  MoreArgs = MoreArgs)

  ## Compute SE
  bf <- numeric(kg)
  SE <- VT1 <- VT2 <- numeric(nnew)
  Sig <- matrix(0, kg-1, kg-1)
  Bet <- Omg <- matrix(0, kg, kg)
  fcall <- .Fortran("bfsecalc", bf, Sig, SE, VT1, VT2, reference,
                    llik1, llik2, llikn, Nout1, Ntot1, Nout2, Ntot2,
                    nnew, kg, imeth, nb1, nb2, ibvmeth, Bet, Omg,
                    PACKAGE = "geoBayes")

  Bet <- fcall[[20]]; Bet[lower.tri(Bet)] <- t(Bet)[lower.tri(Bet)]
  Omg <- fcall[[21]]; Omg[lower.tri(Omg)] <- t(Omg)[lower.tri(Omg)]
  BMP <- chol2inv(chol(Bet + 1/kg)) - 1/kg
  BOB <- t(BMP) %*% Omg %*% BMP

  ## Return
  out <- fcall[1:5]
  names(out) <- c("logbf", "bfSigma", "SE", "VT1", "VT2")
  out$Varlogbf <- BOB
  out
}
