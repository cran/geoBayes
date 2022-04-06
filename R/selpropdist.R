###############################################################################
#### Code for the paper "Selection of Proposal Distributions for
#### Multiple Importance Sampling" by V Roy and E Evangelou
###############################################################################

##' Selection of multiple importance sampling distributions
##'
##' \describe{
##' \item{SEQ}{is a sequential method starting with \code{istart} and
##' additng to it until \code{K} proposals have been selected. At each
##' iteration, the point with the highest (relative?) standard error
##' is added} 
##' \item{MNX}{is the minimax method. The chosen proposal corresponds
##' to the lowest maximum (relative?) standard error.}
##' \item{ENT}{is the entropy method. The chosen proposal corresponds
##' to the highest determinant of the (relative?) covariance matrix at
##' the first stage.}
##' }
##' @title Selection of multiple importance sampling distributions
##' @param pargrid A data frame with components "linkp", "phi", "omg",
##'   "kappa". Each row gives a combination of the parameters to
##'   compute the new standard errors.
##' @param K How many proposal densities in total to choose among the
##'   rows of \code{pargrid}. Needed for SEQ only. For MNX and ENT
##'   this is determined by the length of \code{istart}. 
##' @param istart Start with these rows of \code{pargrid}. A vector of
##'   indices.
##' @param relativeSE Logical. Whether the choice is based on the
##'   standard error (FALSE), or relative standard error (TRUE).
##' @param N1 The sample size for stage 1.
##' @param N2 The sample sie for stage 2.
##' @param Nthin Thinning
##' @param Nbi Burn-in
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
##' @param S1method Which method to use to calculate the Bayes
##'   factors: Reverse logistic or Meng-Wong.
##' @param bvmethod Which method to use for the calculation of the
##'   batch variance. The standard method splits to disjoint batches.
##'   The second and third method use the spectral variance method
##'   with different lag windows.
##' @param transf Whether to use a transformed sample for the
##'   computations. If \code{"no"} or \code{FALSE}, it doesn't. If
##'   \code{"mu"} or \code{TRUE}, it uses the samples for the mean. If
##'   \code{"wo"} it uses an alternative transformation. The latter
##'   can be used only for the families indicated by
##'   \code{.geoBayes_models$haswo}.
##' @return A list with components
##' \describe{
##' \item{selected}{The rows of \code{pargrid} selected.}
##' \item{isel}{The indices of the rows of \code{pargrid} selected.}
##' \item{se}{The standard error corresponding to the selected
##'   parameters.}
##' \item{samples}{A list containing the samples from the selected
##'   parameters.} 
##' }
##' @name select_proposals
##' @references Roy, V., & Evangelou, E. (2018). Selection of proposal
##'   distributions for generalized importance sampling estimators.
##'   arXiv preprint arXiv:1805.00829. 
##' @export
##' @examples
##' \dontrun{
##' data(rhizoctonia)
##' ### Define the model
##' corrf <- "spherical"
##' kappa <- 0
##' ssqdf <- 1
##' ssqsc <- 1
##' betm0 <- 0
##' betQ0 <- .01
##' family <- "binomial.probit"
##' formula <- Infected ~ 1
##' atsample <- ~ Xcoord + Ycoord
##' ### Skeleton points
##' philist <- seq(100, 200, 10)
##' omglist <- 0
##' parlist <- expand.grid(linkp=0, phi=philist, omg=omglist, kappa = kappa)
##' ### MCMC sizes
##' Nout <- 100
##' Nthin <- 1
##' Nbi <- 10
##' ## Select proposals
##' K <- 3                        # Choose 3 proposals
##' istart_SEQ <- 6               # Start with middle
##' istart_MNX <- istart_ENT <- c(6, 2, 10)
##' cooling_MNX <- .05/log((0:24 %/% 5)*5 + exp(1))
##' cooling_ENT <- .3/log((0:49 %/% 10)*10 + exp(1))
##' prop_SEQ <- select_proposals_SEQ(pargrid = parlist, K = K,
##'                                  istart = istart_SEQ,
##'                                  relativeSE = TRUE, 
##'                                  N1 = Nout, N2 = Nout,
##'                                  Nthin = Nthin, Nbi = Nbi,
##'                                  formula = formula, family = family,
##'                                  data = rhizoctonia, weights = Total, 
##'                                  atsample = atsample, corrfcn = corrf,
##'                                  betm0 = betm0, betQ0 = betQ0,
##'                                  ssqdf = ssqdf, ssqsc = ssqsc,
##'                                  dispersion = 1, longlat = FALSE,
##'                                  nbatch1 = 0.5, nbatch2 = 0.5,
##'                                  bvmethod = "TukeyHanning",
##'                                  transf = "mu")
##' prop_MNX <- select_proposals_MNX(pargrid = parlist,
##'                                  istart = istart_MNX, nfix = 1L,
##'                                  cooling = cooling_MNX, 
##'                                  relativeSE = TRUE, 
##'                                  N1 = Nout, N2 = Nout,
##'                                  Nthin = Nthin, Nbi = Nbi,
##'                                  formula = formula, family = family,
##'                                  data = rhizoctonia, weights = Total, 
##'                                  atsample = atsample, corrfcn = corrf,
##'                                  betm0 = betm0, betQ0 = betQ0,
##'                                  ssqdf = ssqdf, ssqsc = ssqsc,
##'                                  dispersion = 1, longlat = FALSE,
##'                                  nbatch1 = 0.5, nbatch2 = 0.5,
##'                                  bvmethod = "TukeyHanning",
##'                                  transf = "mu",
##'                                  verbose = TRUE)
##' prop_ENT <- select_proposals_ENT(pargrid = parlist,
##'                                  istart = istart_ENT, nfix = 1L,
##'                                  cooling = cooling_ENT, 
##'                                  relativeSE = TRUE, 
##'                                  N1 = Nout, 
##'                                  Nthin = Nthin, Nbi = Nbi,
##'                                  formula = formula, family = family,
##'                                  data = rhizoctonia, weights = Total, 
##'                                  atsample = atsample, corrfcn = corrf,
##'                                  betm0 = betm0, betQ0 = betQ0,
##'                                  ssqdf = ssqdf, ssqsc = ssqsc,
##'                                  dispersion = 1, longlat = FALSE,
##'                                  nbatch1 = 0.5, nbatch2 = 0.5,
##'                                  bvmethod = "TukeyHanning",
##'                                  transf = "mu",
##'                                  verbose = TRUE)
##' }
select_proposals_SEQ <- function (pargrid, K, istart,
                                  relativeSE = FALSE, 
                                  N1, N2, Nthin, Nbi,
                                  formula, family = "gaussian",
                                  data, weights, subset, offset, 
                                  atsample, corrfcn = "matern",
                                  betm0, betQ0, ssqdf, ssqsc,
                                  dispersion = 1, longlat = FALSE,
                                  nbatch1 = 0.5, nbatch2 = 0.5,
                                  S1method = c("RL", "MW"),
                                  bvmethod = c("Standard", "TukeyHanning",
                                               "Bartlett"),
                                  transf = c("no", "mu", "wo"))
{
  cl <- match.call()
  pargrid <- check_pargrid(pargrid, family, corrfcn)
  partry <- pargrid[, c("nu", "phi", "omg", "kappa")]
  ngrid <- nrow(pargrid)
  istart <- unique(as.integer(istart))
  if (any(istart < 1L) || any(istart) > ngrid) {
    stop ("istart elements must be between 1 and nrow(pargrid).")
  }
  K <- as.integer(K)
  if (K <= 0L) stop("K must be a positive integer.")
  if (K > ngrid) stop("K cannot be more than the number of elements in pargrid.")
  N1 <- as.integer(N1)
  N2 <- as.integer(N2)
  Nout <- c(N1, N2)
  samples <- list()
  length(samples) <- K
  args_mcsglmm <- names(formals(mcsglmm))
  args_mcsglmm <- args_mcsglmm[args_mcsglmm %in% names(cl)[-1]]
  clsample <- cl[c("", args_mcsglmm)]
  clsample[[1]] <- quote(mcsglmm)
  clsample$corrtuning <- list(phi = 0, omg = 0, kappa = 0)
  clsample$Nout <- Nout
  clsample$test <- FALSE
  samplefun <- function(call, pars) {
    call$linkp <- pars[[1]]
    call$phi <- pars[[2]]
    call$omg <- pars[[3]]
    call$kappa <- pars[[4]]
    eval(call)
  }
  nstart <- length(istart)
  if (nstart > K) stop("Need at most K indices in istart.")
  iin <- integer(K)
  iin[1:nstart] <- istart
  for (i in seq_len(nstart)) {
    samples[[i]] <- samplefun(clsample, partry[iin[i], ])
  }
  se <- bmbfse(pargrid = pargrid,
               runs = samples[1L:nstart],
               bfsize1 = Nout[1L],
               nbatch1 = nbatch1,
               nbatch2 = nbatch2,
               S1method = S1method,
               bvmethod = bvmethod,
               transf = transf)
  relativeSE <- as.logical(relativeSE)
  if (isTRUE(relativeSE)) {
    objfun <- function(s) s$pargrid$SE*exp(-s$pargrid$logbf)
  } else {
    objfun <- function(s) s$pargrid$SE
  }
  if (K > nstart) {
    for (i in (nstart+1L):K) {
      sevals <- objfun(se)
      sevals[iin[1L:(i-1L)]] <- NA
      iin[i] <- which.max(sevals)
      samples[[i]] <- samplefun(clsample, partry[iin[i], ])
      se <- bmbfse(pargrid = pargrid,
                   runs = samples[1L:i],
                   bfsize1 = Nout[1L],
                   nbatch1 = nbatch1,
                   nbatch2 = nbatch2,
                   S1method = S1method,
                   bvmethod = bvmethod,
                   transf = transf)
    }
  }
  list(selected = pargrid[iin, ], isel = iin, se = se,
       samples = samples)
}

##' @param cooling A decreasing sequence of temperature values for the
##'   simulated annealing. All elements must be positive. A suggested
##'   value is \code{Tinit / log(((0:N) \%/\% Tstp)*Tstp + exp(1))} for
##'   \code{N+1} iterations, where \code{Tinit} is the initial
##'   temperature and \code{Tstp} is the number of iterations before
##'   the temperature is reduced.
##' @param nfix In the case of MNX and ENT, the first \code{nfix}
##'   elements of \code{istart} will always be included.
##' @param verbose Logical. Prints information about the simulated
##'   annealing. 
##' @name select_proposals
##' @export
##' @importFrom stats runif
select_proposals_MNX <- function (pargrid, istart, nfix, 
                                  relativeSE = FALSE, 
                                  N1, N2, Nthin, Nbi,
                                  cooling, 
                                  formula, family = "gaussian",
                                  data, weights, subset, offset, 
                                  atsample, corrfcn = "matern",
                                  betm0, betQ0, ssqdf, ssqsc,
                                  dispersion = 1, longlat = FALSE,
                                  nbatch1 = 0.5, nbatch2 = 0.5,
                                  S1method = c("RL", "MW"),
                                  bvmethod = c("Standard", "TukeyHanning",
                                               "Bartlett"),
                                  transf = c("no", "mu", "wo"),
                                  verbose = FALSE)
{
  cl <- match.call()
  verbose <- isTRUE(as.logical(verbose))
  pargrid <- check_pargrid(pargrid, family, corrfcn)
  partry <- pargrid[, c("nu", "phi", "omg", "kappa")]
  ngrid <- nrow(pargrid)
  istart <- unique(as.integer(istart))
  if (any(istart < 1L) || any(istart) > ngrid) {
    stop ("istart elements must be between 1 and nrow(pargrid).")
  }
  K <- length(istart)
  if (K <= 0L) stop("K must be a positive integer.")
  if (K > ngrid) stop("K cannot be more than the number of elements in pargrid.")
  nfix <- as.integer(nfix)
  N1 <- as.integer(N1)
  N2 <- as.integer(N2)
  Nout <- c(N1, N2)
  samples <- list()
  length(samples) <- ngrid
  cooling <- as.numeric(cooling)
  niter <- length(cooling)
  if (niter && !all(cooling > 0)) {
    stop ("Input cooling must contain only positive values.")
  }
  if (niter > 1L && any(cooling[2L:niter] > cooling[1L:(niter-1L)])) {
    stop ("Input cooling must be a decreasing sequence.")
  }
  args_mcsglmm <- names(formals(mcsglmm))
  args_mcsglmm <- args_mcsglmm[args_mcsglmm %in% names(cl)[-1]]
  clsample <- cl[c("", args_mcsglmm)]
  clsample[[1]] <- quote(mcsglmm)
  clsample$corrtuning <- list(phi = 0, omg = 0, kappa = 0)
  clsample$Nout <- Nout
  clsample$test <- FALSE
  samplefun <- function(call, pars) {
    call$linkp <- pars[[1]]
    call$phi <- pars[[2]]
    call$omg <- pars[[3]]
    call$kappa <- pars[[4]]
    eval(call)
  }
  for (i in seq_len(K)) {
    samples[[istart[i]]] <- samplefun(clsample, partry[istart[i], ])
  }
  se <- bmbfse(pargrid = pargrid,
               runs = samples[istart],
               bfsize1 = Nout[1L],
               nbatch1 = nbatch1,
               nbatch2 = nbatch2,
               S1method = S1method,
               bvmethod = bvmethod,
               transf = transf)
  relativeSE <- as.logical(relativeSE)
  if (isTRUE(relativeSE)) {
    objfun <- function(s) -max(s$pargrid$SE*exp(-s$pargrid$logbf))
  } else {
    objfun <- function(s) -max(s$pargrid$SE)
  }
  id_bst <- id_cur <- istart
  se_bst <- se
  s_bst <- s_cur <- objfun(se)
  lsam <- K - nfix
  ltry <- ngrid - K
  set <- seq_len(ngrid)
  ## Simulated annealing
  if (verbose) {
    cat("     Iter   Objective", "\n")
    cat("---------   ---------", "\n")
    cat(format(0L, width = 9), "   ", format(s_cur, width = 9), "\n")
  }
  for (i in seq_len(niter)) {
    set_new <- set[-id_cur]
    j <- nfix+sample.int(lsam, 1L)
    k <- set_new[sample.int(ltry, 1L)]
    id_new <- id_cur
    id_new[j] <- k
    if (is.null(samples[[k]])) {
      samples[[k]] <- samplefun(clsample, partry[k, ])
    }
    se <- bmbfse(pargrid = pargrid,
                 runs = samples[id_new],
                 bfsize1 = Nout[1L],
                 nbatch1 = nbatch1,
                 nbatch2 = nbatch2,
                 S1method = S1method,
                 bvmethod = bvmethod,
                 transf = transf)
    s_new <- objfun(se)
    accept <- runif(1L) < exp((s_new - s_cur)/cooling[i])
    if (accept) {
      id_cur <- id_new
      s_cur <- s_new
      if (s_new > s_bst) {
        id_bst <- id_new
        s_bst <- s_new
        se_bst <- se
      }
    }
    if(verbose) {
      cat(format(i, width = 9), "   ", format(s_cur, width = 9), "\n")
    }
  }
  list(selected = pargrid[id_bst, ], isel = id_bst, bmbfse = se_bst,
       samples = samples[id_bst])
}


##' @name select_proposals
##' @export 
select_proposals_ENT <- function (pargrid, istart, nfix, 
                                  relativeSE = FALSE, 
                                  N1, Nthin, Nbi,
                                  cooling, 
                                  formula, family = "gaussian",
                                  data, weights, subset, offset, 
                                  atsample, corrfcn = "matern",
                                  betm0, betQ0, ssqdf, ssqsc,
                                  dispersion = 1, longlat = FALSE,
                                  nbatch1 = 0.5, nbatch2 = 0.5,
                                  S1method = c("RL", "MW"),
                                  bvmethod = c("Standard", "TukeyHanning",
                                               "Bartlett"),
                                  transf = c("no", "mu", "wo"),
                                  verbose = FALSE)
{
  cl <- match.call()
  verbose <- isTRUE(as.logical(verbose))
  pargrid <- check_pargrid(pargrid, family, corrfcn)
  partry <- pargrid[, c("nu", "phi", "omg", "kappa")]
  ngrid <- nrow(pargrid)
  istart <- unique(as.integer(istart))
  if (any(istart < 1L) || any(istart) > ngrid) {
    stop ("istart elements must be between 1 and nrow(pargrid).")
  }
  K <- length(istart)
  if (K <= 0L) stop("K must be a positive integer.")
  if (K > ngrid) stop("K cannot be more than the number of elements in pargrid.")
  nfix <- as.integer(nfix)
  N1 <- as.integer(N1)
  N2 <- 0L
  Nout <- c(N1, N2)
  samples <- list()
  length(samples) <- ngrid
  cooling <- as.numeric(cooling)
  niter <- length(cooling)
  if (niter && !all(cooling > 0)) {
    stop ("Input cooling must contain only positive values.")
  }
  if (niter > 1L && any(cooling[2L:niter] > cooling[1L:(niter-1L)])) {
    stop ("Input cooling must be a decreasing sequence.")
  }
  args_mcsglmm <- names(formals(mcsglmm))
  args_mcsglmm <- args_mcsglmm[args_mcsglmm %in% names(cl)[-1]]
  clsample <- cl[c("", args_mcsglmm)]
  clsample[[1]] <- quote(mcsglmm)
  clsample$corrtuning <- list(phi = 0, omg = 0, kappa = 0)
  clsample$Nout <- Nout
  clsample$test <- FALSE
  samplefun <- function(call, pars) {
    call$linkp <- pars[[1]]
    call$phi <- pars[[2]]
    call$omg <- pars[[3]]
    call$kappa <- pars[[4]]
    eval(call)
  }
  for (i in seq_len(K)) {
    samples[[istart[i]]] <- samplefun(clsample, partry[istart[i], ])
  }
  se <- bmbfse(pargrid = pargrid,
               runs = samples[istart],
               bfsize1 = Nout[1L],
               nbatch1 = nbatch1,
               nbatch2 = nbatch2,
               S1method = S1method,
               bvmethod = bvmethod,
               transf = transf)
  relativeSE <- as.logical(relativeSE)
  if (isTRUE(relativeSE)) {
    objfun <- function(s)
      sum(log(diag(chol(s$bfSigma/tcrossprod(s$bfEstimate[-1])))))
  } else {
    objfun <- function(s) sum(log(diag(chol(s$bfSigma))))
  }
  id_bst <- id_cur <- istart
  se_bst <- se
  s_bst <- s_cur <- objfun(se)
  lsam <- K - nfix
  ltry <- ngrid - K
  set <- seq_len(ngrid)
  ## Simulated annealing
  if (verbose) {
    cat("     Iter   Objective", "\n")
    cat("---------   ---------", "\n")
    cat(format(0L, width = 9), "   ", format(s_cur, width = 9), "\n")
  }
  for (i in seq_len(niter)) {
    set_new <- set[-id_cur]
    j <- nfix+sample.int(lsam, 1L)
    k <- set_new[sample.int(ltry, 1L)]
    id_new <- id_cur
    id_new[j] <- k
    if (is.null(samples[[k]])) {
      samples[[k]] <- samplefun(clsample, partry[k, ])
    }
    se <- bmbfse(pargrid = pargrid,
                 runs = samples[id_new],
                 bfsize1 = Nout[1L],
                 nbatch1 = nbatch1,
                 nbatch2 = nbatch2,
                 S1method = S1method,
                 bvmethod = bvmethod,
                 transf = transf)
    s_new <- objfun(se)
    accept <- runif(1L) < exp((s_new - s_cur)/cooling[i])
    if (accept) {
      id_cur <- id_new
      s_cur <- s_new
      if (s_new > s_bst) {
        id_bst <- id_new
        s_bst <- s_new
        se_bst <- se
      }
    }
    if(verbose) {
      cat(format(i, width = 9), "   ", format(s_cur, width = 9), "\n")
    }
  }
  list(selected = pargrid[id_bst, ], isel = id_bst, bmbfse = se_bst,
       samples = samples[id_bst])
}
