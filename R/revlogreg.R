##' Perform the reverse logistic regression estimation
##'
##' Estimation is done by maximising the reverse logistic log likelihood.
##' @title Reverse logistic regression estimation
##' @param lglk The value of the loglikelihood at different samples and
##'   different parameters. This should be entered as a matrix where the
##'   rows are the values of the samples and the columns correspond to
##'   the parameters. The
##'   [i,j] element of the matrix is the value of the loglikelihood at
##'   the ith sample when all samples are put together evaluated
##'   at the jth parameter value.
##' @param N A vector of length ncol(lglk) or a scalar corresponding to
##'   the sample sizes from each model. Must sum(N) == nrow(lglk). The
##'   first N[1] samples come from model corresponding to the first set
##'   of parameters, then (N[1]+1):N[2] are from the model corresponding
##'   to the second set of parameters, and so on.
##' @return A vector containing the reverse logistic regression estimates
##'   of the logarithm of the Bayes factors.
##'   The first set of parameters is taken as the reference model so its
##'   estimate is always 0.
##' @references Geyer, C. J. (1994). Estimating normalizing constants
##'   and reweighting mixtures in Markov chain Monte Carlo. Technical
##'   Report 568, School of Statistics, University of Minnesota.
##' @importFrom stats optim
##' @export
revlogreg <- function(lglk, N)
{
  ## Purpose: Reverse logistic regression
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## lglk   The log-likelihood. Each column corresponds to
  ## the value of the log-likelihood evaluated at different parameters.
  ## N   A vector of length ncol(lglk) or a scalar corresponding to
  ## the sample sizes from each model. sum(N) == nrow(lglk)

  lglk <- as.matrix(lglk)
  nruns <- ncol(lglk)
  Ntot <- nrow(lglk)
  N <- rep(as.integer(N), length = nruns)
  if (sum(N) != Ntot) stop ("Inconsistent sample sizes implied by lglk and N")
  ## if (any(!is.finite(lglk))) stop ("Log-likelihood contains non-finite values")

  if (nruns == 1) return (0)

  isum <- matrix(rep(c(TRUE, rep(FALSE, nruns)), length.out = nruns*nruns),
                 nruns)[rep(1:nruns, N), ]

  fun <- function (eta) {
    lliketa <- lglk + matrix(eta, Ntot, nruns, TRUE)
    mxlliketa <- apply(lliketa, 1, max)
    lliketammx <- lliketa - mxlliketa
    lgdenom <- mxlliketa + log(rowSums(exp(lliketammx)))
    lgnum <- lliketa[isum]
    out <- sum(lgnum - lgdenom)
    -out
  }

  gun <- function (eta) {
    lliketa <- lglk + matrix(eta, Ntot, nruns, TRUE)
    mxlliketa <- apply(lliketa, 1, max)
    lliketammx <- lliketa - mxlliketa
    lgp <- lliketammx - log(rowSums(exp(lliketammx)))
    out <- N - colSums(exp(lgp))
    -out
  }

  op <- stats::optim(log(N/Ntot), fun, gun, method = "BFGS")
  out <- op$par - log(N/Ntot)
  out[1] - out
}
