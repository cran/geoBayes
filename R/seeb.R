##' Standard errors for the empirical Bayes estimates of the parameters.
##'
##' @title Empirical Bayes standard errors
##' @param mcrun The output from the function \code{mcsglmm} where the
##'   parameters linkp, phi, omg, kappa are set at their empirical
##'   Bayes estimates (output of \code{\link{bf2optim}}).
##' @param transf The type of transformation to use.
##' @return The precision (inverse covariance) matrix.
##' @importFrom stats var
##' @useDynLib geoBayes llikfcn_dh_tr
##' @references Casella, G. (2001). Empirical Bayes Gibbs sampling.
##'   Biostatistics, 2(4), 485-500.
##'
##' Evangelou, E., & Roy, V. (2019). Estimation and prediction for
##'   spatial generalized linear mixed models with parametric links
##'   via reparameterized importance sampling. Spatial Statistics, 29,
##'   289-315.
##' @export
bf2se <- function(mcrun, transf = c("no", "mu", "wo"))
{
  family <- mcrun$MODEL$family
  ifam <- .geoBayes_family(family)
  corrfcn <- mcrun$MODEL$corrfcn
  icf <- .geoBayes_correlation(corrfcn)
  F <- mcrun$DATA$modelmatrix
  n <- NROW(F)
  p <- NCOL(F)
  loc <- mcrun$DATA$locations
  longlat <- mcrun$DATA$longlat
  dm <- spDists(loc, longlat = longlat)
  y <- mcrun$DATA$response
  l <- mcrun$DATA$weights
  linkp <- mcrun$FIXED$linkp_num[1]
  if (length(linkp) < 1) stop ("This function requires fixed linkp.")
  nu <- .geoBayes_getlinkp(linkp, ifam)
  phi <- mcrun$FIXED$phi[1]
  if (length(phi) < 1) stop ("This function requires fixed phi.")
  omg <- mcrun$FIXED$omg[1]
  if (length(omg) < 1) stop ("This function requires fixed omg.")
  kappa <- mcrun$FIXED$kappa[1]
  if (length(kappa) < 1) stop ("This function requires fixed kappa.")
  betm0 <- mcrun$MODEL$betm0
  betQ0 <- mcrun$MODEL$betQ0
  ssqdf <- mcrun$MODEL$ssqdf
  ssqsc <- mcrun$MODEL$ssqsc
  tsqdf <- mcrun$MODEL$tsqdf
  tsqsc <- mcrun$MODEL$tsqsc
  dispersion <- mcrun$MODEL$dispersion
  if (ifam == 0) {
    tsq <- tsqsc
  } else {
    tsq <- dispersion
    tsqdf <- 0
  }
  getsample <- transfsample(list(mcrun),
                            list(response = y, family = family), transf)
  sample <- matrix(unlist(getsample$sample), n)
  itr <- getsample$itr
  transf <- getsample$transf
  real_transf <- getsample$real_transf
  ifam <- getsample$ifam
  froutine <- paste0("llikfcn_dh_", real_transf)
  froutine <- "llikfcn_dh_tr" # XXX
  nsim <- NCOL(sample)
  dlglk <- array(0, c(4, nsim))
  hlglk <- array(0, c(16, nsim))
  fff <- .Fortran(froutine, dlglk, hlglk, phi, omg, linkp, kappa,
                  sample, nsim, y, l, F, dm, betm0, betQ0, ssqdf, ssqsc,
                  tsqdf, tsq, icf, n, p, ifam, itr, PACKAGE = "geoBayes")
  dlglk <- fff[[1]]
  hlglk <- fff[[2]]
  out <- matrix(0, 4, 4)
  ll <- upper.tri(out, diag = TRUE)
  out[ll] <- rowMeans(hlglk[ll, ])
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  out <- out + var(t(dlglk))
  out <- -out
  out
}
