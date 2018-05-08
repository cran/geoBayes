##' Standard errors for BF estimates
##'
##' Using the formula from Casella.
##' @title SE for BF
##' @param mcrun The output from the function \code{mcsglmm}.
##' @param transf The type of transformation to use.
##' @return The precision matrix
##' @importFrom stats var
##' @useDynLib geoBayes llikfcn_dh_tr
##' @export
bf2se <- function(mcrun, transf = c("no", "mu", "wo"))
{
  family <- mcrun$family
  corrfcn <- mcrun$corrfcn
  F <- mcrun$modelmatrix
  n <- NROW(F)
  p <- NCOL(F)
  loc <- mcrun$locations
  longlat <- mcrun$longlat
  dm <- spDists(loc, longlat = longlat)
  y <- mcrun$response
  l <- mcrun$weights
  linkp <- mcrun$nu[1]
  phi <- mcrun$phi[1]
  omg <- mcrun$omg[1]
  kappa <- mcrun$kappa[1]
  betm0 <- mcrun$betm0
  betQ0 <- mcrun$betQ0
  ssqdf <- mcrun$ssqdf
  ssqsc <- mcrun$ssqsc
  dispersion <- mcrun$dispersion
  tsqdf <- mcrun$tsqdf
  tsqsc <- mcrun$tsqsc
  icf <- .geoBayes_correlation(corrfcn)
  ifam <- .geoBayes_family(family)
  nu <- .geoBayes_getlinkp(linkp, ifam)
  if (ifam == 0) {
    tsq <- tsqsc
  } else {
    tsq <- dispersion
    tsqdf <- 0
  }
  getsample <- transfsample(list(mcrun[c("z", "mu", "nu", "whichobs")]),
                            mcrun[c("response", "family")], transf)
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
