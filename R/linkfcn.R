##' Link function for the exponential family.
##'
##' \code{linkfcn} maps the mean of the response variable \code{mu} to
##' the linear predictor \code{z}. \code{linkinv} is its inverse.
##'
##' For the Gaussian family, if the link parameter is positive, then
##' the extended link is used, defined by \deqn{z =
##' \frac{sign(\mu)|\mu|^\nu - 1}{\nu}}{z = (sign(mu)*abs(mu)^nu -
##' 1)/nu} In the other case, the usual Box-Cox link is used.
##'
## For the Poisson and gamma families, the Box-Cox transformation is
## used, defined by \deqn{z = \frac{\mu^\nu - 1}{\nu}}{z = (mu^nu -
## 1)/nu}
##'
##' For the Poisson and gamma families, if the link parameter is
##' positive, then the link is defined by \deqn{z = \frac{sign(w)
##' (e^{\nu |w|}-1)}{\nu}}{z = sign(w)*expm1(nu*w)/nu} where
##' \eqn{w = \log(\mu)}{w = log(mu)}. In the other case, the usual
##' Box-Cox link is used.
##'
##' For the GEV binomial family, the link function is defined by
##' \deqn{\mu = 1 - \exp\{-\max(0, 1 + \nu z)^{\frac{1}{\nu}}\}}{mu =
##' 1 - exp[-max(0, 1 + nu z)^(1/nu)]} for any real \eqn{\nu}{nu}. At
##' \eqn{\nu = 0}{nu = 0} it reduces to the complementary log-log
##' link.
##'
##' The Wallace binomial family is a fast approximation to the robit
##' family. It is defined as \deqn{\mu =
##' \Phi(\mbox{sign}(z) c(\nu) \sqrt{\nu \log(1 + z^2/\nu)})}{mu =
##' Phi(sign(z) c(nu) sqrt{nu log(1 + z^2/nu)})}
##' where \eqn{c(\nu) = (8\nu+1)/(8\nu+3)}{c(nu) = (8*nu+1)/(8*nu+3)}
##'
##' @title Calculate the link function for exponential families
##' @param mu Numeric. The mean of the response variable.
##' @param z Numeric. The linear predictor.
##' @param linkp The link function parameter. A scalar.
##' @param family The distribution of the response variable from
##'   \code{\link{.geoBayes_models}}. Either an integer or the
##'   family name.
##' @return A numeric array of the same dimension as the function's
##' first argument.
## ##' @seealso \code{\link{comparebinlinks}}
##' @examples \dontrun{
##' mu <- seq(0.1, 0.9, 0.1)
##' linkfcn(mu, 7, "binomial")       # robit(7) link function
##' linkfcn(mu, , "binomial.logit")  # logit link function
##'
##' mu <- seq(-3, 3, 1)
##' linkfcn(mu, 0.5, "gaussian")     # sqrt transformation
##' linkinv(linkfcn(mu, 0.5, "gaussian"), 0.5, "gaussian")
##' curve(linkfcn(x, 0.5, "gaussian"), -3, 3)
##' }
##' @family linkfcn
##' @name linkfcn
##' @rdname linkfcn
##' @importFrom stats qlogis qnorm qt qf plogis pnorm pt
##' @useDynLib geoBayes flinkfcn
##' @references Evangelou, E., & Roy, V. (2019). Estimation and prediction for
##'   spatial generalized linear mixed models with parametric links
##'   via reparameterized importance sampling. Spatial Statistics, 29,
##'   289-315.
##' @export
linkfcn <- function (mu, linkp, family = "gaussian") {
  if (is.numeric(family)) {
    ifam <- as.integer(family)
    if (ifam < 0 || ifam > NROW(.geoBayes_models)) {
      stop ("Cannot deduce family")
    }
  } else {
    ifam <- .geoBayes_family(family)
  }
  if (!.geoBayes_models$needlinkp[ifam]) linkp <- 0
  if (length(linkp) != 1) {
    stop ("The linkp argument must be scalar.")
  }
  linkp <- .geoBayes_getlinkp(linkp, ifam)
  x <- y <- as.double(mu)
  ii <- is.finite(x)
  n <- as.integer(sum(ii))
  mu[ii] <- .Fortran("flinkfcn", x[ii], n, y[ii], linkp, ifam,
                     PACKAGE = "geoBayes")[[1]]
  mu
}

##' @family linkfcn
##' @rdname linkfcn
##' @useDynLib geoBayes flinkinv
##' @export
linkinv <- function (z, linkp, family = "gaussian") {
  if (is.numeric(family)) {
    ifam <- as.integer(family)
    if (ifam < 0 || ifam > NROW(.geoBayes_models)) {
      stop ("Cannot deduce family")
    }
  } else {
    ifam <- .geoBayes_family(family)
  }
  if (!.geoBayes_models$needlinkp[ifam]) linkp <- 0
  if (length(linkp) != 1) {
    stop ("The linkp argument must be scalar.")
  }
  linkp <- .geoBayes_getlinkp(linkp, ifam)
  x <- y <- as.double(z)
  ii <- is.finite(x)
  n <- as.integer(sum(ii))
  z[ii] <- .Fortran("flinkinv", x[ii], n, y[ii], linkp, ifam,
                    PACKAGE = "geoBayes")[[1]]
  z
}
