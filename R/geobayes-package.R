##' Analysis of geostatistical data using Bayes and Empirical Bayes
##' methods.
##'
##' This package provides functions to fit geostatistical data. The
##' data can be continuous, binary or count data and the models
##' implemented are flexible. Conjugate priors are assumed on some
##' parameters while inference on the other parameters can be done
##' through a full Bayesian analysis of by empirical Bayes methods.
##'
##' Some demonstration examples are provided. Type \code{demo(package
##' = "geoBayes")} to examine them.
##' @title The \code{geoBayes} package
##' @name geoBayes
##' @aliases geoBayes-package
## @docType package
##' @useDynLib geoBayes, .registration = TRUE
##' @author Evangelos Evangelou <e.evangelou@@maths.bath.ac.uk> and
##' Vivekananda Roy <vroy@@iastate.edu>
##' @examples \dontrun{
##' demo(package = "geoBayes")
##' demo(rhizoctonia1, package = "geoBayes")
##' demo(rhizoctonia1, package = "geoBayes")
##' }
##' @references Roy, V., Evangelou, E. and Zhu, Z. (2014). Empirical
##' Bayes methods for the transformed Gaussian random fields model
##' with additive measurement errors. In Upadhyay, S. K., Singh, U.,
##' Dey, D. K., and Loganathan, A., editors, \emph{Current Trends in
##' Bayesian Methodology with Applications}, Boca Raton, FL, USA, CRC
##' Press.
##'
##' Roy, V., Evangelou, E., and Zhu, Z. (2015). Efficient estimation
##' and prediction for the Bayesian spatial generalized linear mixed
##' model with flexible link functions. \emph{Biometrics}, 72(1),
##'   289-298.
##'
##' Evangelou, E., & Roy, V. (2019). Estimation and prediction for
##'   spatial generalized linear mixed models with parametric links
##'   via reparameterized importance sampling. \emph{Spatial Statistics}, 29,
##'   289-315.
##'
##' Roy, V., & Evangelou, E. (2024). Selection of proposal
##'   distributions for multiple importance sampling.
##'   \emph{Statistica Sinica}, 34, 27-46.
##' @keywords package
"_PACKAGE"


.onUnload <- function (libpath) {
  library.dynam.unload("geoBayes", libpath)
}
