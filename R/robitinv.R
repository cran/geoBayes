##' Compute the inverse logit, probit, robit transformations.
##'
##' The input argument \code{d} determines the link function. If
##' \code{"logit"} or \code{"probit"}, then the respective link
##' function is used. If a positive integer, then the robit(\code{d}),
##' i.e. the CDF of the Student's t with degrees of freedom \code{d}
##' is used. Note that, contrary to the common, the inverse logit
##' transformation is defined here as the CDF of the logistic
##' distribution with scale \code{0.6458} to make it comparable with
##' the robit link.
##' @title Compute the inverse logit, probit, robit transformations
##' @param x Numerical input at which the inverse link function is to
##' be computed.
##' @param d Defines the link function. See details.
##' @param log.p Whether to return the \code{log} of the inverse link.
##' @return A numerical output of the same lenght and dimension as
##' \code{x}.
##' @export
robitinv <- function (x,d,log.p=TRUE) {
  links <- c('logit','probit')
  stopmsg <- paste('Link function must be either a positive number or',
                   paste(paste('"',links,'"',sep=''),collapse=', '))
  if (is.character(d)) {
    i <- pmatch(d,links,nomatch=length(links)+1)
    out <- switch(i,
                  plogis(x,scale=0.6458,log.p=log.p), # logit
                  pnorm(x,log.p=log.p),               # probit
                  stop (stopmsg))                     # error
  } else if (d <= 0) {
    stop (stopmsg)
  } else {
    out <- pt(x,d,log.p=log.p)
  }
  return(out)
}
