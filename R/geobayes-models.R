#' Models used in the geoBayes package
#'
#' This hidden variable contains a choice of models that can be fit
#' with this package. The model can be chosen in the \code{family}
#' input of the relevant function. This variable cannot be
#' overwritten.
#' @name geoBayes_models
#' @export
.geoBayes_models <-
  data.frame(
    family = c("gaussian",            # 1
               "binomial.robit",      # 2
               "binomial.logit",      # 3
               "binomial.probit",     # 4
               "binomial.Wallace",    # 5
               "poisson.modifiedbc",  # 6
               "poisson.boxcox",      # 7
               "Gamma.modifiedbc",    # 8
               "Gamma.boxcox",        # 9
               "binomial.modifiedGEV",# 10
               "binomial.GEV"         # 11
               ),
    needlinkp = c(TRUE,               #  1 "gaussian",
                  TRUE,               #  2 "binomial.robit",
                  FALSE,              #  3 "binomial.logit",
                  FALSE,              #  4 "binomial.probit",
                  TRUE,               #  5 "binomial.Wallace",
                  TRUE,               #  6 "poisson.modifiedbc",
                  TRUE,               #  7 "poisson.boxcox",
                  TRUE,               #  8 "Gamma.modifiedbc",
                  TRUE,               #  9 "Gamma.boxcox"
                  TRUE,               # 10 "binomial.modifiedGEV"
                  TRUE                # 11 "binomial.GEV"
                  ),
    linkplb = c(-Inf,                 #  1 "gaussian",
                0,                    #  2 "binomial.robit",
                NA,                   #  3 "binomial.logit",
                NA,                   #  4 "binomial.probit",
                0,                    #  5 "binomial.Wallace",
                0,                    #  6 "poisson.modifiedbc",
                -Inf,                 #  7 "poisson.boxcox",
                0,                    #  8 "Gamma.modifiedbc",
                -Inf,                 #  9 "Gamma.boxcox"
                0,                    # 10 "binomial.modifiedGEV"
                -Inf                  # 11 "binomial.GEV"
                ),
    incllb = c(NA,                    #  1 "gaussian",
               FALSE,                 #  2 "binomial.robit",
               NA,                    #  3 "binomial.logit",
               NA,                    #  4 "binomial.probit",
               FALSE,                 #  5 "binomial.Wallace",
               TRUE,                  #  6 "poisson.modifiedbc",
               NA,                    #  7 "poisson.boxcox",
               TRUE,                  #  8 "Gamma.modifiedbc",
               NA,                    #  9 "Gamma.boxcox"
               TRUE,                  # 10 "binomial.modifiedGEV"
               NA                     # 11 "binomial.GEV"
               ),
    linkpub = c(Inf,                  #  1 "gaussian",
                Inf,                  #  2 "binomial.robit",
                NA,                   #  3 "binomial.logit",
                NA,                   #  4 "binomial.probit",
                Inf,                  #  5 "binomial.Wallace",
                Inf,                  #  6 "poisson.modifiedbc",
                Inf,                  #  7 "poisson.boxcox",
                Inf,                  #  8 "Gamma.modifiedbc",
                Inf,                  #  9 "Gamma.boxcox"
                Inf,                  # 10 "binomial.modifiedGEV"
                Inf                   # 11 "binomial.GEV"
                ),
    inclub = c(NA,                    #  1 "gaussian",
               NA,                    #  2 "binomial.robit",
               NA,                    #  3 "binomial.logit",
               NA,                    #  4 "binomial.probit",
               NA,                    #  5 "binomial.Wallace",
               NA,                    #  6 "poisson.modifiedbc",
               NA,                    #  7 "poisson.boxcox",
               NA,                    #  8 "Gamma.modifiedbc",
               NA,                    #  9 "Gamma.boxcox"
               NA,                    # 10 "binomial.modifiedGEV"
               NA                     # 11 "binomial.GEV"
               ),
    haswo = c(FALSE,                  #  1 "gaussian",
              TRUE,                   #  2 "binomial.robit",
              FALSE,                  #  3 "binomial.logit",
              FALSE,                  #  4 "binomial.probit",
              FALSE,                  #  5 "binomial.Wallace",
              FALSE,                  #  6 "poisson.modifiedbc",
              TRUE,                   #  7 "poisson.boxcox",
              FALSE,                  #  8 "Gamma.modifiedbc",
              FALSE,                  #  9 "Gamma.boxcox"
              FALSE,                  # 10 "binomial.modifiedGEV"
              TRUE                    # 11 "binomial.GEV"
              ),
    stringsAsFactors = FALSE)

lockBinding(".geoBayes_models", environment())

.geoBayes_family <- function (family) {
  ## family is a character, partially matching .geoBayes_models before
  ## and after the dot.
  nfamily <- length(family)
  if (nfamily > 1) {
    stop ("You cannot use more than one family.")
  } else if (nfamily == 0) {
    as.integer(1L)
  } else if (pmatch(family, "transformed.gaussian", 0)) {
    as.integer(0L)
  } else {
    fl <- strsplit(family, split = "[[:punct:][:space:]]+")[[1]]
    pattern <- if(length(fl) > 1)
                 paste0("\\.", fl[-1], "[[:alnum:]]*", collapse = "")
    pattern <- paste0("^", fl[1], "[[:alnum:]]*", pattern)
    out <- grep(pattern, .geoBayes_models$family, value = FALSE)
    if (length(out) == 0) stop ("Cannot deduce family from input.")
    out[1]
  }
}



## Checks whether the link parameter is consistent with the family and
## returns a real representation of the link parameter.
.geoBayes_getlinkp <- function (linkp, ifam) {
  lbt <- `>`; ubt <- `<`
  if (ifam == 0 || pmatch(ifam, "transformed.gaussian", 0)) {
    family <- "transformed.gaussian"
    lb <- 0; ub <- Inf; needlinkp <- TRUE
  } else {
    if (!is.numeric(ifam)) ifam <- .geoBayes_family(ifam)
    family <- .geoBayes_models$family[ifam]
    lb <- .geoBayes_models$linkplb[ifam]
    ub <- .geoBayes_models$linkpub[ifam]
    needlinkp <- .geoBayes_models$needlinkp[ifam]
    if (needlinkp && is.finite(lb) && .geoBayes_models$incllb[ifam])
      lbt <- `>=`
    if (needlinkp && is.finite(ub) && .geoBayes_models$inclub[ifam])
      ubt <- `<=`
  }
  if (!needlinkp) return(0.)
  if (is.finite(lb) && !all(lbt(linkp, lb))) {
    stop (paste("Link parameter must be", gsub("[^><=]","",deparse(lbt)), lb,
                "for the", family, "family"))
  }
  if (is.finite(ub) && !all(ubt(linkp, ub))) {
    stop (paste("Link parameter must be", gsub("[^><=]","",deparse(ubt)), ub,
                "for the", family, "family"))
  }
  as.double(linkp)
}

###   valid <- .geoBayes_models_opts$linkp[[family]]
###   for (i in seq_along(valid)) {
###     test <- valid[[i]]
###     nmtest <- names(valid[i])
###     if (length(nmtest) && nmtest != "") { # Character value
###       if (all(as.character(linkp) == nmtest))
###         return (rep(as.double(test), k))
###     } else if (length(test) == 2) { # Interval value
###       if (all(linkp > test[1]) && all(linkp < test[2]))
###         return (as.double(linkp))
###     } else if (length(test) == 1) { # Scalar value
###       return (rep(as.double(test), k))
###     }
###   }
###   stop ("Cannot deduce link parameter from family.")
