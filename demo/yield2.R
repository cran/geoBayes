######################################################################
## 
### Commentary: Estimate the effect of the rhizoctonia disease on
### yield using full Bayesian analysis.
## 
######################################################################

library(geoBayes)

### Load the data
data(rhizoctonia)
rhiz <- na.omit(rhizoctonia)
rhiz$IR <- rhiz$Infected/rhiz$Total # Incidence rate of the
                              # rhizoctonia disease

### Define the model
corrf <- "spherical"
ssqdf <- 4
ssqsc <- 1
tsqdf <- 4
tsqsc <- 1
betm0 <- 0
betQ0 <- diag(.01, 2, 2)
phiprior <- c(200, 1, 1000, 100) # U(100, 300)
phisc <- .7
omgprior <- c(3, 1, 1000, 0) # U(0, 3)
omgsc <- 1
linkp <- 1

## MCMC parameters
Nout <- 1000
Nbi <- 3000
Nthin <- 10

samplt <- mcstrga(Yield ~ IR, data = rhiz, 
                  atsample = ~ Xcoord + Ycoord, corrf = corrf, 
                  Nout = Nout, Nthin = Nthin,
                  Nbi = Nbi, betm0 = betm0, betQ0 = betQ0,
                  ssqdf = ssqdf, ssqsc = ssqsc,
                  tsqdf = tsqdf, tsqsc = tsqsc,
                  corrprior = list(phi = phiprior, omg = omgprior),
                  linkp = linkp,
                  corrtuning = list(phi = phisc, omg = omgsc, kappa = 0), 
                  test=100)

sample <- update(samplt, test = FALSE)

mcsamp <- mcmcmake(sample)

library(mcmcplots)
ipar <- grep(paste(c("phi", "tsq", "ssq", "beta(_[0-9]+)?", "omg"),
                   collapse = "|"), dimnames(mcsamp)[[2]])

traplot(mcsamp[, ipar])
denplot(mcsamp[, ipar])
summary(mcsamp[, ipar])
