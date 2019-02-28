######################################################################
##
### Commentary: MCMC analysis of the rhizoctonia data.
##
######################################################################

library(geoBayes)

data(rhizoctonia)

### Create prediction grid
predgrid <- mkpredgrid2d(rhizoctonia[c("Xcoord", "Ycoord")],
                         par.x = 100, chull = TRUE, exf = 1.2)
rhizdata <- stackdata(rhizoctonia, predgrid$grid)

### Define the model
corrf <- "spherical"
kappa <- 0
ssqdf <- 1
ssqsc <- 1
betm0 <- 0
betQ0 <- .01
phiprior <- c(100, 1, 1000, 100) # U(100, 200)
phisc <- 2.5
omgprior <- c(2, 1, 1, 0)        # Exp(mean = 2)
omgsc <- .3
linkp <- 20

### MCMC sizes
Nout <- 1000
Nthin <- 10
Nbi <- 300

emt <- mcsglmm(Infected ~ 1, 'binomial', rhizdata, weights = Total,
               atsample = ~ Xcoord + Ycoord,
               Nout = Nout, Nthin = Nthin, Nbi = Nbi,
               betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
               corrprior = list(phi = phiprior, omg = omgprior), linkp = linkp,
               corrfcn = corrf, kappa = kappa,
               corrtuning = list(phi = phisc, omg = omgsc, kappa = 0), 
               dispersion = 1, test = TRUE)

emc <- update(emt, test = FALSE)

emcmc <- mcmcmake(emc)

plot.ts(cbind(phi = emc$MCMC$phi, omg = emc$MCMC$omg,
              beta = c(emc$MCMC$beta), ssq = emc$MCMC$ssq),
        nc = 2)

summary(emcmc[, c("phi", "omg", "beta", "ssq")])

plot(emcmc[, c("phi", "omg", "beta", "ssq")])

library(geoR)
z0pred <- rowMeans(emc$MCMC$z[!emc$MCMC$whichobs, ])
geoR:::image.kriging(locations = predgrid$xygrid,
                     borders = predgrid$borders, values = z0pred,
                     x.leg = c(3150, 3450), y.leg = c(1120, 1200),
                     col = gray((64:32)/64))

