###############################################################################
#### Compute SE for expectations
###############################################################################

fse_bm <- function (fsamples, vsamples, Nmcmc, nb1)
### fsamples has dimension Ntot by nh by nf
### vsamples has dimension Ntot by nh
### Nmcmc has length kg and sum(Nmcmc) = Ntot. Nmcmc are the sample sizes
### from each density
### nb1 - integer, the number of batches
{
  vsamples <- as.matrix(vsamples)
  Ntot <- nrow(vsamples)
  nh <- ncol(vsamples)
  dimout <- dimw <- dim(fsamples)
  if (dimw[1] != Ntot) stop ("fsamples and vsamples row mismatch")
  if (dimw[2] != Ntot) stop ("fsamples and vsamples col mismatch")
  if (length(dimw) == 2L) {
    dimw <- dim(fsamples) <- c(dimw, 1L)
  } else {
    dimw <- dim(fsamples) <-
      c(dimw[1:2], prod(dimw[3:length(dimw)]))
  }
  nf <- dimw[3]
  wsamples <- fsamples*c(vsamples)
  
}
