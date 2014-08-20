##' Combine \code{data.frame}s
##'
##' This function combines \code{data.frame}s by filling in missing
##' variables with \code{NA}. This is useful for combining data from
##' sampled locations with prediction locations.
##' @title Combine \code{data.frame}s
##' @param ... \code{data.frame}s or objects that can be coerced to
##' \code{data.frame}s 
##' @return A stacked \code{data.frame}
##' @export 
##' @examples
##' \dontrun{
##'   data(rhizoctonia)
##' 
##'   predgrid <- mkpredgrid2d(rhizoctonia[c("Xcoord", "Ycoord")],
##'                            par.x = 100, chull = TRUE, exf = 1.2)
##'   rhizdata <- stackdata(rhizoctonia, predgrid$grid)
##' }
stackdata <- function (...) {
  fillNA <- function (d, allnames) {
    ii <- !(allnames %in% names(d))
    d[allnames[ii]] <- NA
    d
  }
  input <- lapply(list(...), data.frame)
  nmall <- unique(unlist(lapply(input, names)))
  newdata <- lapply(input, fillNA, allnames = nmall)
  out <- do.call(rbind, newdata)
  out
}

##' This function creates a grid for prediction.
##'
##' If \code{chull} this function first calculates the convex hull of
##' the points. If \code{exf} is not 1 the borders are expanded. Then
##' the function calls \code{\link[sp]{point.in.polygon}} to select
##' points that fall inside the borders.
##' @title Make prediction grid
##' @param pnts.x x coordinate of the domain. Could also be a
##' two-column matrix containing the x and y coordinates
##' @param pnts.y y coordinate of the domain. Should be
##' omitted or set to \code{NULL} if the argument \code{pnts.x} is a
##' two-column matrix.
##' @param par.x A scalar parameter for the x component of the new
##' grid. This parameter corresponds to either the \code{by} or the
##' \code{length.out} arguments of the function
##' \code{\link[base]{seq}}. Could also be a vector of two elements
##' containing the parameter for x and y.
##' @param par.y As in \code{par.x} for the y component of the new
##' grid. Should be omitted or set to \code{NULL} if the argument
##' \code{par.x} is a two-dimensional vector.
##' @param isby If \code{TRUE}, the arguments \code{par.x} and
##' \code{par.y} correspond to the \code{by} argument of the function
##' \code{\link[base]{seq}}, otherwise they correspond to
##' \code{length.out}.
##' @param chull Whether to calculate the convex hull of the points.
##' Set this to \code{TRUE} if \code{pnts.x} and \code{pnts.y} denote
##' the sampled locations. If they correspond to the borders of the
##' domain, it is recommended to set this to \code{FALSE}.
##' @param exf An expansion factor of the convex hull of
##' \code{cbind(pnts.x, pnts.y)}. Must be positive. If larger or
##' smaller than 1, the convex hull is respectively expanded or
##' contracted.
##' @return A list with components
##' \itemize{
##' \item \code{grid} A two-column matrix with the prediction grid
##' \item \code{xycoord} A list with components "x" and "y"
##' containing the sequence of points used to create the grid
##' \item \code{xygrid} A matrix with the full square grid derived
##' from \code{xycoord}
##' \item \code{borders} The (expanded) borders of the domain
##' }
##' @seealso \code{\link[geoR]{pred_grid}}
##' @importFrom sp point.in.polygon
##' @export
##' @examples
##' \dontrun{
##'   data(rhizoctonia)
##' 
##'   predgrid <- mkpredgrid2d(rhizoctonia[c("Xcoord", "Ycoord")],
##'                            par.x = 100, chull = TRUE, exf = 1.2)
##' }
mkpredgrid2d <- function (pnts.x, pnts.y, par.x, par.y, isby = FALSE,
                          chull = FALSE, exf = 1) {
  if (exf <= 0) stop ("Argument exf must be positive")
  if (missing(pnts.y)) pnts.y <- NULL
  if (missing(par.y)) par.y <- NULL
  if (!is.null(pnts.x)) pnts.x <- as.matrix(pnts.x)
  if (!is.null(pnts.y)) pnts.y <- as.matrix(pnts.y)
  ph <- cbind(pnts.x, pnts.y)
  nm <- dimnames(ph)[[2]]
  par <- c(par.x[1], par.y[1])
  d <- NCOL(ph)
  if (d != 2) stop ("Can only generate 2-dimensional grids")
  if (isTRUE(chull)) {
    ph <- ph[chull(ph), ]
  }
  if (exf != 1) { ## Extend the covex hull by a factor exf
    centr <- colMeans(ph)     # Central coordinate of the polygon
    phu <- cbind(ph[, 1] - centr[1], ph[, 2] - centr[2]) # Uncenter
    phu_r <- exf*sqrt(phu[, 2]^2 + phu[, 1]^2) # Covert to polar
    phu_u <- atan2(phu[, 2], phu[, 1])
    phu_x <- phu_r*cos(phu_u) # Convert back to cartesian
    phu_y <- phu_r*sin(phu_u)
    ph[, 1] <- phu_x + centr[1]
    ph[, 2] <- phu_y + centr[2]
  }
  ft <- apply(ph, 2, range)
  if (isTRUE(!isby)) {
    par <- ceiling(par)
    par <- ((ft[2, ] - ft[1, ])/(par - 1))
  } else {
    par <- rep_len(par, d)
  }
  xycoord <- lapply(1:d, function (i) seq(ft[1, i], ft[2, i], par[i]))
  names(xycoord) <- nm
  eg <- as.matrix(expand.grid(xycoord, KEEP.OUT.ATTRS = FALSE))
  dimnames(eg) <- list(NULL, nm)
  iin <- sp::point.in.polygon(eg[, 1], eg[, 2], ph[, 1], ph[, 2]) > 0
  grid <- eg[iin, ]
  out <- list(grid = grid, xycoord = xycoord, xygrid  = eg, borders = ph)
  out
}
