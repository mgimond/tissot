# File houses multiple functions (some scraped from the web, others created in
# house) to generate and process tissot indicatrix ellipses.

# The tissot() function generates Tissot indicatrix ellipse and distortion values
# Original Author:  Bill Huber (as sourced at http://gis.stackexchange.com/a/53452)
# Script was slightly modified to deal with some instabilities via numeric.deriv()

tissot <- function(coord, prj=function(z) z+0, asDegrees=TRUE, 
                   A = 6378137, f.inv=298.257223563, idx_area = FALSE, ...) {
  #
  # Compute properties of scale distortion and Tissot's indicatrix at location `x` = c(`lambda`, `phi`)
  # using `prj` as the projection.  `A` is the ellipsoidal semi-major axis (in meters) and `f.inv` is
  # the inverse flattening.  The projection must return a vector (x, y) when given a vector (lambda, phi).
  # (Not vectorized.)  Optional arguments `...` are passed to `prj`.
  #
  # Source: "Map projections: a working manual", Snyder pp 20-26 (http://pubs.er.usgs.gov/publication/pp1395)
  #         (WGS 84 defaults for the ellipsoidal parameters).

  # Split coordinates into separate values
  lambda <- coord[1]
  phi <- coord[2]

  # All input and output angles are in degrees.  
  to.degrees <- function(x) x * 180 / pi
  to.radians <- function(x) x * pi / 180
  clamp <- function(x) min(max(x, -1), 1)                             # Avoids invalid args to asin
  norm <- function(x) sqrt(sum(x*x))

  # Precomputation.
  if (f.inv==0) {                                                     # Use f.inv==0 to indicate a sphere
    e2 <- 0 
  } else {
    e2 <- (2 - 1/f.inv) / f.inv                                       # Squared eccentricity
  }
  if (asDegrees) phi.r <- to.radians(phi) else phi.r <- phi
  cos.phi <- cos(phi.r)                                               # Convenience term
  e2.sinphi <- 1 - e2 * sin(phi.r)^2                                  # Convenience term
  e2.sinphi2 <- sqrt(e2.sinphi)                                       # Convenience term
  if (asDegrees) units <- 180 / pi else units <- 1                    # Angle measurement units per radian
  #
  # Lengths (the metric).
  #
  radius.meridian <- A * (1 - e2) / e2.sinphi2^3                      # (4-18) radius of curvature along meridian
  length.meridian <- radius.meridian                                  # (4-19) Length of a radian
  radius.normal <- A / e2.sinphi2                                     # (4-20) radius of curvature perpendicular to meridian
  length.normal <- radius.normal * cos.phi                            # (4-21)
  #
  # The projection and its first derivatives, normalized to unit lengths.
  #
  x <- c(lambda, phi)
 # d1 <- numericDeriv(quote(prj(x, ...)), theta="x")
  d <- numeric.deriv(expr=prj2(lambda,phi, ...), theta=c("lambda","phi"))
  z <- d[1:2]                                                         # Projected coordinates
  names(z) <- c("x", "y")
  g <- attr(d, "gradient")                                            # First derivative (matrix)
  g <- g %*% diag(units / c(length.normal, length.meridian))          # Unit derivatives
  dimnames(g) <- list(c("x", "y"), c("lambda", "phi"))
  g.det <- det(g)                                                     # Equivalent to (4-15)
  if(idx_area == TRUE){
    return(g.det)
  }
  #
  # Computation.
  #
  h   <- norm(g[, "phi"])                                             # (4-27)
  k   <- norm(g[, "lambda"])                                          # (4-28)
  a.p <- sqrt(max(0, h^2 + k^2 + 2 * g.det))                          # (4-12) (intermediate)
  b.p <- sqrt(max(0, h^2 + k^2 - 2 * g.det))                          # (4-13) (intermediate)
  a   <- (a.p + b.p)/2                                                # (4-12a) semi-major axis
  b   <- (a.p - b.p)/2                                                # (4-13a) semi-minor axis
  omega <- 2 * asin(clamp(b.p / a.p))                                 # (4-1a) max deviation from correct angle
  theta.p <- asin(clamp(g.det / (h * k)))                             # (4-14)
  conv <- (atan2(g["y", "phi"], g["x","phi"]) + pi / 2) %% (2 * pi) - pi # Middle of p. 21
  #
  # The indicatrix itself.
  # `svd` essentially redoes the preceding computation of `h`, `k`, and `theta.p`.
  #
  m <- svd(g)
  axes <- zapsmall(diag(m$d) %*% apply(m$v, 1, function(x) x / norm(x))) # major and minor axes distortion
  dimnames(axes) <- list(c("major", "minor"), NULL)
  
  return(list(location=c(lambda, phi), projected=z, 
              meridian_radius=radius.meridian, meridian_length=length.meridian,
              normal_radius=radius.normal, normal_length=length.normal,
              scale.meridian=h, scale.parallel=k, scale.area=g.det, max.scale=a, min.scale=b, 
              to.degrees(zapsmall(c(angle_deformation=omega, convergence=conv, intersection_angle=theta.p))),
              axes=axes, derivatives=g))
}

# Following function reprocesses the output of `tissot` into convenient 
# geometrical data.
# Original Author:  Bill Huber (as sourced at http://gis.stackexchange.com/a/53452)

indicatrix <- function(x, scale=1, ...) {
  o <- x$projected
  base <- ellipse(o, matrix(c(1,0,0,1), 2), scale=scale, ...)             # A reference circle
  outline <- ellipse(o, x$axes, scale=scale, ...)
  axis.major <- rbind(o + scale * x$axes[1, ], o - scale * x$axes[1, ])
  axis.minor <- rbind(o + scale * x$axes[2, ], o - scale * x$axes[2, ])
  d.lambda <- rbind(o + scale * x$derivatives[, "lambda"], o - scale * x$derivatives[, "lambda"])
  d.phi <- rbind(o + scale * x$derivatives[, "phi"], o - scale * x$derivatives[, "phi"])
  return(list(center=x$projected, base=base, outline=outline, 
              axis.major=axis.major, axis.minor=axis.minor,
              d.lambda=d.lambda, d.phi=d.phi))
}

ellipse <- function(center, axes, scale=1, n=36, from=0, to=2*pi) {
  # Vector representation of an ellipse at `center` with axes in the *rows* of `axes`.
  # Returns an `n` by 2 array of points, one per row.
  theta <- seq(from=from, to=to, length.out=n)
  t((scale * t(axes))  %*% rbind(cos(theta), sin(theta)) + center)
}

# NOTE: The prj() function is substituted for prj2(). It's currently
# not used in any of the  other function
# prj <- function(z, proj.in = CRS("+proj=longlat +datum=WGS84"),
#                 proj.out) {
#   z.pt <- SpatialPoints(coords=matrix(z, ncol=2), proj4string=proj.in)
#   w.pt <- spTransform(z.pt, CRS=proj.out)
#   return(w.pt@coords[1, ])
# }

# Project x,y coord values
prj2 <- function(x, y=NULL, proj.in = 4326, proj.out) {
  names1 <-  NULL
  if(is.null(y) & length(x) == 2){
    x1 <- x
    y <- x1[2]
    x <- x1[1]
    names1 <- names(x1)
  }
  z.pt <-  st_sfc( st_multipoint( matrix(c(x,y), ncol=2)), crs=proj.in)
  w.pt <- sf::st_transform(z.pt, proj.out)
 # out <- w.pt@coords[1, ]
  out <- sf::st_coordinates(w.pt)[,1:2]
  if (!is.null(names1) ){
    names(out) <- names1
  }
  return(out)
}

# Following function was added by Manny Gimond. The default precision in the 
# numericDeriv() function forced some instabilities in the gradient output 
# (e.g. try c(1-60,-40) in a "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84"  
# coordinate system). This function reduces the precision to alleviate the 
# aforementioned  problem.
numeric.deriv <- function(expr, theta, rho=sys.frame(sys.parent()))
{
  # eps <- sqrt(.Machine$double.eps)
  eps <- (.Machine$double.eps)^0.3
  ans <- eval(substitute(expr), rho)
  grad <- matrix(, length(ans), length(theta),
                 dimnames=list(NULL, theta))
  
  for (i in seq_along(theta)) {
    old <- get(theta[i], envir=rho) 
    delta <- eps * max(1, abs(old)) 
    assign(theta[i], old+delta, envir=rho) 
    ans1 <- eval(substitute(expr), rho) 
    assign(theta[i], old, envir=rho) 
    grad[, i] <- (ans1 - ans)/delta 
  }
  attr(ans, "gradient") <- grad
  ans
}

# Run the tissot function and extract its parameters
ti  <- function(coord, prj, proj.in = "+proj=longlat +datum=WGS84" , 
                proj.out, scale=500000, n = 61){
  TI <- tissot(coord, prj, 
               proj.in = proj.in , 
               proj.out = proj.out)
  i <- indicatrix(TI, scale= scale, n=n)
  # Rounding may be needed for values approaching 0
  i <- lapply(i, round, digits=4)
  i[["axes"]] <- TI$axes
  i[["area"]] <- TI$scale.area
  return(i)
}

# Extract tissot area values
ti_area  <- function(coord, prj, proj.in = "+proj=longlat +datum=WGS84" , 
                     proj.out, idx_area = TRUE){
  i <- tissot(coord, prj, 
              proj.in = proj.in , 
              proj.out = proj.out, 
              idx_area = TRUE)
  return(i)
}

# Create an sfc poly from point data
tis2poly <- function(poly, proj.out){
  p <- st_multipolygon(poly)
  return(st_sfc(p, crs = proj.out))
}

# Create an sfc poyline from point data
tis2line <- function(polyline, proj.out){
  p <- lapply(polyline, function(x) st_linestring(matrix(x$poly[1:2], ncol=2)))
  return(st_sfc(p, crs = proj.out))
}

# Generate a local Indicatrix circle
local_TI <- function(long, lat, proj.in = "+proj=longlat +datum=WGS84", 
                     proj.out){
  TI <- tissot(c(long, lat), prj, proj.in=proj.in, proj.out=proj.out)
  TI$max.scale
  TI$min.scale
  TI$max.scale / TI$min.scale
  
  i <- indicatrix(TI, scale=10^4, n=71)
  plot(i$outline, type="n", asp=1, xlab="x", ylab="y")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
  grid(col="white")
  polygon(i$base, col="bisque", border="Gray")
  if(TI$min.scale != TI$max.scale){
    lines(i$axis.major, lwd=2, col=rgb(.45, 1, .45))
    lines(i$axis.minor, lwd=2, col=rgb(1, .35, .35))
  }
  lines(i$outline, asp=1, lwd=2, col=rgb(.35, .35, 1))
  lines(i$d.lambda, lwd=2, col=rgb(.35, .35, .35), lty=3)
  lines(i$d.phi, lwd=2, col=rgb(.35, .35, .35), lty=2)
  title.out<- sprintf("Long=%0.2f, Lat=%0.3f, Area scale=%0.3f, 
                      Max scale=%0.3f, Min scale=%0.3f\n %s ",
                      long,lat,TI$scale.area, TI$max.scale, 
                      TI$min.scale, proj.out)
  mtext(title.out,3,col="blue")
  return(TI)
}

# Generate sf objects
tissot_sf <- function(i.lst, proj.out){
  bases <- lapply(i.lst, function(x) list(poly = x$base))
  base_sf <- tis2poly(bases, proj.out)
  
  # Extract indicatrix
  ind <- lapply(i.lst, function(x) list(poly = x$outline))
  ind_sf <- tis2poly(ind, proj.out)
  ind_sf <- st_cast(ind_sf, "POLYGON")
  #st_is_valid(ind_sf)
  
  # Extract minor/major axes
  maja <- lapply(i.lst, function(x) st_linestring(x$axis.major))
  maja_sf <- st_sfc(maja, crs = proj.out)
  
  mina <- lapply(i.lst, function(x) st_linestring(x$axis.minor))
  mina_sf <- st_sfc(mina, crs = proj.out)
  
  # Extract lambda and phi
  lam <- lapply(i.lst, function(x) st_linestring(x$d.lambda))
  lam_sf <- st_sfc(lam, crs = proj.out)
  
  phi <- lapply(i.lst, function(x) st_linestring(x$d.phi))
  phi_sf <- st_sfc(phi, crs = proj.out)
  
  # Extract area values
  ind_sf <- st_as_sf(ind_sf)
  ind_sf$area <- sapply(i.lst, function(x) (x$area))
  
  # Output list
  list(base = base_sf, ind = ind_sf , maja = maja_sf, mina = mina_sf,
       lam = lam_sf, phi = phi_sf, area = ind_sf$area)
}


# Function to recenter and dissolve a map along a different longitude
# Code originally written by Michael Summner 
# (https://github.com/AustralianAntarcticDivision/SOmap/issues/34)
c_recenter <- function(x, clon = NULL, ..., tryfix = TRUE) {
  if (is.null(clon)) return(x)
  if (!st_is_longlat(x)) stop("recentring not appropriate for non longlat data")
  ## save the crs while we do our munging
  crs <- st_crs(x)
  x <- st_set_crs(x, NA)
  
  
  ## try to fix problematic geometry
  if (tryfix) {
    if (all(grepl("POLYGON", st_geometry_type(x)))) x <- suppressWarnings(st_buffer(sf::st_as_sf(x), 0))
    x <- st_wrap_dateline(x, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
  }
  wbox <- st_bbox(c(xmin = -180, ymin = -90, xmax = (clon)%%360 - 180, ymax = 90))
  west <- suppressWarnings(st_crop(x, wbox))
  west <- st_set_geometry(west, st_geometry(west) + c(360, 0))
  east <- suppressWarnings(st_crop(x, st_bbox(c(xmin = (clon)%%360 - 180, xmax = 180, ymin = -90, ymax = 90))))
  xx <- rbind(
    west, east
  ) 
  ## ensure geometries are of consistent type
  xx <- sf::st_cast(xx)
  bb <- st_bbox(xx)
  st_set_crs(xx, crs)
}


# Check for valid coordinate extent. If projected coord values are Inf, 
# remove that coordinate value from the returned list.

coord_check <- function(lonlat, proj.in = "+proj=latlong +ellps=WGS84",
                        proj.out){
  coord.prj <- sf::sf_project(lonlat, from = proj.in, to = proj.out)
  coord.chk <- t(apply(coord.prj[,1:2], 1, function(x) !is.infinite(x)))
  coord.bool <- apply(coord.chk, 1, function(x) as.logical(min(x)))
  lonlat[coord.bool, ]
}
