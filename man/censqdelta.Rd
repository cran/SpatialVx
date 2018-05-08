\name{censqdelta}
\alias{censqdelta}
\title{
Centered on Square Domain Baddeley's Delta Metric
}
\description{
Baddeley's delta metric is sensitive to the position of non-zero grid points within the domain, as well as to the size of the domain.  In order to obtain consistent values of the metric across cases, it is recommended to first position the sets to be compared so that they are centered with respect to one another on a square domain; where the square domain is the same for all comparison sets.
}
\usage{
censqdelta(x, y, N, const = Inf, p = 2, ...)
}
\arguments{
  \item{x, y}{ Matrices representing binary images to be compared.  If they are not binary, then they will be forced to binary by setting anything above zero to one.
}
  \item{N}{
The size of the square domain.  If missing, it will be the size of the largest side, and if it is even, one will be added to it.
}
  \item{const}{ single numeric giving the \code{c} argument to \code{deltametric}, which is the constant value over which the distance map is reduced to this value.}
  \item{p}{ single numeric giving the \code{p} argument to \code{deltametric}, which specifies the type of Lp norm used to calculate the delta maetric.}
  \item{\dots}{ Not used.
}
}
\details{
Baddeley's delta metric (Baddeley, 1992a,b) is the L_p norm over the absolute difference of distance maps for two binary images, A and B.  A concave function (e.g., f(t) = min(t, constant)) may first be applied ot each distance map before taking their absolute differences, which makes the result less sensitive to small changes in one or both images than other similar metrics.  The metric is sensitive to size, shape and location differences, which make it very practical for comparing forecasts to observations in terms of the position, area extent, and area shape  errors.  However, its sensitivity to domain size and position within the domain are undesirable, but are easily fixed by calculating the metric over a consistent, square domain with the combined verification set centerd on that domain.  See the example section below to see the issue.

This function essentially takes a window of size N by N and moves so that the centroid of each pair of sets, A and B, is the center of the window before calculating the metric.

Centering and squaring is recommended for carrying out a procedure such as that proposed in Gilleland et al. (2008).  Centering on a square domain alleviates the problems discovered by Schwedler and Baldwin (2011) who suggested using a small value of the constant in f(t) = min(t, constant) applied to the distance maps.  This solution is not very appealing because of the sensitivity in choice of the constant that generally diminishes as it approaches the domain size (Gilleland, 2011).

After centering the sets on a square domain, the function \code{deltametric} from package \pkg{spatstat} is used to calculate the metric.
}
\value{
A single numeric value is returned.
}
\references{
Baddeley, A. (1992a)  An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, W. Forstner and S. Ruwiedel, Eds., Wichmann, 59--78.

Baddeley, A. (1992b)  Errors in binary images and an Lp version of the Hausdorff metric.  \emph{Nieuw Arch. Wiskunde}, \bold{10}, 157--183.

Gilleland, E. (2011) Spatial Forecast Verification: Baddeley's Delta Metric Applied to the ICP Test Cases. \emph{Weather Forecast.}, \bold{26} (3), 409--415.

Gilleland, E. (2017) A new characterization in the spatial verification framework for false alarms, misses, and overall patterns. \emph{Weather Forecast.}, \bold{32} (1), 187--198, DOI: 10.1175/WAF-D-16-0134.1.

Gilleland, E., Lee, T. C. M.,  Halley Gotway, J., Bullock, R. G. and Brown, B. G. (2008) Computationally efficient spatial forecast verification using Baddeley's delta image metric.  \emph{Mon. Wea. Rev.}, \bold{136}, 1747--1757.

Schwedler, B. R. J. and Baldwin, M. E. (2011) Diagnosing the sensitivity of binary image measures to bias, location, and event frequency within a forecast verification framework. \emph{Weather Forecast.}, \bold{26}, 1032--1044.

}
\author{
Eric Gilleland
}
\seealso{
\code{\link{deltametric}}, \code{\link{locmeasures2d} }
}
\examples{
x <- y <- matrix( 0, 100, 200 )
x[ 45, 10 ] <- 1
x <- kernel2dsmooth( x, kernel.type = "disk", r = 4 )

y[ 50, 60 ] <- 1
y <- kernel2dsmooth( y, kernel.type = "disk", r = 10 )

censqdelta( x, y )

\dontrun{
# Example form Gilleland (2017).
#
# I1 = circle with radius = 20 centered at 100, 100
# I2 = circle with radius = 20 centered at 140, 100
# I3 = circle with radius = 20 centered at 180, 100
# I4 = circle with radius = 20 centered at 140, 140

I1 <- I2 <- I3 <- I4 <- matrix( 0, 200, 200 )

I1[ 100, 100 ] <- 1
I1 <- kernel2dsmooth( I1, kernel.type = "disk", r = 20 )
I1[ I1 > 0 ] <- 1
if( any( I1 < 0 ) ) I1[ I1 < 0 ] <- 0

I2[ 140, 100 ] <- 1
I2 <- kernel2dsmooth( I2, kernel.type = "disk", r = 20 )
I2[ I2 > 0 ] <- 1
if( any( I2 < 0 ) ) I2[ I2 < 0 ] <- 0

I3[ 180, 100 ] <- 1
I3 <- kernel2dsmooth( I3, kernel.type = "disk", r = 20 )
I3[ I3 > 0 ] <- 1
if( any( I3 < 0 ) ) I3[ I3 < 0 ] <- 0

I4[ 140, 140 ] <- 1
I4 <- kernel2dsmooth( I4, kernel.type = "disk", r = 20 )
I4[ I4 > 0 ] <- 1
if( any( I4 < 0 ) ) I4[ I4 < 0 ] <- 0

image( I1, col = c("white", "darkblue") )
contour( I2, add = TRUE )
contour( I3, add = TRUE )
contour( I4, add = TRUE )

# Each circle is the same size and shape, and the domain is square.
# I1 and I2, I2 and I3, and I2 and I4 are all the same distance
# away from each other.  I1 and I4 and I3 and I4 are also the same distance
# from each other.  I3 touches the edge of the domain.
# 

# First, calculate the Baddeley delta metric on each 
# comparison.

I1im <- as.im( I1 )
I2im <- as.im( I2 )
I3im <- as.im( I3 )
I4im <- as.im( I4 )

I1im <- solutionset( I1im > 0 )
I2im <- solutionset( I2im > 0 )
I3im <- solutionset( I3im > 0 )
I4im <- solutionset( I4im > 0 )

deltametric( I1im, I2im )
deltametric( I2im, I3im )
deltametric( I2im, I4im )

# Above are all different values.
# Below, they are all 28.84478.
censqdelta( I1, I2 )
censqdelta( I2, I3 )
censqdelta( I2, I4 )

# Similarly for I1 and I4 vs I3 and I4.
deltametric( I1im, I4im )
deltametric( I3im, I4im )

censqdelta( I1, I4 )
censqdelta( I3, I4 )

# To see why this problem exists.
dm1 <- distmap( I1im )
dm1 <- as.matrix( dm1 )
dm2 <- distmap( I2im )
dm2 <- as.matrix( dm2 )

par( mfrow = c( 2, 2 ) )
image.plot( dm1 )
contour( I1, add = TRUE, col = "white" )
image.plot( dm2 )
contour( I2, add = TRUE, col = "white" )

image.plot( abs( dm1 ) - abs( dm2 ) )
contour( I1, add = TRUE, col = "white" )
contour( I2, add = TRUE, col = "white" )

}

}
\keyword{ arith }% use one of  RShowDoc("KEYWORDS")
\keyword{ manip }% __ONLY ONE__ keyword per line
