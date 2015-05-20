\name{Fint2d}
\alias{Fint2d}
\title{
2-d Interpolation
}
\description{
Interpolate a function of two variables by rounding (i.e. taking the nearest value), bilinear or bicubic interpolation.
}
\usage{
Fint2d(X, Ws, s, method = c("round", "bilinear", "bicubic"), derivs = FALSE, ...)
}
\arguments{
  \item{X}{
A numeric n by m matrix giving the value of the function at the old coordinates.
}
  \item{Ws}{
A numeric k by 2 matrix of new grid coordinates where k <= m * n.
}
  \item{s}{
A numeric k by 2 matrix of old grid coordinates where k <= m * n.
}
  \item{method}{
character naming one of \dQuote{round} (default), \dQuote{bilinear}, or \dQuote{bicubic} giving the specific interpolation method to use.
}
  \item{derivs}{
logical, should the gradient interpolatants be returned?
}
  \item{\dots}{
Not used.
}
}
\details{
Method round simply returns the values at each grid point that correspond to the nearest points in the old grid.

Interpolation of a function, say H, is achieved by the following formula (cf. Gilleland et al 2010, sec. 3), where r and s represent the fractional part of their respective coordinate.  that is, r = x - g( x ) and s = y - g( y ), where g( x ) is the greatest integer less than x. 

    sum_k sum_l b_k( r ) * b_l( s ) * H(g( x ) + l,  g( y ) + k).

The specific choices for the values of b_l and b_k and their ranges depends on the type of interpolation.  For bilinear interpolation, they both range from 0 to 1, and are given by: b_0( x ) = 1 - x and b_1( x ) = x.  for bicubic interpolation, they both range from -1 to 2 and are given by:

    b_{-1}( t ) = (2 * t^2 - t^3 - t) / 2

    b_0( t ) = (3 * t^3 - 5 * t^2 + 2) / 2

    b_1( t ) = (4 * t^2 - 3 * t^3 + t) / 2

    b_2( t ) = ((t - 1) * t^2) / 2.

}
\value{
If deriv is FALSE, then a matrix is returned whose values correspond to the new coordinates.  Otherwise a list is returned with components:
  \item{xy }{matrix whose values correspond to the new coordinates.}
  \item{dx, dy }{ matrices giving the x and y direction gradients of the interpolation.}
}
\references{
Gilleland and co-authors (2010) Spatial forecast verification: Image warping.  \emph{NCAR Technical Note}, NCAR/TN-482+STR, DOI: 10.5065/D62805JJ.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{rigider}}, \code{\link{rigidTransform}}
}
\examples{

# see rigider for an example.

}
\keyword{ manip }
\keyword{ arith }
\keyword{ math }
