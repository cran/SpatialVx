\name{imomenter}
\alias{imomenter}
\alias{imomenter.im}
\alias{imomenter.matrix}
\alias{print.imomented}
\title{
Image Moments
}
\description{
Calculate some of the raw image moments, as well as some useful image characteristics.
}
\usage{
imomenter(x, loc = NULL, ...)

\method{imomenter}{im}(x, loc = NULL, ...)

\method{imomenter}{matrix}(x, loc = NULL, ...)

\method{print}{imomented}(x, ...)
}
\arguments{
  \item{x}{
\code{imomenter}: matrix or object of class \dQuote{im} (from package \pkg{spatstat}).

\code{print}: object of class \dQuote{imomenter}.
}
  \item{loc}{
A two-column matrix giving the location coordinates.  May be missing in which case they are assumed to be integers giving the row and column numbers.
}
  \item{\dots}{
Not used.
}
}
\details{
Calculates Hu's image moments (Hu 1962).  Calculates the raw moments: M00 (aka area), M10, M01, M11, M20, and M02, as well as the (normalized) central moments: mu11', mu20', and mu02', which are returned as the image covariance matrix: rbind(c(mu20', mu11'), c(mu11', mu02')).  In addition, the image centroid and orientation angle are returned, as calculated using the image moments.  It should be noted that while the centroid is technically defined for the null case (all zero-valued grid points), the way it is calculated using image moments means that it will be undefined because of division by zero in the formulation.

The orientation angle calculated here is that which is used by MODE, although not currently used in the MODE analyses in this package (smatr is used instead to find the major axis, etc).  The eigenvalues of the image covariance correspond to the major and minor axes of the image.

For more information on image moments, see http://en.wikipedia.org/wiki/Image_moments and references therein.
}
\value{
A list object of class \dQuote{imomented} is returned with components:
  \item{area}{Same as M00.}
  \item{centroid}{numeric with named components \dQuote{x} and \dQuote{y} giving the x- and y- coordinates of the centroid as calculated by the image moment method.}
  \item{orientation.angle}{The orientation angle of the image as calculated by image moments.}
  \item{raw.moments}{named numeric vector with the raw image moments: M00, M10, M01, M11, M20 and M02 used in calculating the other returned values.}
  \item{cov}{2 by 2 image covariance as calculated by the image moment method.}
}
\references{
Hu, M. K. (1962) Visual Pattern Recognition by Moment Invariants. \emph{IRE Trans. Info. Theory}, \bold{IT-8}, 179--187.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{Mij}}, \code{\link{FeatureAxis}}
}
\examples{

look <- matrix(0, 10, 10)

look[3:5, 7:8] <- rnorm(6)

imomenter(look)

\dontrun{
data( "geom000" )
data( "ICPg240Locs" )

imomenter( geom000 )
imomenter( geom000, loc = ICPg240Locs )

data( "geom004" )

imomenter( geom004 )

imomenter( geom004, loc = ICPg240Locs )
}

}
\keyword{manip}
\keyword{arith}
