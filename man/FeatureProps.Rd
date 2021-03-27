\name{FeatureProps}
\alias{FeatureProps}
\title{
Single Feature Properties
}
\description{
Calculate properties for an identified feature.
}
\usage{
FeatureProps(x, Im = NULL, which.props = c("centroid", "area", "axis", "intensity"),
    areafac = 1, q = c(0.25, 0.9), loc = NULL, ...)
}
\arguments{
  \item{x}{
object of class \dQuote{owin} containing a binary image matrix defining the feature.
}
  \item{Im}{
Matrix giving the original values of the field from which the feature was extracted.  Only needed if the feature intensity is desired.
}
  \item{which.props}{
character vector giving one or more of \dQuote{centroid}, \dQuote{area}, \dQuote{axis} and \dQuote{intensity}.  If \dQuote{axis} is given, then a call to \code{FeatureAxis} is made.
}
  \item{areafac}{
numeric, in determining the lengths of the axes, they are multiplied by a factor of \code{fac} (e.g., if the grid points are k by k km each, then one could set this to k so that the resulting lengths are in terms of km rather than grid points.
}
  \item{q}{
numeric vector of values between 0 and 1 inclusive giving the quantiles for determining the intensity of the feature.
}
   \item{loc}{optional argument giving a two-column matrix of grid locations for finding the centroid.  If NULL, indices based on the dimension of x are used.}
  \item{\dots}{
additional arguments to \code{FeatureAxis}.
}
}
\details{
This function takes an \code{owin} image and returns several property values for that image, including: centroid, spatial area, major and minor axis angle/length, as well as the overall intensity of the field (cf., Davis et al., 2006a, b).
}
\value{
list object with components depending on the which.props argument.  One or more of:
   \item{centroid}{list with components x and y giving the centroid of the object.}
   \item{area}{numeric giving the area of the feature.}
   \item{axis}{list object of class FeatureAxis as returned by the same-named function.}
}
\references{
Davis, C. A., Brown, B. G. and Bullock, R. G. (2006a) Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.

Davis, C. A., Brown, B. G. and Bullock, R. G. (2006b) Object-based verification of precipitation forecasts, Part II: Application to convective rain systems. \emph{Mon. Wea. Rev.}, \bold{134}, 1785--1795.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{FeatureAxis}}, \code{\link{owin}}, \code{\link{convexhull}}, \code{\link{sma}}, \code{\link{as.psp}}, \code{\link{angles.psp}}, \code{\link{rotate.owin}}, \code{\link{rotate.psp}}, \code{\link{boundingbox}}, \code{\link{midpoints.psp}}, \code{\link{lengths_psp}}, \code{\link{infline}}, \code{\link{clip.infline}}, \code{\link{deltamm}}, \code{\link{FeatureFinder}}, \code{\link{disjointer}}, \code{\link{connected}}, \code{\link{tiles}}, \code{\link{tess}}, \code{\link{solutionset}}
}
\examples{
data( "ExampleSpatialVxSet" )

x <- ExampleSpatialVxSet$vx

look <- disk2dsmooth(x,5)
u <- quantile(look,0.99)
sIx <- matrix(0, 100, 100)
sIx[ look > u] <- 1
look2 <- disjointer(sIx)[[1]]
look2 <- flipxy(look2)

FeatureProps(look2,
    which.props=c("centroid", "area", "axis"))

}
\keyword{ math }
