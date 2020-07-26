\name{saller}
\alias{saller}
\alias{print.saller}
\alias{summary.saller}
\title{
Feature-based Analysis of a Field (Image)
}
\description{
Feature-based analysis of a field (image) 
}
\usage{
saller(x, d = NULL, distfun = "rdist", ...)

\method{print}{saller}(x, ...)

\method{summary}{saller}(object, ...)
}
\arguments{
  \item{x}{
\code{saller}: \code{x} is a list object returned by \code{FeatureFinder} or other feature identification function that returns a list with components \code{X.feats}, \code{Y.feats} (themselves lists with owin class objects defining separate features in the verification and forecast fields, resp.), and \code{X.labeled}, \code{Y.labeled} (fields with the numbers from 0 to the number of features also defining the separate feature locations (e.g., as returned by the \code{connected} function of package \pkg{spatstat}.

\code{print}: list object returned by \code{saller}.
}
  \item{object}{

\code{summary}: object the returned by \code{saller}.
}
  \item{d}{
(optional) the SAL (\code{saller}) method requires division by the longest distance between two border points.  If NULL, this is taken to be simply the length of the longest side.
}
  \item{distfun}{Function with which to calculate centroid distances.  Default uses straight Euclidean.  To do great-circle distance, use \code{rdist.earth} and be sure that \code{object} has a loc attribute with lon/lat coordinates.}
\item{\dots}{Optional arguments to \code{distfun}.  Not used by \code{print} or \code{summary}.}
}
\details{
saller: Computes S, A, and L of the SAL method introduced by Wernli et al. (2008).
}
\value{
saller returns a list with components:
  \item{A }{numeric giving the amplitude component.}
  \item{L }{numeric giving the lcoation component.}
  \item{S}{numeric giving the structure component.}
  \item{L1,L2}{numeric giving the values that sum together to give L.}
  \item{L1.alt, L.alt}{numeric giving an alternative L1 component, and subsequently alternative L where it is calculated using the centroid of the field containing only defined features rather than the original raw field.}

print invisibly returns a named vector with S, A and L.

summary does not return anything.

}
\references{
Wernli, H., Paulat, M., Hagen, M. and Frei, C. (2008) SAL--A novel quality measure for the verification of quantitative precipitation forecasts.  \emph{Mon. Wea. Rev.}, \bold{136}, 4470--4487, doi:10.1175/2008MWR2415.1.
}
\author{
Eric Gilleland
}
\note{
There are several ways to identify features, and some are provided by this package, but only a few.  For example, the method for identifying features in the SAL method as introduced by Wernli et al. (2008) utilizes information from a contour field of a particular variable, and is therefore not currently included in this package.  Users are encouraged to write their own such functions, and should feel free to contribute them to this package by contacting the maintainer.

The SAL method typically looks at a small domain, and it is up to the user to set this up before calling these functions, as they are not designed to handle such a situation.

}

\seealso{
\code{\link{centroid.owin}}, \code{\link{connected}}, \code{\link{tiles}}, \code{\link{tess}}, \code{\link{deltamm}}, \code{\link{make.SpatialVx}}
}
\examples{

data( "ExampleSpatialVxSet" )

x <- ExampleSpatialVxSet$vx
xhat <- ExampleSpatialVxSet$fcst

q <- mean( c(c(x[x>0]),c(xhat[xhat>0])), na.rm=TRUE)

hold <- make.SpatialVx( x, xhat, field.type="contrived", units="none",
    data.name = "Example", obs.name = "x", model.name = "xhat" )

hold2 <- FeatureFinder(hold, smoothpar=5, thresh=q)
\dontrun{plot(hold2)}

look <- saller(hold2)   
summary(look)


}
\keyword{ math }
