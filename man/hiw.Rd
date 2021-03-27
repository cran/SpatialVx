\name{hiw}
\alias{hiw}
\alias{distill.hiw}
\alias{plot.hiw}
\alias{print.hiw}
\alias{summary.hiw}

\title{
Spatial Forecast Verification Shape Analysis
}

\description{
Shape analysis for spatial forecast verification (hiw is OE for shape; yields MnE hue).
}

\usage{
hiw(x, simplify = 0, A = pi * c(0, 1/16, 1/8, 1/6, 1/4, 1/2, 9/16, 5/8, 2/3, 3/4),
    verbose = FALSE, ...)

\method{distill}{hiw}(x, ...)

\method{plot}{hiw}(x, ..., which = c("X", "Xhat"), ftr.num = 1, zoom = TRUE, 
    seg.col = "darkblue")

\method{print}{hiw}(x, ...)

\method{summary}{hiw}(object, ..., silent = FALSE)

}
\arguments{
  \item{x,object}{
\code{hiw}: object of class \dQuote{features}.

\code{distill}, \code{plot}, \code{summary}: object of class \dQuote{hiw}.
}
  \item{simplify}{
\code{dmin} argument in call to \code{simplify.owin} from \pkg{spatstat}.  If 0 (default), then no call is made to \code{simplify.owin}.
}
  \item{A}{
numeric vector of angles for which to apply shape analysis.  Note that this vector will be rounded to 6 digits.  If values are less than that, might be prudent to add 1e-6 to them.
}
  \item{verbose}{logical, should progress information be printed to the screen?}
  \item{which}{character string naming whether to plot a feature from the obsevation field (default) or the forecast field.}
  \item{ftr.num}{integer stating which feature number to plot.}
  \item{zoom}{logical, should the feature be plotted within its original domain, or a blow-up of the feature (default)?}
  \item{seg.col}{color for the line segments.}
  \item{silent}{logical, should the summary information be printed to the screen?}
  \item{\dots}{
Not used by \code{hiw}, \code{distill}, \code{summary}.  Optional arguments to \code{plot}.
}
}
\details{
This function is an attempt to approximate the technique described first in Micheas et al. (2007) and as modified in Lack et al. (2010).  It will only find the centroids, rays extending from them to the boundaries, and the boundary points.  Use \code{distill} to convert this output into an object readable by, for example, \code{procGPA} from package \pkg{shapes}.

First, identified features (which may be identified by any feature identification function that yields an object of class \dQuote{features}) are taken, the centroid is found (the centroid is found via \code{centroid.owin} so that the x- and y- coordinates are fliped from what you might expect) and very long line segments are found radiating out in both directions from the center.  They are then clipped by where they cross the boundaries of the features.

The \pkg{spatstat} package is used heavily by this function.  In particular, the function \code{as.polygonal} is applied to the \code{owin} objects (possible after first calling \code{simplify.owin}).  Line segments are created using the feature centroids, as found by \code{centroid.owin}, and the user-supplied angles, along with a very long length (equal to the domain size).  Boundary crossings are found using \code{crossing.psp}, and new line segment patterns are created using the centroids and boundary crossing information (extra points along line segments are subsequently removed through a painstaking process, and \code{as.psp} is called again, and any missing line segments are subsequently accounted for, for later calculations).  Additionally, lengths of line segments are found via \code{lengths_psp}.  Angles must also be re-determined and corresponded to the originally passed angles.  Therefore, it is necessary to round the angles to 6 digits, or \dQuote{equal} angles may not be considered equal, which will cause problems.

The \code{hiw} function merely does the above step, as well as finds the lengths of the resulting line segments.  For non-convex objects, the longest line segment is returned, and if the boundary crossings do not lie on opposite sides of the centroid, then the negative of the shortest segment is returned for that particular value.  Also returned are the mean, min and maximum intensities for each feature, as well as the final angles returned.  It is possible to have missing values for some of these components.

The \code{summary} function computes SSloc, SSavg, SSmin, and SSmax between each pair of features between fields.  \code{distill} may be used to create an object that can be further analyzed (for shape) using the \pkg{shapes} package.

While any feature identification function may be used, it is recommended to throw out small sized features as the results may be misleading (e.g., comparisons between features consisting of single points, etc.).
}
\value{
A list object of class \dQuote{hiw} is returned with components the same as in the original \dQuote{features} class object, as well as:
  \item{radial.segments }{a list with components X and Xhat each giving lists of the \dQuote{psp} class (i.e., line segment) object for each feature containing the radial segments from the feature centroids to the boundaries.}
  \item{centers}{list with components X and Xhat giving two-column matrices containing the x- and y- coordinate centroids for each feature (as determined by centroid.owin).}
  \item{intensities}{list with components X and Xhat giving three-column matrices that contain the mean, min and max intensities for each feature.}
  \item{angles,lengths}{list with components X and Xhat each giving lists containing the lengths of the line segments and their respective angles.  Missing values are possible here.}

distill returns an array whose dimensions are the number of landmarks (i.e., boundary points) by two by the number of observed and forecast features.  An attribute called \dQuote{field.identifier} is also given that is a character vector containing repeated \dQuote{X} and \dQuote{Xhat} values identiying which of the third dimension are associated with the observed field (X) and those identified with the forecast field (Xhat).  Note that missing values may be present, which may need to be dealt with (by the user) before using functions from the shapes package.

summary invisibly returns a list object with components:

  \item{X,Xhat}{matrices whose rows represent features and whose columns give their centroids (note that x refers to the columns and y to the rows), as well as the average, min and max intensities.}
  \item{SS}{matrix with four rows and columns equal to the number of possible combinations of feature matchings between fields.  Gives the sum of square translation/location error (i.e., squared centroid distance), as well as the average, min and max squared differences between each combination of features.}
  \item{ind}{two-column matrix whose rows indicate the feature numbers from each field being compared; corresponding to the columns of SS above.}
}
\references{
Lack, S. A., Limpert, G. L. and Fox, N. I. (2010) An object-oriented multiscale verification scheme.  \emph{Wea. Forecasting}, \bold{25}, 79--92.

Micheas, A. C., Fox, N. I., Lack, S. A., and Wikle, C. K. (2007) Cell identification and verification of QPF ensembles using shape analysis techniques.  \emph{J. Hydrology}, \bold{343}, 105--116.
}
\author{
Eric Gilleland
}

\seealso{

To indentify features and create objects of class \dQuote{features}, see, for example: \code{\link{FeatureFinder}}

\code{\link{centroid.owin}}, \code{\link{as.psp}}, \code{\link{psp}}, \code{\link{crossing.psp}}, \code{\link{lengths_psp}}, \code{\link{angles.psp}}
}
\examples{

data( "geom000" )
data( "geom001" )
data( "geom004" )
data( "ICPg240Locs" )

hold <- make.SpatialVx( geom000, geom001, thresholds = c(0.01, 50.01),
    projection = TRUE, map = TRUE, loc = ICPg240Locs, loc.byrow = TRUE,
    field.type = "Geometric Objects Pretending to be Precipitation",
    units = "mm/h", data.name = "ICP Geometric Cases", obs.name = "geom000",
    model.name = "geom001" )

look <- FeatureFinder(hold, do.smooth = FALSE, thresh = 2, min.size = 200)

look <- hiw(look)

distill.hiw(look)

# Actually, you just need to type:
# distill(look)

summary(look)

# Note: procGPA will not allow missing values.

par(mfrow=c(1,2))
plot(look)
plot(look, which = "Xhat")

}
\keyword{ misc }
