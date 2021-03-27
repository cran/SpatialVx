\name{FeatureAxis}
\alias{FeatureAxis}
\alias{plot.FeatureAxis}
\alias{summary.FeatureAxis}
\title{
Major and Minor Axes of a Feature
}
\description{
Calculate the major and minor axes of a feature and various other properties such as the aspect ratio.
}
\usage{
FeatureAxis(x, fac = 1, flipit = FALSE, twixt = FALSE)

\method{plot}{FeatureAxis}(x, ..., zoom = FALSE)

\method{summary}{FeatureAxis}(object, ...)

}
\arguments{
  \item{x}{
For \code{FeatureAxis} this is an object of class \dQuote{owin} containing a binary image matrix defining the feature.  In the case of \code{plot.FeatureAxis}, this is the value returned from \code{FeatureAxis}.
}

  \item{object}{list object of class \dQuote{FeatureAxis} as returned by \code{FeatureAxis}.}
  \item{fac}{
numeric, in determining the lengths of the axes, they are multiplied by a factor of \code{fac} (e.g., if the grid points are k by k km each, then one could set this to k so that the resulting lengths are in terms of km rather than grid points.
}
  \item{flipit}{logical, should the objects be flipped over x and y?  The disjointer function results in images that are flipped, this would flip them back.
}
  \item{twixt}{
logical, should the major axis angle be forced to be between +/- 90 degrees?
}
   \item{zoom}{logical, should the object be plotted on its bounding box (TRUE) or on the original grid (FALSE, default)?  Useful if the feature is too small to be seen well on the original gid.}
   \item{\dots}{For \code{plot.FeatureAxis} these are additional arguments to the \code{plot} function.  Not used by \code{summary.FeatureAxis}.}
}
\details{
This function attempts to identify the major and minor axes for a pre-defined feature (sometimes referred to as an object).  This function relies heavily on the \pkg{spatstat} and \pkg{smatr} packages.  First, the convex hull of the feature is determined using the \code{convexhull} function from the \pkg{spatstat} package.  The major axis is then found using the \code{sma} function from package \pkg{smatr}, which is then converted into a \code{psp} object (see \code{as.psp} from \pkg{spatstat}) from which the axis angle and length are found (using \code{angles.psp} and \code{lengths_psp}, resp., from \pkg{spatstat}).

The minor axis anlge is easily found after rotating the major axis 90 degrees using \code{rotate.psp} from \pkg{spatstat}.  The length of the minor axis is more difficult.  Here, it is found by rotating the convex hull of the feature by the major axis angle (so that it is upright) using \code{rotate.owin} from \pkg{spatstat}, and then computing the bounding box (using \code{boundingbox} from \pkg{spatstat}).  The differnce is then taken between the range of x- coordinates of the bounding box.  This seems to give a reasonable value for the length of the minor axis.  A \code{psp} object is then created using the mid point of the major axis (which should be close to the centroid of the feature) using \code{as.psp} and \code{midpoints.psp} from \pkg{spatstat} along with the length and angle already found for the minor axis.

See the help files for the above mentioned functions for references, etc.
}
\value{
FeatureAxis: A list object of class \dQuote{FeatureAxis} is returned with components:
   \item{z }{same as the argument x passed in.}
   \item{MajorAxis,MinorAxis }{a psp object with one segment that is the major (minor) axis.}
   \item{OrientationAngle}{list with two components: MajorAxis (the angle in degrees of the major axis wrt the abscissa), MinorAxis (the angle in degrees wrt the abscissa).}
   \item{aspect.ratio}{numeric giving the ratio of the length of the minor axis to that of the major axis (always between 0 and 1).}
   \item{MidPoint}{an object of class \dQuote{ppp} giving the mid point of the major (minor) axis.}
   \item{lengths}{list object with components: MajorAxis giving the length (possibly multiplied by a factor) of the major axis, and MinorAxis same as MajorAxis but for the minor axis.}
   \item{sma.fit}{The fitted object returned by the sma function.  This is useful, e.g., if confidence intervals for the axis are desired.  See the sma help file for more details.}

No value is returned from the \code{plot} or \code{summary} method functions.

}
\author{
Eric Gilleland
}

\seealso{
\code{\link{owin}}, \code{\link{convexhull}}, \code{\link{sma}}, \code{\link{as.psp}}, \code{\link{angles.psp}}, \code{\link{rotate.owin}}, \code{\link{rotate.psp}}, \code{\link{boundingbox}}, \code{\link{midpoints.psp}}, \code{\link{lengths_psp}}, \code{\link{infline}}, \code{\link{clip.infline}}, \code{\link{deltamm}}, \code{\link{FeatureFinder}}, \code{\link{disjointer}}, \code{\link{connected}}, \code{\link{tiles}}, \code{\link{tess}}, \code{\link{solutionset}}
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
tmp <- FeatureAxis(look2)
plot(tmp)
summary(tmp)

\dontrun{
data( "pert000" )
data( "pert004" )
data( "ICPg240Locs" )

hold <- make.SpatialVx( pert000, pert004,
    loc = ICPg240Locs, projection = TRUE, map = TRUE,
    loc.byrow = TRUE,
    field.type = "Precipitation", units = "mm/h",
    data.name = "Perturbed ICP Cases", obs.name = "pert000",
    model.name = "pert004" )

look <- FeatureFinder(hold, smoothpar=10.5)
par(mfrow=c(1,2))
plot(look)

par(mfrow=c(2,2))
image.plot(look$X.labeled)
image.plot(look$Y.labeled)

# The next line will likely be very slow.
look2 <- deltamm(x=look, verbose=TRUE)
image.plot(look2$X.labeled)
image.plot(look2$Y.labeled)

look2$mm.new.labels # the first seven features are matched.

ang1 <- FeatureAxis(look2$X.feats[[1]])
ang2 <- FeatureAxis(look2$Y.feats[[1]])
plot(ang1)
plot(ang2)
summary(ang1)
summary(ang2)

ang3 <- FeatureAxis(look2$X.feats[[4]])
ang4 <- FeatureAxis(look2$Y.feats[[4]])
plot(ang3)
plot(ang4)
summary(ang3)
summary(ang4)
   }

}
\keyword{ math }
\keyword{ graphs }
