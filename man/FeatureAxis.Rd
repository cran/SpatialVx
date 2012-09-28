\name{FeatureAxis}
\alias{FeatureAxis}
\alias{FeatureProps}
\alias{plot.FeatureAxis}
\title{
Calculate the major (and minor) axis of a feature
}
\description{
Calculate the major and minor axes of a feature and various other properties such as the aspect ratio.
}
\usage{
FeatureAxis(x, fac = 1, flipit = FALSE, twixt = FALSE)
\method{plot}{FeatureAxis}(x, ...)
FeatureProps(x, Im = NULL, which.props=c("centroid", "area", "axis", "intensity"),
	areafac=1, q = c(0.25, 0.9), ...)
}
\arguments{
  \item{x}{
For \code{FeatureAxis} this is an object of class "owin" containing a binary image matrix defining the feature.  In the case of \code{plot.FeatureAxis}, this is the value returned from \code{FeatureAxis}.
}
   \item{Im}{Matrix giving the original values of the field from which the feature was extracted.  Only needed if the feature intensity is desired.}
  \item{fac,areafac}{
numeric, in determining the lengths of the axes, they are multiplied by a factor of \code{fac} (e.g., if the grid points are k X k km each, then one could set this to k so that the resulting lengths are in terms of km rather than grid points.
}
   \item{q}{numeric vector of values between 0 and 1 inclusive giving the quantiles for determining the intensity of the feature.}
  \item{flipit}{logical, should the objects be flipped over x and y?  The disjointer function results in images that are flipped, this would flip them back.
}
  \item{twixt}{
logical, should the major axis angle be forced to be between +/- 90 degrees?
}
   \item{which.props}{character vector giving one or more of "centroid", "area", "axis" and "intensity".  If "axis" is given, then a call to \code{FeatureAxis} is made.}
   \item{\dots}{For \code{plot.FeatureAxis} these are additional arguments to the \code{plot} function.  For \code{FeatureProps}, these are additional arguments to \code{FeatureAxis}.}
}
\details{
This function attempts to identify the major and minor axes for a pre-defined feature (sometimes referred to as an object).  This function relies heavily on the spatstat and smatr packages.  First, the convex hull of the feature is determined using the convexhull function from the spatstat package.  The major axis is then found using the sma funciton from package smatr, which is then converted into a psp object (see as.psp from spatstat) from which the axis angle and length are found (using angles.psp and lengths.psp, resp., from spatstat).

The minor axis anlge is easily found after rotating the major axis 90 degrees using rotate.psp from spatstat.  The length of the minor axis is more difficult.  Here, it is found by rotating the convex hull of the feature by the major axis angle (so that it is upright) using rotate.owin from spatstat, and then computing the bounding box (using bounding.box from spatstat).  An infinite line (using infline from spatstat) is created, and then clipped (using clip.infline from spatstat) by the upright bounding box.  This seems to give a reasonable value for the length of the minor axis.  A psp object is then created using the mid point of the major axis (which should be close to the centroid of the feature) using as.psp and midpoints.psp from spatstat along with the length and angle already found for the minor axis.

See the help files for the above mentioned functoins for references, etc.
}
\value{
FeatureAxis: A list object of class "FeatureAxis" is returned with components:
   \item{z }{same as the argument x passed in.}
   \item{MajorAxis,MinorAxis }{a psp object with one segment that is the major (minor) axis.}
   \item{OrientationAngle}{list with two components: MajorAxis (the angle in degrees of the major axis wrt the abscissa), MinorAxis (the angle in degrees wrt the abscissa).}
   \item{aspect.ratio}{numeric giving the ratio of the length of the minor axis to that of the major axis (always between 0 and 1).}
   \item{MidPoint}{an object of class "ppp" giving the mid point of the major (minor) axis.}
   \item{lengths}{list object with components: MajorAxis giving the length (possibly multiplied by a factor) of the major axis, and MinorAxis same as MajorAxis but for the minor axis.}
   \item{sma.fit}{The fitted object returned by the sma function.  This is useful, e.g., if confidence intervals for the axis are desired.  See the sma help file for more details.}

No value is returned from plot.FeatureAxis, a plot is created.

For FeatureProps, a list object with components depending on the which.props argument.  One or more of:
   \item{centroid}{list with components x and y giving the centroid of the object.}
   \item{area}{numeric giving the area of the feature.}
   \item{axis}{list object of class FeatureAxis as returned by the same-named function.}
}
\author{
Eric Gilleland
}

\section{Warning }{The length of the minor axis is determined by the clipping of an infinite line with a bounding box.  Therefore, it will be slightly longer than it should be.  In practice, this does not seem to be horrible.  Note, however, that for circular shaped features, this may lead to a minor axis that is apparently longer than the major axis (i.e., aspect ratio will be greater than 1).  Perhaps this will be improved in the future.}

\seealso{
\code{\link{owin}}, \code{\link{convexhull}}, \code{\link{sma}}, \code{\link{as.psp}}, \code{\link{angles.psp}}, \code{\link{rotate.owin}}, \code{\link{rotate.psp}}, \code{\link{bounding.box}}, \code{\link{midpoints.psp}}, \code{\link{lengths.psp}}, \code{\link{infline}}, \code{\link{clip.infline}}, \code{\link{deltamm}}, \code{\link{FeatureSuite}}, \code{\link{FeatureSuitePrep}}, \code{\link{convthresh}}, \code{\link{disjointer}}, \code{\link{connected}}, \code{\link{tiles}}, \code{\link{tess}}, \code{\link{solutionset}}
}
\examples{
grid<- list( x= seq( 0,5,,100), y= seq(0,5,,100)) 
obj <- stationary.image.cov( grid=grid, V=cbind(c(0.25,0.15),c(0.45,0.75)), setup=TRUE)
look<- sim.rf( obj)
look <- disk2dsmooth(look,5)
u <- quantile(look,0.99)
sIx <- matrix(0, 100, 100)
sIx[ look > u] <- 1
look2 <- disjointer(sIx)[[1]]
look2 <- flipxy(look2)
tmp <- FeatureAxis(look2)
plot(tmp)
tmp$OrientationAngle
tmp$lengths

\dontrun{
data(pert000)
data(pert004)
hold <- FeatureSuitePrep("pert004", "pert000")
look <- convthresh( hold, smoothpar=10.5)
par( mfrow=c(2,2))
image.plot(look$X.labeled)
image.plot(look$Y.labeled)

look2 <- deltamm( look, verbose=TRUE)
image.plot(look2$X.labeled)
image.plot(look2$Y.labeled)

look2$mm.new.labels # the first seven features are matched.

ang1 <- FeatureAxis(look2$X.feats[[1]])
ang2 <- FeatureAxis(look2$Y.feats[[1]])
plot(ang1)
plot(ang2)

ang3 <- FeatureAxis(look2$X.feats[[4]])
ang4 <- FeatureAxis(look2$Y.feats[[4]])
plot(ang3)
plot(ang4)
   }

}
\keyword{ math }
\keyword{ graphs }
