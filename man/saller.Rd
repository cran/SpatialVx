\name{saller}
\alias{saller}
\alias{print.saller}
\alias{summary.saller}
\title{
Analyze Features in a Field (Image)
}
\description{
Functions to analyze features in a field (image) in conjunction, e.g., with the FeatureSuite function.
}
\usage{
saller(object, x, y = NULL, matches = NULL, d = NULL, time.point = 1, 
    model = 1)

\method{print}{saller}(x, ...)

\method{summary}{saller}(object, ...)
}
\arguments{
  \item{x}{
\code{saller}: \code{x} is a list object returned by \code{convthresh} or other feature identification function that returns a list with components \code{X.feats}, \code{Y.feats} (themselves lists with owin class objects defining separate features in the verification and forecast fields, resp.), and \code{X.labeled}, \code{Y.labeled} (fields with the numbers from 0 to the number of features also defining the separate feature locations (e.g., as returned by the \code{connected} function of package \pkg{spatstat}.

\code{print}: list object returned by \code{saller}.
}
  \item{object}{
\code{saller}: A list object of class \dQuote{SpatialVx} as returned by \code{make.SpatialVx} providing information about the original fields.

\code{summary}: object the returned by \code{saller}.
}
  \item{y}{
(optional) if not NULL, it is used instead of \code{x}, and otherwise has the same stucture thereof.  When used with \code{FeatureSuite}, one might want to merge objects in one or both fields before calling the analysis function.  If so, \code{y} is the means by which \code{FeatureSuite} can ensure that this happens.
}
  \item{matches}{
Not used by \code{saller} function.  If matching of features between fields is required, this will be a two-column matrix describing which features from the forecast (first column) match to those of the verificaiton field (second column).
}
  \item{d}{
(optional) the SAL (\code{saller}) method requires division by the longest distance between two border points.  If NULL, this is taken to be simply the length of the longest side.
}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{model}{numeric indicating which forecast model to select for the analysis.}
\item{\dots}{not used.}
}
\details{
saller: Computes S, A, and L of the SAL method introduced by Wernli et al. (2008).
}
\value{
saller returns a list with components:
  \item{A }{numeric giving the amplitude component.}
  \item{L }{numeric giving the lcoation component.}
  \item{S}{numeric giving the structure component.}

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
\code{\link{FeatureSuite}}, \code{\link{centroid.owin}}, \code{\link{connected}}, \code{\link{tiles}}, \code{\link{tess}}, \code{\link{deltamm}}, \code{\link{make.SpatialVx}}
}
\examples{

data(ExampleSpatialVxSet)

x <- ExampleSpatialVxSet$vx
xhat <- ExampleSpatialVxSet$fcst

q <- mean( c(c(x[x>0]),c(xhat[xhat>0])), na.rm=TRUE)

hold <- make.SpatialVx(x, xhat, field.type="contrived", units="none",
    data.name=c("Example", "x", "xhat"))

hold2 <- convthresh(hold, smoothpar=5, thresh=q)
\dontrun{plot(hold2, set.pw=TRUE)}

look <- saller(hold, hold2)   
summary(look)

\dontrun{
data(pert000)
data(pert004)
data(ICPg240Locs)

hold <- make.SpatialVx(pert000, pert004, loc=ICPg240Locs,
    projection=TRUE, map=TRUE,
    field.type="Precipitation", units="mm/h",
    data.name=c("Perturbed ICP Cases", "pert000", "pert004"))

look <- convthresh(hold, smoothpar=10.5)
summary(look)
plot(look, set.pw=TRUE)

saller(hold, look)

   }

}
\keyword{ math }
