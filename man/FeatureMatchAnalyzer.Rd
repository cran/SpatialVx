\name{FeatureMatchAnalyzer}
\alias{FeatureMatchAnalyzer}
\alias{FeatureMatchAnalyzer.matched.centmatch}
\alias{FeatureMatchAnalyzer.matched.deltamm}
\alias{plot.FeatureMatchAnalyzer}
\alias{print.FeatureMatchAnalyzer}
\alias{summary.FeatureMatchAnalyzer}
\alias{FeatureComps}
\alias{distill.FeatureComps}
\title{
Analyze Features of a Verification Set
}
\description{
Analyze matched features of a verification set.
}
\usage{
FeatureMatchAnalyzer(x, which.comps=c("cent.dist", "angle.diff", "area.ratio", "int.area",
                    "bdelta", "haus", "ph", "med", "msd", "fom", "minsep",
		    "bearing"), sizefac=1, alpha=0.1, k=4, p=2, c=Inf,
		    distfun="distmapfun", ...)

\method{FeatureMatchAnalyzer}{matched.centmatch}(x, which.comps=c("cent.dist", "angle.diff",
		    "area.ratio", "int.area", "bdelta", "haus", "ph", "med",
		    "msd", "fom", "minsep", "bearing"), sizefac=1, alpha=0.1, k=4, p=2,
		    c=Inf, distfun="distmapfun", ...)

\method{FeatureMatchAnalyzer}{matched.deltamm}(x, which.comps = c("cent.dist", "angle.diff",
		    "area.ratio", "int.area", "bdelta", "haus", "ph", "med", "msd",
		    "fom", "minsep", "bearing"), sizefac = 1, alpha = 0.1, k = 4, p = 2,
		    c = Inf, distfun = "distmapfun", ..., y = NULL, matches = NULL,
		    object = NULL)

\method{summary}{FeatureMatchAnalyzer}(object, ...)

\method{plot}{FeatureMatchAnalyzer}(x, ..., type = c("all", "ph", "med", "msd", 
    "fom", "minsep", "cent.dist", "angle.diff", "area.ratio",
    "int.area", "bearing", "bdelta", "haus"))

\method{print}{FeatureMatchAnalyzer}(x, ...)

FeatureComps(Y, X, which.comps=c("cent.dist", "angle.diff", "area.ratio", "int.area",
    "bdelta", "haus", "ph", "med", "msd", "fom", "minsep", "bearing"),
    sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", deg = TRUE,
    aty = "compass", loc = NULL, ...)

\method{distill}{FeatureComps}(x, ...)
}
\arguments{
  \item{x,y,matches}{
\code{x}, \code{y} and \code{matches} are list objects with components as output by \code{deltamm} or similar function.  Only one is used, and it first checks for \code{matches}, then \code{y}, and finally \code{x}.  It expects a component named \code{mm.new.labels} that gives the number of matched objects.  In the case of the \code{plot} and \code{print} method functions, \code{x} is a list object as returned by \code{FeatureMatchAnalyzer}.

\code{distill}: output from \code{FeatureComps}.
}
  \item{X,Y}{list object giving a pixel image as output from \code{solutionset} from package \pkg{spatstat} for the verification and forecast fields, resp.  These arguments are passed directly to the \code{locperf} function.}
  \item{object}{
list object returned of class \dQuote{FeatureMatchAnalyzer}, this is the returned value from the self-same function.
}
  \item{which.comps,type}{
character vector indicating which properties of the features are to be analyzed (\code{which.comps}) or plotted (\code{type}).
}
  \item{sizefac}{
single numeric by which area calculations should be multiplied in order to get the desired units.  If unity (default) results are in terms of grid squares.
}
  \item{alpha}{
numeric value for the FOM measure (see the help file for \code{locperf}.
}
  \item{k}{
numeric indicating which quantile to use if the partial Hausdorff measure is to be used.
}
  \item{p}{
numeric giving the value of the parameter p for the Baddeley metric.
}
  \item{c}{
numeric giving the cut-off value for the Baddeley metric.
}
  \item{distfun}{
character naming a distance functions to use in calculating the various binary image measures.  Default is Euclidean distance.
}
  \item{deg, aty}{optional arguments to the \code{bearing} function.}
    \item{loc}{two-column matrix giving location coordinates for centroid distance.  If NULL, uses an indices based on the dimension of the field.}
  \item{\dots}{
Additional arguments to \code{deltametric} from package \pkg{spatstat}.  In case of the \code{summary} method function, additional optional arguments may be passed, which include \code{silent} (logical, should the information be printed to the screen or not?), \code{interest} (numeric vector defining an interest value for calculating total interest for each matched object, if NULL, this is not performed), \code{con} (name of function that takes three arguments, the first two are matrices whose rows are objects and columns are matched feature properties, where the former is a matrix of matched feature property values (e.g., angle difference) and the latter is a matrix of interest values determined by the \code{interest} argument (whereby each row is identical), the third argument to \code{con} must be called \code{which.comps}, and it gives the short-form feature property names (i.e., same as which.comps argument); see details section).  In the case of the \code{plot} method function, these are optional arguments to the function \code{barplot}.

Not used by \code{distill}.
}
}
\details{
\code{FeatureMatchAnalyzer} operates on objects of class \dQuote{matched}.  It is set up to calculate the values discussed in sec. 4 of Davis et al. (2006) for a single verification set (i.e., mean and standard deviation are not computed because it is only a single case).  If criteria is 1, then features separated by a distance D < the sum of the sizes of the two features (size of a feature is defined as the square root of its area) are considered a match.  If criteria is 2, then a match is made if D < the average of the sizes of the two features.  Finally, criteria 3 decides a match as being anything less than a pre-determined constant.

\code{FeatureComps} is the primary function called by \code{FeatureMatchAnalyzer}, and is designed as a more stand-alone type of function.  Several of the measures that can be calculated are simply the binary image measures/metrics available via, e.g., \code{locperf}.  It calculates comparisons between two matched features (i.e., between the verification and forecast fields).

\code{distill} reduces a \dQuote{FeatureComps} list object to a named numeric vector containing (in this order) the components that exist from "cent.dist", "angle.diff", "area.ratio", "int.area", "bdelta", "haus", "ph", "med", "msd", "fom", and "minsep".  This is used, for example, by \code{interester}, which is why the order is important.

The \code{summary} method function for \code{FeatureMatchAnalyzer} allows for passing a function, con, to determine confidence for each interest value.  The idea being to set the interest to zero when the particular interest value does not make sense.  For example, angle difference makes no sense if both objects are circles.  Currently, no functions are included in this package for actually doing this, and so the functionality itself has not been tested. 

The \code{print} method function for \code{FeatureMatchAnalyzer} first converts the object to a simple named matrix, then prints the matrix out.  The resulting matrix is returned invisibly.

}
\value{
FeatureMatchAnalyzer returns a list of list objects.  The specific components depend on the 'which.comps' argument, and are the same as those returned by FeatureComps.  These can be any of the following.
  \item{cent.dist }{numeric giving the centroid (Euclidean) distance.}
  \item{angle.diff }{numeric giving the orientation (major axis) angle difference.}
  \item{area.ratio}{numeric giving the area ratio, which is always between 0 and 1 because this is defined by Davis et al. (2006) to be the area of the smaller feature divided by that of the larger feature regardless of which field the feature belongs to.}
   \item{int.area}{numeric giving the intersection area of the features.}
   \item{bdelta}{numeric giving Baddeley's delta metric between the two features.}
   \item{haus, ph, med, msd, fom, minsep}{numeric, see locperf for specific information.}
   \item{bearing}{numeric giving the bearing from the forecast object centroid to the observed object centroid.}
The summary method for FeatureMatchAnalyzer invisibly returns a matrix with the same information, but where each matched object is a row and each column is the specific statistic.  Or, if optional interest argument is passed, a list with components:

\code{print} returns a named vector invisibly.
}
\references{
Davis, C. A., Brown, B. G. and Bullock, R. G. (2006) Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.
}
\author{
Eric Gilleland
}
\seealso{
Functions to identify features: \code{\link{FeatureFinder}}

Functions to merge and/or match objects: \code{\link{deltamm}}, \code{\link{centmatch}}, \code{\link{MergeForce}}

Functions to compute feature properties: \code{\link{locperf}}, \code{\link{deltametric}}, \code{\link{bearing}}

Function to calculate fuzzy logic interest values: \code{\link{interester}}
}
\examples{
data( "ExampleSpatialVxSet" )

x <- ExampleSpatialVxSet$vx
xhat <- ExampleSpatialVxSet$fcst

hold <- make.SpatialVx( x, xhat, field.type="Example",
    units = "units", data.name = "Example", 
    obs.name = "x", model.name = "xhat" )

look <- FeatureFinder(hold, smoothpar=1.5)
look2 <- centmatch(look)

tmp <- FeatureMatchAnalyzer(look2)
tmp
summary(tmp)
plot(tmp)


}
\keyword{ math }
