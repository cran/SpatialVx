\name{FeatureMatchAnalyzer}
\alias{FeatureMatchAnalyzer}
\alias{FeatureComps}
\title{
Analyze Features of a Verification Set
}
\description{
Analyze matched features of a verification set.
}
\usage{
FeatureMatchAnalyzer(x, y = NULL, matches = NULL, object = NULL, which.comps = c("cent.dist", "angle.diff", "area.ratio", "int.area", "bdelta", "haus", "ph", "mhd", "med", "msd", "fom", "minsep"), sizefac = 1, alpha = 0.1, k = 4, p = 2, c = Inf, distfun = "distmapfun", ...)
FeatureComps(Y, X, which.comps=c("cent.dist", "angle.diff", "area.ratio", "int.area", "bdelta", "haus", "ph", "mhd", "med", "msd", "fom", "minsep"), sizefac=1, alpha=0.1, k=4, p=2, c=Inf, distfun="distmapfun", ...)
}
\arguments{
  \item{x,y,matches}{
\code{x}, \code{y} and \code{matches} are list objects with components as output by \code{deltamm} or similar function.  Only one is used, and it first checks for \code{matches}, then \code{y}, and finally \code{x}.  It expects a component names \code{mm.new.labels} that gives the number of matched objects.
}
  \item{X,Y}{list object giving a pixel image as output from \code{solutionset} from package \pkg{spatstat} for the verification and forecast fields, resp.  These arguments are passed directly to the \code{locperf} function.}
  \item{object}{
Obligatory argument for use with \code{FeatureSuite} function.  Not used by \code{FeatureMatchAnalyzer} (so far).
}
  \item{which.comps}{
character vector indicating which properties of the features are to be analyzed.
}
  \item{sizefac}{
single numeric by which area calculations should be multiplied in order to get the desired units.  If unity (default) results are in terms of grid squares.
}
  \item{alpha}{
numeric value for teh FOM measure (see the help file for \code{locperf}.
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
  \item{\dots}{
Additional arguments to \code{deltametric}.
}
}
\details{
FeatureMatchAnalyzer is designed to be used with 'FeatureSuite'.  It is set up to calculate the values discussed in sec. 4 of Davis et al. (2006) for a single verification set (i.e., mean and standard deviation are not computed because it is only a single case).  If criteria is 1, then features separated by a distance D < the sum of the sizes of the two features (size of a feature is defined as the square root of its area) are considered a match.  If criteria is 2, then a match is made if D < the average of the sizes of the two features.  Finally, criteria 3 decides a match as being anything less than a pre-determined constant.

FeatureComps is the primary function called by FeatureMatchAnalyzer, and is designed as a more stand-alone type of function.  Several of the measures that can be calculated are simply the binary image measures/metrics available via, e.g., locperf.  It calculates comparisons between two matched features (i.e., between the verification and forecast fields).

}
\value{
FeatureMatchAnalyzer returns a list of list objects.  The specific components depend on the 'which.comps' argument, and are the same as those returned by FeatureComps.  These can be any of the following.
  \item{cent.dist }{numeric giving the centroid (Euclidean) distance.}
  \item{angle.diff }{numeric giving the orientation (major axis) angle difference.}
  \item{area.ratio}{numeric giving the area ratio, which is always between 0 and 1 because this is defined by Davis et al. (2006) to be the area of the smaller feature divided by that of the larger feature regardless of which field the feature belongs to.}
   \item{int.area}{numeric giving the intersection area of the features.}
   \item{bdelta}{numeric giving Baddeley's delta metric between the two features.}
   \item{haus, ph, mhd, med, msd, fom, minsep}{numeric, see locperf for specific information.}
}
\references{
Davis CA, BG Brown, and RG Bullock, 2006. Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{locperf}}, \code{\link{FeatureSuite}}, \code{\link{convthresh}}, \code{\link{deltamm}}, \code{\link{deltametric}}
}
\examples{
x <- y <- matrix(0, 10, 12)
x[2:3,c(3:6, 8:10)] <- 1
y[c(1:2, 9:10),c(3:6)] <- 1

hold <- FeatureSuitePrep("y", "x")
look <- convthresh( hold, smoothpar=1.5)
look2 <- centmatch(look, object=hold)
FeatureMatchAnalyzer(matches=look2)

\dontrun{
data(pert000)
data(pert004)
hold <- FeatureSuitePrep("pert004", "pert000")
look <- convthresh( hold, smoothpar=10.5)
look2 <- centmatch(look, object=hold)
FeatureMatchAnalyzer(matches=look2)
   }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
