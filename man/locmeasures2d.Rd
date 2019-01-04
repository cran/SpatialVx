\name{locmeasures2d}
\alias{locmeasures2d}
\alias{locmeasures2d.default}
\alias{locmeasures2d.SpatialVx}
\alias{summary.locmeasures2d}
\alias{print.locmeasures2d}
\title{
Binary Image Measures
}
\description{
Calculate some binary image measures between two fields.
}
\usage{
locmeasures2d(object, which.stats = c("bdelta", "haus", "qdmapdiff", 
    "med", "msd", "ph", "fom"), distfun = "distmapfun", distfun.params = NULL, 
    k = NULL, alpha = 0.1, bdconst = NULL, p = 2, ...)

\method{locmeasures2d}{default}(object, which.stats = c("bdelta", "haus", "qdmapdiff", 
    "med", "msd", "ph", "fom"), distfun = "distmapfun", distfun.params = NULL, 
    k = NULL, alpha = 0.1, bdconst = NULL, p = 2, ..., Y, thresholds=NULL)

\method{locmeasures2d}{SpatialVx}(object, which.stats = c("bdelta", "haus", "qdmapdiff", 
    "med", "msd", "ph", "fom"), distfun = "distmapfun", distfun.params = NULL, 
    k = NULL, alpha = 0.1, bdconst = NULL, p = 2, ..., time.point = 1, 
    obs = 1, model = 1)

\method{print}{locmeasures2d}(x, ...)

\method{summary}{locmeasures2d}(object, ...)
}
\arguments{
  \item{object}{
For \code{locmeasures2d}, an object of class \dQuote{SpatialVx} or a valid matrix in which case \code{Y} must be explicitly provided.  For \code{summary} method function, a list object output from \code{locmeasures2d}.
}
    \item{x}{returned object from \code{locmeasures2d}.}
  \item{which.stats}{
character vector telling which measures should be calculated.
}
  \item{distfun}{character naming a function to calculate the shortest distances between each point x in the grid and the set of events.  Default is the Euclidean distance metric. Must take \code{x} as an argument, which is the event field for which the distances are to be calculated.  Must return a matrix of the same dimension as x.}
   \item{distfun.params}{list with named components giving any additional arguments to the \code{distfun} function.}
    \item{k}{
numeric vector for use with the partial Hausdorff distance.  For k that are whole numerics or integers >= 1, then the k-th highest value is returned by \code{locmeasures2d}.  If 0 <= k < 1, then the corresponding quantile is returned.
}
  \item{alpha}{numeric giving the alpha parameter for Pratt's Figure of Merit (FOM).  See the help file for \code{locperf} for more details.}
  \item{bdconst}{
numeric giving the cut-off value for Baddeley's delta metric.
}
  \item{p}{
numeric vector giving one or more values for the parameter p in Baddeley's delta metric.  Usually this is just 2.
}
    \item{Y}{m X n matrix giving the forecast field.}
    \item{thresholds}{numeric or two-column matrix giving the threshold to be applied to the verification (column one) and forecast (column two) fields.  If a vector, same thresholds are applied to both fields.}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
  \item{\dots}{
optional arguments to \code{deltametric} and \code{distmap} from package spatstat.  Not used by the \code{summary} or \code{print} methods here.  }
}
\details{
It is useful to introduce some notation.  Let d(x,A) be the shortest distance from a point x, anywhere in the grid, to a set A contained in the grid.  Here, Euclidean distance is used (default) for d(x,A), but note that some papers (e.g., Venugopal et al., 2005) use other distances, such as the taxi-cab distance (use \code{distfun} argument to change the distance method).

The Hausdorff distance between two sets A and B contained in the finite grid is given by max( max( d(x,A), x in B), max( d(x,B), x in A)), and can be re-written as H(A,B) = max( abs( d(x,A) - d(x,B))), where x is taken over all points in the grid.  Several of the distances here are modifications of the Hausdorff distance.  The Baddeley metric, for example, is the Lp norm of abs( w(d(x,A)) - w(d(x,B))), where again x is taken from over the entire grid, and w is any concave continuous function that is strictly increasing at zero.  Here, w(t) = min( t, c), where c is some constant given by the \code{bdconst} argument.

Calculates one or more of the following binary image measures:

\dQuote{bdelta} Baddeley delta metric (Baddeley, 1992a,b; Gilleland, 2011; Schwedler and Baldwin, 2011)

\dQuote{haus} Hausdorff distance (Baddeley, 1992b; Schwedler and Baldwin, 2011)

\dQuote{qdmapdiff} Quantile (or rank) of the differences in distance maps.  See the help file for \code{locperf}.

\dQuote{med} Mean Error Distance (Peli and Malah, 1982; Baddeley, 1992a).  See the help file for \code{locperf}.

\dQuote{msd} Mean Square Error Distance (Peli and Malah, 1982; Baddeley, 1992a).  See the help file for \code{locperf}.

\dQuote{ph} Partial Hausdorff distance.  See the help file for \code{locperf}.

\dQuote{fom} Pratt's Figure of Merit (Peli and Malah, 1982; Baddeley, 1992a, Eq (1)).  See the help file for \code{locperf}.

These distances are summaries in and of themselves, so the summary method function simply displays the results in an easy to read manner.
}
\value{
A list with at least one of the following components depending on the argument which.stats
\item{bdelta }{p by q matrix giving the Baddeley delta metric for each desired value of p (rows) and each threshold (columns)}
\item{haus}{numeric vector giving the Hausdorff distance for each threshold}
\item{qdmapdiff}{k by q matrix giving the difference in distance maps for each of the k-th largest value(s) or quantile(s) (rows) for each threshold (columns).}
\item{medMiss, medFalseAlarm, msdMiss, msdFalseAlarm}{two-row matrix giving the mean error (or square error) distance as (Forecast, Observation) or misses and (Observation, Forecast) or false alarms.}
\item{ph}{k by q matrix giving the k-th largest value(s) or quantile(s) (rows) for each threshold (columns) of the maximum between the distances from one field to the other.}
\item{fom}{numeric vector giving Pratt's figure of merit.}
}
\references{
Baddeley, A. (1992a)  An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, W. Forstner and S. Ruwiedel, Eds., Wichmann, 59--78.

Baddeley, A. (1992b)  Errors in binary images and an Lp version of the Hausdorff metric.  \emph{Nieuw Arch. Wiskunde}, \bold{10}, 157--183.

Gilleland, E. (2011)  Spatial forecast verification: Baddeley's delta metric applied to the ICP test cases.  \emph{Wea. Forecasting}, \bold{26}, 409--415, doi:10.1175/WAF-D-10-05061.1.  

Peli, T. and Malah, D. (1982) A study on edge detection algorithms.  \emph{Computer Graphics and Image Processing}, \bold{20}, 1--21.

Schwedler, B. R. J. and Baldwin, M. E. (2011)  Diagnosing the sensitivity of binary image measures to bias, location, and event frequency within a forecast verification framework.  \emph{Wea. Forecasting}, \bold{26}, 1032--1044, doi:10.1175/WAF-D-11-00032.1.

Venugopal, V., Basu, S. and Foufoula-Georgiou, E. (2005) A new metric for comparing precipitation patterns with an application to ensemble forecasts.  \emph{J. Geophys. Res.}, \bold{110}, D08111, doi:10.1029/2004JD005395, 11pp.
}
\author{
Eric Gilleland
}
\note{
Binary fields are determined by having values >= the thresholds.
}


\seealso{
\code{\link{deltametric}}, \code{\link{distmap}}, \code{\link{make.SpatialVx}}
}
\examples{
x <- y <- matrix(0, 10, 12)
x[2,3] <- 1
y[4,7] <- 1

hold <- make.SpatialVx(x, y, thresholds = 0.1,
    field.type = "random", units = "grid squares")
locmeasures2d(hold, k = 1)

# Alternatively ...
locmeasures2d(x, thresholds = 0.1, k = 1, Y = y)

\dontrun{
data( "geom000" )
data( "geom001" )
data( "ICPg240Locs" )

hold <- make.SpatialVx( geom000, geom001, thresholds = c(0.1, 50.1), 
    projection = TRUE, map=TRUE, loc = ICPg240Locs, loc.byrow = TRUE,
    field.type = "Precipitation", units = "in/100",
    data.name= "ICP Geometric Cases", obs.name = "geom000",
    model.name = "geom001" )

hold2 <- locmeasures2d(hold, k=c(4, 0.975), alpha=c(0.1,0.9))
summary(hold2)
}
}

\keyword{math}
