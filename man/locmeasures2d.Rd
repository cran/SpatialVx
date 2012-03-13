\name{locmeasures2d}
\alias{locmeasures2d}
\alias{summary.locmeasures2d}
\title{
Binary Image Measures
}
\description{
Calculate some binary image measures between two fields.
}
\usage{
locmeasures2d(object, which.stats = c("baddeley", "hausdorff", "ph", "mhd", "med", "msd", "fom"), distfun = "distmapfun", distfun.params = NULL, ...)
\method{summary}{locmeasures2d}(object, ...)
}
\arguments{
  \item{object}{
list object returned from \code{locmeasures2dPrep} for \code{locmeasures2d} function, and a list object output from \code{locmeasures2d} for the summary method.
}
  \item{which.stats}{
character vector telling which measures should be calculated.
}
  \item{distfun}{character naming a function to calculate the shortest distances between each point x in teh grid and the set of events.  Default is the Euclidean distance metric. Must take x as an argument, which is the event field for which the distances are to be calculated.  Must return a matrix of the same dimension as x.}
   \item{distfun.params}{list with named components giving any additional arguments to the \code{distfun} function.}
  \item{\dots}{
optional arguments to \code{deltametric} and \code{distmap} from package spatstat.  Not used by the summary method here.
}
}
\details{
It is useful to introduce some notation.  Let d(x,A) be the shortest distance from a point x, from the entire grid, to a set A contained in the grid.  Here, Euclidean distance is used for d(x,A), but note that some papers (e.g., Venugopal et al., 2005) use other distances, such as the taxi-cab distance.

The Hausdorff distance between two sets A and B contained in the finite grid is given by max( max( d(x,A), x in B), max( d(x,B), x in A)), and can be re-written as H(A,B) = max( abs( d(x,A) - d(x,B))), where x is taken over all points in teh grid.  Several of the distances here are modifications of the Hausdorff distance.  The Baddeley metric, for example, is the Lp norm of abs( w(d(x,A)) - w(d(x,B))), where again x is taken from over the entire grid, and w is any concave continuous function that is strictly increasing at zero.  Here, w(t) = min( t, c), where c is some constant given by the \code{bdconst} argument.

Calculates one or more of the following binary image measures:

"baddeley" Baddeley delta metric (Baddeley, 1992a,b; Gilleland, 2011; Schwedler and Baldwin, 2011)

"hausdorff" Hausdorff distance (Baddeley, 1992b; Schwedler and Baldwin, 2011)

"ph" Partial Hausdorff distance (Schwedler and Baldwin, 2011).  See the help file for \code{locperf}.

"mhd" Modified Hausdorff Distance (Schwedler and Baldwin, 2011).  See the help file for \code{locperf}.

"med" Mean Error Distance (Peli and Malah, 1982; Baddeley, 1992a).  See the help file for \code{locperf}.

"msd" Mean Square Error Distance (Peli and Malah, 1982; Baddeley, 1992a).  See the help file for \code{locperf}.

"fom" Pratt's Figure of Merit (Peli and Malah, 1982; Baddeley, 1992a, Eq (1)).  See the help file for \code{locperf}.

These distances are summaries in and of themselves, so the summary method function simply displays the results in an easy to read manner.
}
\value{
A list with at least one of the following components depending on the argument which.stats
\item{baddeley }{p X q matrix giving the Baddeley delta metric for each desired value of p (rows) and each threshold (columns)}
\item{hausdorff}{numeric vector giving the Hausdorff distance for each threshold}
\item{ph}{k X q matrix giving the partial Hausdorff distance for each k-th largest value (rows) and each threshold (columns)}
\item{mhd}{numeric vector giving the modified Hausdorff distance for each threshold.}
\item{med}{numeric vector giving the mean error distance.}
\item{msd}{numeric vector giving the mean square error distance.}
\item{fom}{numeric vector giving Pratt's figure of merit.}
}
\references{
Baddeley, A., 1992a.  An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, W. Forstner and S. Ruwiedel, Eds., Wichmann, 59--78.

Baddeley, A., 1992b.  Errors in binary images and an Lp version of the Hausdorff metric.  \emph{Nieuw Arch. Wiskunde}, \bold{10}, 157--183.

Gilleland, E., 2011.  Spatial forecast verification: Baddeley's delta metric applied to the ICP test cases.  \emph{Wea. Forecasting}, \bold{26}, 409--415, DOI: 10.1175/WAF-D-10-05061.1.  

Peli, T. and D. Malah, 1982: A study on edge detection algorithms.  \emph{Computer Graphics and Image Processing}, \bold{20}, 1--21.

Schwedler, B. R. J. and M. E. Baldwin, 2011.  Diagnosing the sensitivity of binary image measures to bias, location, and event frequency within a forecast verification framework.  Accepted to \emph{Wea. Forecasting}, doi: 10.1175/WAF-D-11-00032.1.

Venugopal, V., S. Basu, and E. Foufoula-Georgiou, 2005: A new metric for comparing precipitation patterns with an application to ensemble forecasts.  \emph{J. Geophys. Res.}, \bold{110}, D08111, doi:10.1029/2004JD005395, 11pp.
}
\author{
Eric Gilleland
}
\note{
Binary fields are determined by having values >= the thresholds.
}


\seealso{
\code{\link{deltametric}}, \code{\link{distmap}}, \code{\link{locmeasures2dPrep}}
}
\examples{
x <- y <- matrix(0, 10, 12)
x[2,3] <- 1
y[4,7] <- 1
hold <- locmeasures2dPrep("y", "x", thresholds=0.1, k=1, units="grid squares")
locmeasures2d( hold)

\dontrun{
data( geom000)
data(geom001)
hold <- locmeasures2dPrep("geom001", "geom000", thresholds=c(0.1,50.1), k=c(4,0.975), alpha=c(0.1,0.9), units="in/100")
hold2 <- locmeasures2d( hold)
summary( hold2)
}
}

\keyword{math}
