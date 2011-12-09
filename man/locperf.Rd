\name{locperf}
\alias{locperf}
\alias{distob}
\alias{distmapfun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Localization performance measures
}
\description{
Some localization performance (distance) measures for binary images.
}
\usage{
locperf(X, Y, which.stats = c("ph", "mhd", "med", "msd", "fom", "minsep"), alpha = 0.1, k = 4, distfun = "distmapfun", ...)
distob(X, Y, distfun = "distmapfun", ...)
distmapfun(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
list object giving a pixel image as output from \code{solutionset} from package \pkg{spatstat}.  This corresponds to the set B in the Details section below.
}
  \item{Y}{
list object giving a pixel image as output from \code{solutionset} from package \pkg{spatstat}.  This corresponds to the set A in the Details section below.
}
   \item{x}{list object of class "owin" as returned by \code{solutionset} from package \pkg{spatstat}.}
  \item{which.stats}{
character vector stating which localization performance measure to calculate.
}
  \item{alpha}{
numeric giving the scaling constant for Pratt's figure of merit (FOM).  Only used for \code{which.stat} method "fom".
}
  \item{k}{
single numeric giving the order for the partial Hausdorff distance.  If 0 <= k < 1, this is assumed to be a quantile for use with the \code{quantile} function.  Otherwise, k should be a whole number such that 1 <= k <= Nxy, where nxy is the total number of grid points in the set.
}
  \item{distfun}{character specifying a distance metric that returns a matrix of same dimension as \code{X} yielding, at each point x, the shortest distances from x to the set of events in the field.  Default is \code{distmapfun}, which returns the Euclidean distances.}
  \item{\dots}{
Optional arguments to the \code{distfun} function.  In the case of \code{distmapfun}, these are the optional arguments to \code{distmap} from package \pkg{spatstat}.
}
}
\details{
This function computes localization performance (or distance) measures detailed in Peli and Malah (1982) and Baddeley (1992), as well as a modification of one of these distances detailed in Zhu et al. (2011); distob.

First, it is helpful to establish some notation.  Suppose a distance rho(x,y) is defined between any two pixels x and y in the entire raster of pixels/grid (If distfun is distmapfun (default), then rho is the Euclidean distance) that satisfies the formal mathematical axioms of a metric.  Let d(x,A) denote the shortest distance (smallest value of rho) from the point x in the entire raster to the the set A contained in the raster.  That is, d(x,A) = min(rho(x,a): a in A contained in the raster) [formally, the minimum should be the infimum], with d(x, empty set) defined to be infinity.  Note that the distfun argument is a function that returns d(x,A) for all x in the raster.

The mean error distance ("med") is the mean of d(x,A) over the points in B.  That is e.bar = mean( d(x,A)), over all x in B.

The mean square error distance ("msd") is the mean of the squared d(x,A) over the points in B.  That is, e2.bar = mean( d(x,A)^2), over all x in B.

Pratt's figure of merit ("fom") is given by: FOM(A,B) = sum( 1/(1+alpha*d(x,A)^2))/max(N(A),N(B)), where x in B, and N(A) (N(B)) is the number of points in the set A (B) and alpha is a scaling constant (see, e.g., Pratt, 1977; Abdou and Pratt, 1979).  The scaling constant is typically set to 1/9 when rho is normalized so that the smallest nonzero distance between pixel neighbors is 1.  The default (0.1) here is approximately 1/9.

Minimum separation distance between boundaries ("minsep") is just the smallest value of the distance map of one field over the subset where events occur in the other.  This is mainly for when single features within the fields are being compared.

distob is a modification of the mean error distance where if there are no events in either field, the value is 0, and if there are no events in one field only, the value is something large (in this case the length of the longest side of the grid).

The Hausdorff distance for a finite grid is given by max( max( d(x,B); x in A), max( d(x,A); x in B)), and can be written as max( abs(d(x,A) - d(x,B)), over all x in the raster).  The partial Hausdorff distance ("ph") is a modification that replaces the maximum in the latter equation with a k-th order statistic (or quantile).  The modified Hausdorff distance ("mhd") is given by mhd(A,B) = max( e.bar(A,B), e.bar(B,A)).  See, e.g., Baddeley, 1992; Schwedler and Baldwin (2011).

For computational efficiency, the distance transform method is used via \code{distmap} from package \pkg{spatstat} for calculating d(x,A) x in the raster.
}
\value{
\code{locperf} returns a list object with components depending on \code{which.stats}: one or more of the following, each of which is a single numeric.
\item{ph }{numeric giving the value of the partial }
\item{comp2 }{Description of 'comp2'}

\code{distob} returns a single numeric.

\code{distmapfun} returns a matrix of same dimension as the input argument's field.
}
\references{
Abdou, I. E. and W. K. Pratt, 1979: Quantitative design and evaluation of enhancement/thresholding edge detectors.  \emph{Proc. IEEE}, \bold{67}, 753--763.

Baddeley, A., 1992.  An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, W. Forstner and S. Ruwiedel, Eds., Wichmann, 59--78.

Peli, T. and D. Malah, 1982: A study on edge detection algorithms.  \emph{Computer Graphics and Image Processing}, \bold{20}, 1--21.

Pratt, W. K., 1977: \emph{Digital Image Processing}.  John Wiley and Sons, New York.

Schwedler, B. R. J. and M. E. Baldwin, 2011.  Diagnosing the sensitivity of binary image measures to bias, location, and event frequency within a forecast verification framework.  \emph{Wea. Forecasting}, \bold{26}, 1032--1044, doi: 10.1175/WAF-D-11-00032.1.

Zhu, M., V. Lakshmanan, P. Zhang, Y. Hong, K. Cheng, and S. Chen, 2011: Spatial verification using a true metric.  \emph{Atmos. Res.}, \bold{102}, 408--419, doi:10.1016/j.atmosres.2011.09.004.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{distmap}}, \code{\link{solutionset}}, \code{\link{im}}, \code{\link{bounding.box}}, \code{\link{as.rectangle}}, \code{\link{metrV}}, \code{\link{locmeasures2d}}, \code{\link{locmeasures2dPrep}}
}
\examples{
x <- y <- matrix( 0, 10, 12)
x[2,3] <- 1
y[4,7] <- 1
x <- im( x)
y <- im( y)
x <- solutionset( x > 0)
y <- solutionset( y > 0)
locperf( x, y)
par( mfrow=c(1,2))
image.plot( distmapfun(x))
image.plot( distmapfun(y))
}
\keyword{ math }
