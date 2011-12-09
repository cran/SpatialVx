\name{convthresh}
\alias{convthresh}
\alias{salIDfun}
\alias{disjointer}
\alias{plot.features}
\alias{summary.features}
\alias{plot.summary.features}
\title{
Convolution Threshold Feature Identification
}
\description{
Identify features within a field via thresholding a convolution-smoothed version of the field.
}
\usage{
convthresh(object, smoothfun = "disk2dsmooth", smoothpar = 1, smoothfunargs = NULL,
          thresh = 1e-08, idfun = "disjointer", zero.down = FALSE, ...)

salIDfun(object, fac = 0.06666667, q = 0.95, wash.out = NULL, thresh = NULL,
          idfun = "disjointer", ...)

disjointer(x, method = "C")

\method{plot}{features}(x, ...)

\method{summary}{features}(object, ...)

\method{plot}{summary.features}(x, ...)
}
\arguments{
  \item{object}{
list object of class \dQuote{FeatureSuitePrep} as returned by \code{FeatureSuitePrep}.  In the case of \code{summary.features}, this is an object returned by \code{convthresh}, \code{salIDfun} (or other similar function that returns an object of class \dQuote{features}).
}
   \item{x}{binary matrix identifying event points.  For \code{plot} method functions, this is the output object returned from either \code{convthresh} or \code{salIDfun} (or other feature identifying function that returns objects of class \dQuote{features}) or by their \code{summary} method function.}
  \item{smoothfun}{
character naming a convolution smoothing function that takes argument \code{x} (the field to be smoothed, not the same as the argument to \code{disjointer}) and a smoothing parameter as the first two arguments.  Other arguments are allowed, and are included through \code{smoothfunargs}.  Default is a disk kernel smoothing function.
}
  \item{smoothpar}{
single numeric giving the smoothing parameter for \code{smoothfun}.
}
  \item{smoothfunargs}{
list object with named additional arguments to \code{smoothfun}.
}
  \item{fac}{numeric giving a factor by which to multiply the R quantile in determining the threshold from the fields.  Default is to multiply by ~ 1/15 as in Wernli et al (2008, 2009).  This is not used if \code{thresh} is supplied.}
  \item{thresh}{
numeric vector of length one or two giving the threshold over which (and including, i.e., >=) events are defined.  If different thresholds are used for the forecast and verification fields, then the first element is the threshold for the forecast, and the second is that for the verification field.  For \code{salIDfun}, if this value is NULL, the thresholds are determined by fac*R_q, where R_q is the \code{q}-th quantile of the field (different thresholds for each field are used).
}
  \item{q}{If \code{thresh} is NULL, then this is the quantile used to deterimne the thresholds for the two fields.}
  \item{idfun}{
character naming the function used to identify (and label) individual features in the convolved and thresholded field.  Must take argument \code{x}, which is the convolved and thresholded (binary) field.
}
  \item{wash.out}{(optional) numeric giving a lower threshold over which the quantile for determining the threshold should be found.  If NULL, this is not performed.  Only used if \code{thresh} is NULL.}
  \item{zero.down}{
logical, should negative values and relatively very small values be set to zero after smoothing the fields?  For thresholds larger than such values, this argument is moot.  \code{zapsmall} is used to set the very small positive values to zero.
}
  \item{method}{ character giving the connected components algorithm to be used in function \code{connected} from \pkg{spatstat}.}
  \item{\dots}{
\code{convthresh} and \code{salIDfun}: additional arguments to \code{idfun}.

Not used by the \code{plot} method function.

The \code{summary} method function can take the argument: \code{silent}--logical, should information be printed to the screen (FALSE) or not (TRUE).
}
}
\details{
The convolution threshold technique is a simple way of identifying features (referred to as objects in image analysis circles where features are the attributes thereof) within a field.  First, the field is smoothed using a convolution smoother, and then it is set to a binary image where everything above a given threshold is set to one.  Individual features are identified through any choice of function given by \code{idfun}.  The default is to use a connected components algorithm using the \pkg{spatstat} function \code{connected}.

This is the method used by Davis et al. (2006a,b).
}
\value{
convthresh and salIDfun return a list object of class \dQuote{features} with comopnents:
  \item{data.name}{character vector naming the verification and forecast (R object) fields, resp.}
  \item{X.feats,Y.feats}{The identified features for the verification and forecast fields as returned by the idfun function.}
  \item{X.labeled,Y.labeled}{matrices of same dimension as the forecast and verification fields giving the images of the convolved and thresholded verification and forecast fields, but with each individually identified object labeled 1 to the number of objects in each field.}
  \item{identifier.function,identifier.label}{character strings naming the function and giving the long name (for use with plot method function).}
 
The plot method functions do not return anything.

The summary method function for objects of class \dQuote{features} returns a list with components:

\item{X,Y}{matrices whose rows are objects and columns are properties: centroidX and centroidY (the x- and y- coordinates for the feature centroids), area (the area of each feature in squared grid points), the orientation angle for the fitted major axis, the aspect ratio, Intensity0.25 and Intensity0.9 (the lower quartile and 0.9 quantile of intensity values for each feature).}
}
\references{
Davis, C. A., Brown, B. G. and Bullock, R. G. (2006a) Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.

Davis, C. A., Brown, B. G. and Bullock, R. G. (2006b) Object-based verification of precipitation forecasts, Part II: Application to convective rain systems. \emph{Mon. Wea. Rev.}, \bold{134}, 1785--1795.

Wernli, H., Paulat, M. Hagen, M. and Frei, C. (2008)  SAL--A novel quality measure for the verification of quantitative precipitation forecasts.  \emph{Mon. Wea. Rev.}, \bold{136}, 4470--4487.

Wernli, H., Hofmann, C. and Zimmer, M. (2009)  Spatial forecast verification methods intercomparison project: Application of the SAL technique.  \emph{Wea. Forecasting}, \bold{24}, 1472--1484, doi:10.1175/2009WAF2222271.1.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{connected}}, \code{\link{as.im}}, \code{\link{tess}}, \code{\link{tiles}}, \code{\link{owin}}, \code{\link{FeatureSuite}}, \code{\link{FeatureSuitePrep}}
}
\examples{
x <- y <- matrix(0, 10, 12)
x[2:3,c(3:6, 8:10)] <- 1
y[c(4:7, 9:10),c(7:9, 11:12)] <- 1

hold <- FeatureSuitePrep("x", "y")
look <- convthresh( hold, smoothpar=0.5)
length( look$X.feats) # two objects in x.
length( look$Y.feats) # four ovjects in y.
plot( look$X.feats)
plot( look$Y.feats)
par(mfrow=c(1,2))
plot(look) 
look2 <- summary(look)
plot(look2)

\dontrun{
data(pert000)
data(pert004)
hold <- FeatureSuitePrep("pert000", "pert004")
look <- convthresh( hold, smoothpar=3.5)
length(look$X.feats)
length(look$Y.feats)
image.plot(look$sIx)
image.plot(look$sIy)
par(mfrow=c(1,2))
plot(look)
look2 <- summary(look, silent=TRUE)
plot(look2)
   }
}
\keyword{ math }
\keyword{ manip }
\keyword{ graphs }
