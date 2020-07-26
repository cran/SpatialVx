\name{FQI}
\alias{FQI}
\alias{UIQI}
\alias{ampstats}
\alias{print.fqi}
\alias{summary.fqi}
\title{
Forecast Quality Index
}
\description{
Functions for calculating the Forecast Quality Index (FQI) and its components.
}
\usage{
FQI(object, surr = NULL, k = 4, time.point = 1, obs = 1, model = 1, ...)

UIQI(X, Xhat, ...)

ampstats(X, Xhat, only.nonzero = FALSE)

\method{print}{fqi}(x, ...)

\method{summary}{fqi}(object, ...)
}
\arguments{
  \item{object}{
list object of class \dQuote{SpatialVx}.  In the case of the \code{summary} method, \code{object} is the list object returned by \code{FQI}.
}
  \item{X,Xhat}{numeric matrices giving the fields for the verification set.}
   \item{x}{list object of class \dQuote{fqi} as returned by \code{FQI}.}
  \item{surr}{
three-dimesnional array containing surrogate fields for \code{X}, e.g. as returned by \code{surrogater2d}.  If NULL, these will be calculated using \code{surrogater2d}.
}
   \item{only.nonzero}{logical, should the means and variances of only the non-zero values of the fields be calculated (if so, the covariance is returned as NA)?}
    \item{k}{
numeric vector for use with the partial Hausdorff distance.  For k that are whole numerics or integers >= 1, then the k-th highest value is returned by \code{locmeasures2d}.  If 0 <= k < 1, then the corresponding quantile is returned.
}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which forecast model to select for the analysis.}
  \item{\dots}{
In the case of \code{FQI}, additional arguments to \code{surrogater2d}.  Only used if \code{surr} is NULL.  In the case of \code{UIQI}, additional arguments to \code{ampstats}.  In the case of \code{summary.fqi}, these are not used.
}
}
\details{
The FQI was proposed as a spatial verification metric (a true metric in the mathematical sense) by Venugopal et al. (2005) to combine amplitude and displacement error information in a single summary statistic.  It is given by

FQI = (PHD_k(X, Xhat)/mean( PHD_k(X, surr_i); i in 1 to number of surrogates)) / (brightness * distortion)

where the numerator is a normalized partial Hausdorff distance (see help file for locperf), brightness (also called bias) is given by 2*(mu1*mu2)/(mu1^2+mu2^2), where mu1 (mu2) is the mean value of X (Xhat), and the distortion term is given by 2*(sig1*sig2)/(sig1^2+sig2^2), where sig1^2 (sig2^2) is the variance of X (Xhat) values.  The denominator is a modified UIQI (Universal Image Quality Index; Wang and Bovik, 2002), which itself is given by

UIQI = cor(X,Xhat)*brightness*distortion.

Note that if \code{only.nonzero} is \code{TRUE} in the call to \code{UIQI}, then the modified UIQI used in the FQI formulation is returned (i.e., without multiplying by the correlation term).

The \code{print} method so far just calls the \code{summary} method.
}
\value{
FQI returns a list with with the following components:
   \item{phd.norm }{matrix of normalized partial Hausdorff distances for each value of k (rows) and each threshold (columns).}
   \item{uiqi.norm }{numeric vector of modified UIQI values for each threshold.}
   \item{fqi}{matrix of FQI values for each value of k (rows) and each threshold (columns).}

It will also have the same attributes as the \dQuote{SpatialVx} object with additional attributes defining the arguments specific to parameters used by the function.

UIQI returns a list with components:
   \item{data.name}{character vector giving the names of the two fields.}
   \item{cor}{single numeric giving the correlation between the two fields.}
   \item{brightness.bias}{single numeric giving the brightness (bias) value.}
   \item{distortion.variability}{single numeric giving the distortion (variability) value.}
   \item{UIQI}{single numeric giving the UIQI (or modified UIQI if only.nonzero is set to TRUE) value.}

ampstats returns a list object with components:
   \item{mean.fcst,mean.vx}{single numerics giving the mean of Xhat and X, resp.}
   \item{var.fcst,var.vx}{single numerics giving the variance of Xhat and X, resp.}
   \item{cov}{single numeric giving the covariance between Xhat and X (if only.nonzero is TRUE, this will be NA).}
}
\references{
Venugopal, V., Basu, S. and Foufoula-Georgiou, E. (2005) A new metric for comparing precipitation patterns with an application to ensemble forecasts.  \emph{J. Geophys. Res.}, \bold{110}, D08111, 11 pp., doi:10.1029/2004JD005395.

Wang, Z. and Bovik, A. C. (2002) A universal image quality index.  \emph{IEEE Signal Process. Lett.}, \bold{9}, 81--84.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{locperf}}, \code{\link{surrogater2d}}, \code{\link{locmeasures2d}}
}
\examples{
data( "ExampleSpatialVxSet" )

 x <- ExampleSpatialVxSet$vx
 xhat <- ExampleSpatialVxSet$fcst

 # Now, find surrogates of the simulated field.
 z <- surrogater2d(x, zero.down=TRUE, n=10)

 u <- list( X = cbind( quantile( c(x), c(0.75, 0.9)) ),
	Xhat = cbind( quantile( c(xhat), c(0.75, 0.9) ) ) )

hold <- make.SpatialVx(x, xhat, thresholds = u,
			field.type = "Example", units = "none",
			data.name = "ExampleSpatialVxSet",
			obs.name = "X", model.name = "Xhat" )

FQI(hold, surr = z, k = c(4, 0.75) )
}
\keyword{ math }
