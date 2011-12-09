\name{FQI}
\alias{FQI}
\alias{UIQI}
\alias{ampstats}
\alias{summary.fqi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions for calculating the Forecast Quality Index
}
\description{
Functions for calculating the Forecast Quality Index (FQI) and its components.
}
\usage{
FQI(object, surr = NULL, ...)
UIQI(Fcst, Vx, ...)
ampstats(Fcst, Vx, only.nonzero = FALSE)
\method{summary}{fqi}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
list object as returned by \code{locmeasures2dPrep}.  In the case of \code{summary.fqi}, \code{object} is the list object returned by \code{FQI}.
}
  \item{Fcst,Vx}{numeric matrices giving the fields for the verification set.}
  \item{surr}{
three-dimesnional array containing surrogate fields for \code{Vx}, e.g. as returned by \code{surrogater2d}.  If NULL, these will be calculated using \code{surrogater2d}.
}
   \item{only.nonzero}{logical, should the means and variances of only the non-zero values of the fields be calculated (if so, the covariance is returned as NA)?}
  \item{\dots}{
In the case of \code{FQI}, additional arguments to \code{surrogater2d}.  Only used if \code{surr} is NULL.  In the case of \code{UIQI}, additional arguments to \code{ampstats}.  In the case of \code{summary.fqi}, these are not used.
}
}
\details{
The FQI was proposed as a spatial verification metric (a true metric in the mathematical sense) by Venugopal et al. (2005) to combine amplitude and displacement error information in a single summary statistic.  It is given by

FQI = (PHD_k(Vx, Fcst)/mean( PHD_k(Vx, surr_i); i in 1 to number of surrogates)) / (brightness * distortion)

where the numerator is a normalized partial Hausdorff distance (see help file for locperf), brightness (also called bias) is given by 2*(mu1*mu2)/(mu1^2+mu2^2), where mu1 (mu2) is the mean value of Vx (Fcst), and the distortion term is given by 2*(sig1*sig2)/(sig1^2+sig2^2), where sig1^2 (sig2^2) is the variance of Vx (Fcst) values.  The denominator is a modified UIQI (Universal Image Quality Index; Wang and Bovik, 2002), which itself is given by

UIQI = cor(Vx,Fcst)*brightness*distortion.

Note that if only.nonzero is TRUE in the call to UIQI, then the modified UIQI used in the FQI formulation is returned (i.e., without multiplying by the correlation term).
}
\value{
FQI returns a list with components:
   \item{prep.object}{character naming the locmeasures2dPrep object used in the call.}
   \item{phd.norm }{matrix of normalized partial Hausdorff distances for each value of k (rows) and each threshold (columns).}
   \item{uiqi.norm }{numeric vector of modified UIQI values for each threshold.}
   \item{fqi}{matrix of FQI values for each value of k (rows) and each threshold (columns).}

UIQI returns a list with components:
   \item{Fcst.name,Vx.name}{characters giving the names of the two fields.}
   \item{cor}{single numeric giving the correlation between the two fields.}
   \item{brightness.bias}{single numeric giving the brightness (bias) value.}
   \item{distortion.variability}{single numeric giving the distortion (variability) value.}
   \item{UIQI}{single numeric giving the UIQI (or modified UIQI if only.nonzero is set to TRUE) value.}

ampstats returns a list object with components:
   \item{mean.fcst,mean.vx}{single numerics giving the mean of Fcst and Vx, resp.}
   \item{var.fcst,var.vx}{single numerics giving the variance of Fcst and Vx, resp.}
   \item{cov}{single numeric giving the covariance between Fcst and Vx (if only.nonzero is TRUE, this will be NA).}
}
\references{
Venugopal, V., S. Basu, and E. Foufoula-Georgiou, 2005: A new metric for comparing precipitation patterns with an application to ensemble forecasts.  \emph{J. Geophys. Res.}, \bold{110}, D08111, doi:10.1029/2004JD005395, 11pp.

Wang, Z. and A. C. Bovik, 2002: A universal image quality index.  \emph{IEEE Signal Process. Lett.}, \bold{9}, 81--84.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{locperf}}, \code{link{surrogater2d}}, \code{\link{locmeasures2d}}, \code{\link{locmeasures2dPrep}}
}
\examples{
# Simulate a precip-like field.
grid<- list( x= seq( 0,5,,100), y= seq(0,5,,100))
obj<-Exp.image.cov( grid=grid, theta=.5, setup=TRUE)
look<- sim.rf( obj)
look[ look < 0] <- 0
look <- zapsmall( look)

look2 <- sim.rf( obj)
look2[ look2 < 0] <- 0
look2 <- zapsmall( look2)

# Now, find surrogates of the simulated field.
look3 <- surrogater2d( look, zero.down=TRUE, n=10)

u <- cbind( quantile( c(look2), c(0.75, 0.9)), quantile( c(look), c(0.75, 0.9)))
hold <- locmeasures2dPrep("look2", "look", thresholds=u, k=c(4, 0.75))
FQI( hold, surr=look3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
