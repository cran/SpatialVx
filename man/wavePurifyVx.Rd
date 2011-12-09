\name{wavePurifyVx}
\alias{wavePurifyVx}
\alias{plot.wavePurifyVx}
\alias{summary.wavePurifyVx}
\title{
Apply traditional forecast verification after wavelet denoising
}
\description{
Apply traditional forecast verification after wavelet denoising ala Briggs and Levine (1997).
}
\usage{
wavePurifyVx(x, y, object = NULL, climate=NULL, which.stats = c("bias", "ts", "ets", "pod", "far", "f", "hk", "mse"), thresholds = NULL, return.fields = FALSE, verbose = FALSE, ...)
\method{plot}{wavePurifyVx}(x, ...)
\method{summary}{wavePurifyVx}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y}{
For \code{wavePurifyVx}, m X n matrices giving the verification and forecast fields, resp.  If either is not supplied, then \code{object} must be supplied, and the fields given by \code{object} are used instead (even if one of \code{x} or \code{y} are supplied.  If one of these is missing, and \code{object} not supplied, the function will error out.  For \code{plot.wavePurifyVx}, list object as output from \code{wavePurifyVx}.
}
  \item{object}{
For \code{wavePurifyVx}, a list object as returned by \code{hoods2dPrep}.  Not used if both \code{x} and \code{y} are given.  For \code{summary.wavePurifyVx}, list object as returned by \code{wavePurifyVx}.
}
  \item{climate}{m X n matrix defining a climatology field.  If not NULL, then the anamoly correlation coefficient will be applied to the wavelet denoised fields.}
  \item{which.stats}{
character describing which traditional verification statistics to calculate on the wavelet denoised fields.  This is the argument passed to the argument of the same name in \code{vxstats}.
}
  \item{thresholds}{
If \code{object} is used, then this argument is ignored.  Otherwise, it is used to define events for all of the verification statistics except MSE.  However, if supplied or other statistics are to be computed, then MSE will be calculated for the fields at values >= \code{thresholds}.  If only MSE is to be computed, and \code{thresholds} is NULL, then no thresholding is applied.  If NULL, and any of the statistics besides MSE are to be calculated, then default values of the 0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9 and 0.95 quantiles (for each field, so that the thresholds may differ between fields) are used.  The same holds for anamoly correlation coefficient.
}
  \item{return.fields}{logical, should the denoised fields be returned (e.g., for subsequent plotting)?}
  \item{verbose}{logical, should progress information (including total run time) be printed to the screen?}
  \item{\dots}{
For \code{wavePurifyVx}, optional additional arguments to \code{denoise.dwt.2d} (or \code{denoise.modwt.2d}) from package \pkg{waveslim}.  Note that if the argument \code{J} is not passed, then it will be determined as J=log2(min(m,n)).  If the fields are dyadic, then the usual DWT is used, otherwise the maximal overlap DWT is used instead.  For the plot and summary method functions, these are not used.
}
}
\details{
If the fields are dyadic, then the denoise.dwt.2d function from package waveslim is applied to each field before calculating the chosen verification statistics.  Otherwise denoise.modwt.2d from the same package is used.  The result is that high-frequency fluctuations in the two fields are removed before calculating verification statistics so that the resulting statistics are less susceptible to small-scale errors (see Briggs and Levine, 1997).  See Percival and Guttorp (1994) and Lindsay et al. (1996) for more on this type of wavelet analysis including maximal overlap DWT.
}
\value{
A list object of class "wavePurifyVx" is returned with possible components (depending on what is supplied in the arguments, etc.):
  \item{X2,Y2}{m X n matrices of the denoised verification and forecast fields, resp. (only if 'return.fields' is TRUE).}
  \item{thresholds }{q X 2 matrix of thresholds applied to the forecast (first column) and verification (second column) fields, resp.  If 'climate' is not NULL, then the same thresholds for the forecast field are applied to the climatology.}
  \item{qs }{If object and thresholds are NULL, and statistics other than MSE or ACC are desired, then this will be created along with the thresholds, and is just a character version of the trhesholds. }
  \item{args}{list object containing all the optional arguments passed into ..., and the value of J used (e.g., even if not passed into ...).}
  \item{bias,ts,ets,pod,far,f,hk,mse,acc}{numeric vectors of length q (i.e., the number of thresholds) giving the associated verification statistics.}
}
\references{
Briggs, W. M. and R. A. Levine, 1997. Wavelets and field forecast verification. \emph{Mon. Wea. Rev.}, \bold{125}, 1329--1341.

Lindsay, R. W., D. B. Percival, and D. A. Rothrock, 1996.  The discrete wavelet transform and the scale analysis of the surface properties of sea ice.
  \emph{IEEE Transactions on Geoscience and Remote Sensing}, \bold{34} (3), 771--787.

Percival, D. B. and P. Guttorp, 1994.  Long-memory processes, the Allan variance and wavelets.  In \emph{Wavelets in Geophysics}, E. Foufoula-Georgiou
 and P. Kumar, Eds., New York: Academic, pp. 325--343.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{denoise.dwt.2d}}, \code{\link{denoise.modwt.2d}}, \code{\link{waverify2d}}, \code{\link{mowaverify2d}}, \code{\link{waveIS}}, \code{\link{hoods2dPrep}}, \code{\link{vxstats}}
}
\examples{
data(UKobs6)
data(UKfcst6)
data(UKloc)
hold <- surrogater2d(UKobs6, n=1, maxiter=50, verbose=TRUE)
hold <- matrix(hold, 256, 256)
UKobj <- hoods2dPrep("UKfcst6", "UKobs6", thresholds=c(0.1,2,5,10), loc=UKloc, units="mm/h")
look <- wavePurifyVx(object=UKobj, climate=hold, return.fields=TRUE, verbose=TRUE)
plot(look)
summary(look)

\dontrun{
data(pert004)
data(pert000)
hold <- hoods2dPrep( "pert004", "pert000", thresholds=c(1,10,50), loc=ICPg240Locs, units="mm/h")
look <- wavePurifyVx(object=hold, return.fields=TRUE, verbose=TRUE)
summary( look)
}

}
\keyword{ math }
