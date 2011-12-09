\name{waveIS}
\alias{waveIS}
\alias{plot.waveIS}
\alias{summary.waveIS}

\title{
Intensity Scale (IS) verification
}
\description{
Intensity Scale (IS) verification based on Casat et al (2004) and Casati (2010).
}
\usage{
waveIS(obj, J = NULL, wavelet.type = "haar", verbose = FALSE)
\method{plot}{waveIS}(x, main1 = "X", main2 = "Y", which.plots = c("all", "mse", "ss", "energy"), level.label = NULL, ...)
\method{summary}{waveIS}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A list object as returned by \code{hoods2dPrep}.
}
   \item{x,object}{list object returned by \code{waveIS}.}
  \item{main1,main2}{character giving labels for the plots where \code{main1} refers to the verification field and \code{main2} to the forecast field.}
  \item{which.plots}{character vector naming one or more specific plots to do.}
  \item{level.label}{optional character vector to use for level names on the plot(s).}
  \item{J}{
numeric integer giving the number of levels to use.  If NULL and the field is dyadic, this will be log2(min(dim(X))), where X is a field from teh verification set.  If NULL and the field is not dyadic, then \code{J} is set equal to 4.  Note that if the fields are not dyadic, the function will be much slower.
}
  \item{wavelet.type}{
character giving the name of the wavelet type to use as accepted by \code{dwt.2d} and \code{modwt.2d}.
}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{\dots}{Not used by \code{plot.waveIS}.  Only sort of used by \code{summary.waveIS}.  One can put the argument silent=TRUE so that nothing is printed to the screen (useful if you just want the values calculated and stored without writing to the screen).}
}
\details{
This function applies various statistics to the detail fields (in wavelet space) of a discrete wavelet decomposition (DWT) of the binary error fields for a verification set.  In particular, the statistics described in Casati et al (2004) and Casati (2010) are calculated.  This function depends on the waverify2d or mowaverify2d function (depending on whether the fields are dyadic or not, resp.), which themselves depend on the dwt.2d and idwt.2d or modwt.2d and imodwt.2d functions.

See the references herein and the help files and references therein for dwt.2d and modwt.2d for more information on this approach, as well as Percival and Guttorp (1994) and Lindsay et al. (1996).
}
\value{
A list object of class "waveIS" that contains the entire prep object passed in by obj, as well as additional components:
  \item{EnVx,EnFcst}{J X q matrices giving the energy for the verification and forecast fields, resp., for each threshold (columns) and scale (rows).}
  \item{MSE,SS }{J X q matrices giving the mean square error and IS skill score for each threshold (column) and scale (rows).}
  \item{Bias}{numeric vector of length q giving the frequency bias of the original fields for each threshold.}
plot.waveIS does not return any value.  A plot is created on the current graphic device.
summary.waveIS returns a list invisibly with the same components as returned by waveIS along with extra components:
  \item{MSEu,SSu,EnVx.u,EnFcst.u}{length q numeric vectors giving the MSE, SS, and Vx and Fcst energy for each threshold (i.e., ignoring the wavelet decomposition).}
  \item{MSEperc,EnVx.perc,EnFcst.perc}{J X q numeric matrices giving percentage form of MSE, Vx Energy and Fcst Energy values, resp.}
  \item{EnRelDiff}{J X q numeric matrix giving the energy relative difference.}
}
\references{
Casati B, G Ross, and DB Stephenson, 2004. A new intensity-scale approach for the verification of spatial precipitation forecasts. \emph{Meteorol. Appl.} \bold{11}, 141--154.

Casati, B., 2010: New Developments of the Intensity-Scale Technique within the Spatial Verification Methods Inter-Comparison Project. \emph{Wea. Forecasting} \bold{25}, (1), 113--143, DOI: 10.1175/2009WAF2222257.1.

Lindsay, R. W., D. B. Percival, and D. A. Rothrock, 1996.  The discrete wavelet transform and the scale analysis of the surface properties of sea ice.
  \emph{IEEE Transactions on Geoscience and Remote Sensing}, \bold{34} (3), 771--787.

Percival, D. B. and P. Guttorp, 1994.  Long-memory processes, the Allan variance and wavelets.  In \emph{Wavelets in Geophysics}, E. Foufoula-Georgiou
 and P. Kumar, Eds., New York: Academic, pp. 325--343.
}
\author{
Eric Gilleland
}

\seealso{
\code{IS}, \code{int.scale.verify} from package \pkg{verification},

\code{\link{dwt.2d}}, \code{\link{modwt.2d}}, \code{\link{idwt.2d}}, \code{\link{imodwt.2d}}, \code{\link{hoods2dPrep}}, \code{\link{hoods2d}}
}
\examples{
data(UKobs6)
data(UKfcst6)
data(UKloc)
hold <- hoods2dPrep( "UKfcst6", "UKobs6", thresholds=c(0.1,0.2,0.5,1,2,5,10,20,50), levels=2^(8-1:8), loc=UKloc, units="mm/h")
look <- waveIS(hold, J=8, verbose=TRUE)
plot(look, which.plots="mse") 
plot(look, which.plots="ss")
plot(look, which.plots="energy")
summary(look)

\dontrun{
data(pert004)
data(pert000)
hold <- hoods2dPrep( "pert004", "pert000", thresholds=c(1,10,50), levels=1:4, loc=ICPg240Locs, units="mm/h")
look <- hoods2d( hold, verbose=TRUE)
plot( look)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }