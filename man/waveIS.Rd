\name{waveIS}
\alias{waveIS}
\alias{waveIS.default}
\alias{waveIS.SpatialVx}
\alias{plot.waveIS}
\alias{summary.waveIS}

\title{
Intensity Scale (IS) Verification
}
\description{
Intensity Scale (IS) verification based on Casat et al (2004) and Casati (2010).
}
\usage{
waveIS(x, th = NULL, J = NULL, wavelet.type = "haar", levels
                 = NULL, max.n = NULL, smooth.fun = "hoods2dsmooth",
                 smooth.params = NULL, rule = ">=", verbose = FALSE,
                 ...)

\method{waveIS}{SpatialVx}(x, th = NULL, J = NULL, wavelet.type = "haar", levels
                 = NULL, max.n = NULL, smooth.fun = "hoods2dsmooth",
                 smooth.params = NULL, rule = ">=", verbose = FALSE,
                 ..., time.point = 1, obs = 1, model = 1 )

\method{waveIS}{default}(x, th = NULL, J = NULL, wavelet.type = "haar", levels
                 = NULL, max.n = NULL, smooth.fun = "hoods2dsmooth",
                 smooth.params = NULL, rule = ">=", verbose = FALSE,
                 ...)

\method{plot}{waveIS}(x, main1 = "X", main2 = "Y",
    which.plots = c("all", "mse", "ss", "energy"),
    level.label = NULL, ...)

\method{summary}{waveIS}(object, ...)
}
\arguments{
   \item{x}{For \code{waveIS} either a list object of class \dQuote{SpatialVx}, a two-component list giving the two fields to be compared (the verification field is assumed to be the first one) or a named list with components \dQuote{X} and \dQuote{Xhat} giving the two fields to be compared.  list object returned by \code{waveIS}.}
    \item{object}{ list object returned by \code{waveIS}.  }
  \item{main1,main2}{character giving labels for the plots where \code{main1} refers to the verification field and \code{main2} to the forecast field.}
  \item{which.plots}{character vector naming one or more specific plots to do.}
  \item{level.label}{optional character vector to use for level names on the plot(s).}
  \item{J}{
numeric integer giving the number of levels to use.  If NULL and the field is dyadic, this will be log2(min(dim(X))), where X is a field from the verification set.  If NULL and the field is not dyadic, then \code{J} is set equal to 4.  Note that if the fields are not dyadic, the function will be much slower.
}
  \item{wavelet.type}{
character giving the name of the wavelet type to use as accepted by \code{dwt.2d} and \code{modwt.2d}.
}
    \item{th}{ list object with named components \dQuote{X} and \dQuote{Xhat} giving the thresholds to use for each field.  If null, taken from teh thresholds attribute for \dQuote{SpatialVx} objects. }
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
   \item{levels}{numeric vector giving the successive values of the smoothing parameter.  For example, for the default method, these are the neighborhood lengths over which the levels^2 nearest neighbors are averaged for each point.  Values should make sense for the specific smoothing function.  For example, for the default method, these should be odd integers.}
    \item{max.n}{(optional) single numeric giving the maximum neighborhood length to use.  Only used if levels are NULL.}
    \item{smooth.fun}{character giving the name of a smoothing function to be applied.  Default is an average over the n^2 nearest neighbors, where n is taken to be each value of the \code{levels} argument.}
    \item{smooth.params}{list object containing any optional arguments to \code{smooth.fun}.  Use NULL if none.}
  \item{rule}{If \code{type} is \dQuote{binary}, return 0 when the rule applied to a grid point's value is not true in relation to the threshold value, and 1 elsewhere.  If \code{type} is \dQuote{replace.below}, then return \code{replace.with} wherever the rule is not true and return the original value otherwise.  By default, it replaces values below the threshold with zero (hence its name), but if \code{rule} is, e.g., \dQuote{<=}, then it will replace values above with zero; or whatever value is chosen for \code{replace.with}.}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{\dots}{Not used by \code{waveIS} (or its method functions) or \code{plot.waveIS}.  Only sort of used by \code{summary.waveIS}.  One can put the argument silent=TRUE so that nothing is printed to the screen (useful if you just want the values calculated and stored without writing to the screen).}
}
\details{
This function applies various statistics to the detail fields (in wavelet space) of a discrete wavelet decomposition (DWT) of the binary error fields for a verification set.  In particular, the statistics described in Casati et al (2004) and Casati (2010) are calculated.  This function depends on the \code{waverify2d} or \code{mowaverify2d} function (depending on whether the fields are dyadic or not, resp.), which themselves depend on the \code{dwt.2d} and \code{idwt.2d} or \code{modwt.2d} and \code{imodwt.2d} functions.

See the references herein and the help files and references therein for \code{dwt.2d} and \code{modwt.2d} for more information on this approach, as well as Percival and Guttorp (1994) and Lindsay et al. (1996).
}
\value{
A list object of class \dQuote{waveIS} that contains the entire prep object passed in by obj, as well as additional components:
  \item{EnVx,EnFcst}{J by q matrices giving the energy for the verification and forecast fields, resp., for each threshold (columns) and scale (rows).}
  \item{MSE,SS }{J by q matrices giving the mean square error and IS skill score for each threshold (column) and scale (rows).}
  \item{Bias}{numeric vector of length q giving the frequency bias of the original fields for each threshold.}
plot.waveIS does not return any value.  A plot is created on the current graphic device.
summary.waveIS returns a list invisibly with the same components as returned by waveIS along with extra components:
  \item{MSEu,SSu,EnVx.u,EnFcst.u}{length q numeric vectors giving the MSE, SS, and Vx and Fcst energy for each threshold (i.e., ignoring the wavelet decomposition).}
  \item{MSEperc,EnVx.perc,EnFcst.perc}{J by q numeric matrices giving percentage form of MSE, Vx Energy and Fcst Energy values, resp.}
  \item{EnRelDiff}{J by q numeric matrix giving the energy relative difference.}
}
\references{
Casati, B., Ross, G. and Stephenson, D. B. (2004) A new intensity-scale approach for the verification of spatial precipitation forecasts. \emph{Meteorol. Appl.} \bold{11}, 141--154.

Casati, B. (2010) New Developments of the Intensity-Scale Technique within the Spatial Verification Methods Inter-Comparison Project. \emph{Wea. Forecasting} \bold{25}, (1), 113--143, doi:10.1175/2009WAF2222257.1.

Lindsay, R. W., Percival, D. B. and Rothrock, D. A. (1996)  The discrete wavelet transform and the scale analysis of the surface properties of sea ice.  \emph{IEEE Transactions on Geoscience and Remote Sensing}, \bold{34} (3), 771--787.

Percival, D. B. and Guttorp, P. (1994)  Long-memory processes, the Allan variance and wavelets.  In \emph{Wavelets in Geophysics}, Foufoula-Georgiou, E. and Kumar, P., Eds., New York: Academic, 325--343.
}
\author{
Eric Gilleland
}

\seealso{
\code{IS}, \code{int.scale.verify} from package \pkg{verification},

\code{\link{dwt.2d}}, \code{\link{modwt.2d}}, \code{\link{idwt.2d}}, \code{\link{imodwt.2d}}, \code{\link{hoods2d}}

\code{\link{thresholder}}
}
\examples{
data( "UKobs6" )
data( "UKfcst6" )
data( "UKloc" )

hold <- make.SpatialVx( UKobs6, UKfcst6,
    thresholds = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
    loc = UKloc, map = TRUE, field.type = "Rainfall", units = "mm/h",
    data.name = "Nimrod", obs.name = "UKobs6", model.name = "UKfcst6" )

look <- waveIS(hold, J=8, levels=2^(8-1:8), verbose=TRUE)
plot(look, which.plots="mse") 
plot(look, which.plots="ss")
plot(look, which.plots="energy")
summary(look)

}
\keyword{ math }
