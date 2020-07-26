\name{wavePurifyVx}
\alias{wavePurifyVx}
\alias{wavePurifyVx.default}
\alias{wavePurifyVx.SpatialVx}
\alias{plot.wavePurifyVx}
\alias{summary.wavePurifyVx}
\title{
Apply Traditional Forecast Verification After Wavelet Denoising
}
\description{
Apply traditional forecast verification after wavelet denoising ala Briggs and Levine (1997).
}
\usage{
wavePurifyVx( x, climate = NULL, which.stats = c("bias",
    "ts", "ets", "pod", "far", "f", "hk", "mse"), thresholds = NULL,
    rule = ">=", return.fields = FALSE, verbose = FALSE, ...)

\method{wavePurifyVx}{SpatialVx}( x, climate = NULL, which.stats = c("bias",
    "ts", "ets", "pod", "far", "f", "hk", "mse"), thresholds = NULL,
    rule = ">=", return.fields = FALSE, verbose = FALSE, ...,
    time.point = 1, obs = 1, model = 1 )

\method{wavePurifyVx}{default}( x, climate = NULL, which.stats = c("bias",
    "ts", "ets", "pod", "far", "f", "hk", "mse"), thresholds = NULL,
    rule = ">=", return.fields = FALSE, verbose = FALSE, ...)

\method{plot}{wavePurifyVx}(x, ..., col = c("gray", tim.colors(64)), zlim, mfrow,
                 horizontal = TRUE, type = c("stats", "fields") )

\method{summary}{wavePurifyVx}(object, ...)
}
\arguments{
  \item{x}{
For \code{wavePurifyVx}, either a list object of class \dQuote{SpatialVx}, or a list with only two components consisting of m by n matrices giving the verification and forecast fields, resp., or a list with named components \dQuote{X} and \dQuote{Xhat}.

For \code{plot.wavePurifyVx}, list object as output from \code{wavePurifyVx}.

}
  \item{object}{ list object as returned by \code{wavePurifyVx}.  }
  \item{climate}{m by n matrix defining a climatology field.  If not NULL, then the anamoly correlation coefficient will be applied to the wavelet denoised fields.}
  \item{which.stats}{
character describing which traditional verification statistics to calculate on the wavelet denoised fields.  This is the argument passed to the argument of the same name in \code{vxstats}.
}
  \item{thresholds}{
Either a numeric vector or a list with components named \dQuote{X} and \dQuote{Xhat} giving thresholds  used to define events for all of the verification statistics except MSE.  However, if supplied or other statistics are to be computed, then MSE will be calculated for the fields at values >= \code{thresholds}.  If only MSE is to be computed, and \code{thresholds} is NULL, then no thresholding is applied.  If NULL, and any of the statistics besides MSE are to be calculated, then default values of the 0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9 and 0.95 quantiles (for each field, so that the thresholds may differ between fields) are used.  The same holds for anamoly correlation coefficient.  The exception is that if the argument is null and \code{x} is a \dQuote{SpatialVx} class object, then the thresholds are taken to be those associated with the \dQuote{thresholds} attribute of this object.
}
  \item{rule}{character string giving the threshold rule to be applied.  See help file for \code{thresholder} function for more information.}
  \item{return.fields}{logical, should the denoised fields be returned (e.g., for subsequent plotting)?}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
  \item{verbose}{logical, should progress information (including total run time) be printed to the screen?}
   \item{col, zlim, horizontal}{optional arguments to \code{image}, and/or \pkg{fields} functions \code{poly.image} and \code{image.plot}}
   \item{mfrow}{ optionally set the plotting panel via \code{mfrow} (see help file for \code{par}).  Default sets a rgion that will show all plots in one set of panels. }
   \item{type}{character string stating whether to plot the resulting statistics or the original fields along with their de-noised counter parts.}
  \item{\dots}{
For \code{wavePurifyVx}, optional additional arguments to \code{denoise.dwt.2d} (or \code{denoise.modwt.2d}) from package \pkg{waveslim}.  Note that if the argument \code{J} is not passed, then it will be determined as J=log2(min(m,n)).  If the fields are dyadic, then the usual DWT is used, otherwise the maximal overlap DWT is used instead.  For the plot and summary method functions, these are not used.  Also passed to \code{poly.image} and \code{image.plot} for plotting routine for \dQuote{fields} type plots when the \dQuote{maps} attribute from the \dQuote{SpatialVx} object is TRUE.
}
}
\details{
If the fields are dyadic, then the \code{denoise.dwt.2d} function from package \pkg{waveslim} is applied to each field before calculating the chosen verification statistics.  Otherwise \code{denoise.modwt.2d} from the same package is used.  The result is that high-frequency fluctuations in the two fields are removed before calculating verification statistics so that the resulting statistics are less susceptible to small-scale errors (see Briggs and Levine, 1997).  See Percival and Guttorp (1994) and Lindsay et al. (1996) for more on this type of wavelet analysis including maximal overlap DWT.
}
\value{
A list object of class \dQuote{wavePurifyVx} is returned with possible components (depending on what is supplied in the arguments, etc.):
  \item{X2,Y2}{m by n matrices of the denoised verification and forecast fields, resp. (only if return.fields is TRUE).}
  \item{thresholds }{q by 2 matrix of thresholds applied to the forecast (first column) and verification (second column) fields, resp.  If climate is not NULL, then the same thresholds for the forecast field are applied to the climatology.}
  \item{qs }{If object and thresholds are NULL, and statistics other than MSE or ACC are desired, then this will be created along with the thresholds, and is just a character version of the trhesholds. }
  \item{args}{list object containing all the optional arguments passed into ..., and the value of J used (e.g., even if not passed into ...).}
  \item{bias,ts,ets,pod,far,f,hk,mse,acc}{numeric vectors of length q (i.e., the number of thresholds) giving the associated verification statistics.}
}
\references{
Briggs, W. M. and Levine, R. A. (1997) Wavelets and field forecast verification. \emph{Mon. Wea. Rev.}, \bold{125}, 1329--1341.

Lindsay, R. W., Percival, D. B. and Rothrock, D. A. (1996)  The discrete wavelet transform and the scale analysis of the surface properties of sea ice.  \emph{IEEE Transactions on Geoscience and Remote Sensing}, \bold{34} (3), 771--787.

Percival, D. B. and Guttorp, P. (1994)  Long-memory processes, the Allan variance and wavelets.  In \emph{Wavelets in Geophysics}, Foufoula-Georgiou, E. and Kumar, P., Eds., New York: Academic, pp. 325--343.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{denoise.dwt.2d}}, \code{\link{denoise.modwt.2d}}, \code{\link{waverify2d}}, \code{\link{mowaverify2d}}, \code{\link{waveIS}}, \code{\link{vxstats}}
}
\examples{
grid <- list( x= seq( 0,5,,100), y= seq(0,5,,100))
obj <- Exp.image.cov( grid=grid, theta=.5, setup=TRUE)
look <- sim.rf( obj)
look[ look < 0] <- 0
look <- zapsmall( look)
     
look2 <- sim.rf( obj)
look2[ look2 < 0] <- 0
look2 <- zapsmall( look2)

look3 <- sim.rf( obj) 
look3[ look3 < 0] <- 0 
look3 <- zapsmall( look3)

hold <- make.SpatialVx( look, look2, thresholds = c(0.1, 1),
            field.type = "random", units = "units")

plot( hold )

res <- wavePurifyVx( hold, climate = look3, return.fields = TRUE, verbose = TRUE )

plot(res, type="fields")
plot(res, type="stats")

summary(res)

\dontrun{
data( "UKobs6" )
data( "UKfcst6" )
data( "UKloc" )

hold <- surrogater2d( UKobs6, n=1, maxiter=50, verbose=TRUE)
hold <- matrix(hold, 256, 256)

UKobj <- make.SpatialVx( UKobs6, UKfcst6, thresholds = c(0.1, 2, 5, 10),
    loc = UKloc, map = TRUE, field.type = "Rainfall", units = "mm/h",
    data.name = "Nimrod", obs.name = "obs 6", model.name = "fcst 6" )

plot(UKobj ) 

look <- wavePurifyVx( object = UKobj, climate = hold,
    return.fields = TRUE, verbose = TRUE)

plot(look, type = "fields" )
plot(look, type = "stats" )

summary( look )

}

}
\keyword{ math }
