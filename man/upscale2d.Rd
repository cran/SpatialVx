\name{upscale2d}
\alias{upscale2d}
\alias{plot.upscale2d}
\alias{print.upscale2d}
\title{
Upscaling Neighborhood Verification on a 2-d Verification Set
}
\description{
Perform upscaling neighborhood verification on a 2-d verification set.
}
\usage{
upscale2d(object, time.point = 1, obs = 1,
                 model = 1, levels = NULL, max.n = NULL, smooth.fun =
                 "hoods2dsmooth", smooth.params = NULL, rule = ">=",
                 verbose = FALSE)

\method{plot}{upscale2d}(x, ... )

\method{print}{upscale2d}(x, ...)
}
\arguments{
  \item{object}{
list object of class \dQuote{SpatialVx}.
}

    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
   \item{levels}{numeric vector giving the successive values of the smoothing parameter.  For example, for the default method, these are the neighborhood lengths over which the levels^2 nearest neighbors are averaged for each point.  Values should make sense for the specific smoothing function.  For example, for the default method, these should be odd integers.}
    \item{max.n}{(optional) single numeric giving the maximum neighborhood length to use.  Only used if levels are NULL.}
    \item{smooth.fun}{character giving the name of a smoothing function to be applied.  Default is an average over the n^2 nearest neighbors, where n is taken to be each value of the \code{levels} argument.}
    \item{smooth.params}{list object containing any optional arguments to \code{smooth.fun}.  Use NULL if none.}
   \item{rule}{character string giving the threshold rule to be applied.  See help file for \code{thresholder} function for more information.}
   \item{verbose}{logical, should progress information be printed to the screen?}
   \item{x}{list object of class \dQuote{upscale2d} as returned by \code{upscale2d}.}
   \item{\dots}{optional arguments to the \code{image.plot} function from package \pkg{fields}.  Can also include the argument \code{type}, which must be one of \dQuote{all}, \dQuote{gss}, \dQuote{ts}, \dQuote{bias} or \dQuote{rmse}.}
}
\details{
Upscaling is performed via neighborhood smoothing.  Here, a boxcar kernel is convolved (using the convolution theorem with FFT's) to obtain an average over the nearest n^2 grid squares at each grid point.  This is performed on the raw forecast and verification fields.  The root mean square error (RMSE) is taken for each threshold (Yates et al., 2006; Ebert, 2008).  Further, binary fields are obtained for each smoothed field via thresholding, and frequency bias, threat score ts) and equitable threat score (ets) are calculated (Zepeda-Arce et al., 2000; Ebert, 2008).
}
\value{
\code{upscale2d} returns a list of class \dQuote{upscale2d} with components: 
   \item{rmse }{numeric vector giving the root mean square error for each neighborhood size provided by object.}
   \item{bias,ts,ets }{numeric matrices giving the frequency bias, ts and ets for each neighborhood size (rows) and threshold (columns).}
}
\references{
Ebert, E. E. (2008) Fuzzy verification of high resolution gridded forecasts: A review and proposed framework.  \emph{Meteorol. Appl.}, \bold{15}, 51--64. doi:10.1002/met.25

Yates, E., Anquetin, S., Ducrocq, V., Creutin, J.-D., Ricard, D. and Chancibault, K. (2006) Point and areal validation of forecast precipitation fields.  \emph{Meteorol. Appl.}, \bold{13}, 1--20.

Zepeda-Arce, J., Foufoula-Georgiou, E., Droegemeier, K. K. (2000) Space-time rainfall organization and its role in validating quantitative precipitation forecasts.  \emph{J. Geophys. Res.}, \bold{105}(D8), 10,129--10,146.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{hoods2d}}, \code{\link{kernel2dsmooth}}, \code{\link{kernel2dmeitsjer}}, \code{\link{fft}}
}
\examples{
x <- matrix( 0, 50, 50)
x[ sample(1:50,10), sample(1:50,10)] <- rexp( 100, 0.25)
y <- kernel2dsmooth( x, kernel.type="disk", r=6.5)
x <- kernel2dsmooth( x, kernel.type="gauss", nx=50, ny=50, sigma=3.5)

hold <- make.SpatialVx( x, y, thresholds = seq(0.01,1,,5), field.type = "random")

look <- upscale2d( hold, levels=c(1, 3, 20) )
look

par( mfrow = c(4, 2 ) )
plot( look )

\dontrun{
data( "geom001" )
data( "geom000" )
data( "ICPg240Locs" )

hold <- make.SpatialVx( geom000, geom001, thresholds = c(0.01, 50.01),
    loc = ICPg240Locs, projection = TRUE, map = TRUE, loc.byrow = TRUE,
    field.type = "Precipitation", units = "mm/h",
    data.name = "Geometric", obs.name = "geom000", model.name = "geom001" )

look <- upscale2d(hold, levels=c(1, 3, 9, 17, 33, 65, 129, 257),
    verbose=TRUE)

par( mfrow = c(4, 2 ) )

plot(look )
look <- upscale2d(hold, q.gt.zero=TRUE, verbose=TRUE)
plot(look)
look <- upscale2d(hold, verbose=TRUE)
plot(look)

}
}
\keyword{ math }
