\name{upscale2d}
\alias{upscale2d}
\alias{plot.upscale2d}
\title{
Perform upscaling neighborhood verification on a 2-d verification set
}
\description{
Perform upscaling neighborhood veriication on a 2-d verification set.
}
\usage{
upscale2d(object, thresholds = NULL, verbose = FALSE)
\method{plot}{upscale2d}(x, ...)
}
\arguments{
  \item{object}{
list object output from the \code{hoods2dPrep} function.
}
  \item{thresholds}{
(optional) alternative thresholds to use on the upscaled fields.  If NULL, the thresholds from the \code{hoods2dPrep} object are used.
}
   \item{verbose}{logical, should progress information be printed to the screen?}
   \item{x}{list object of class "upscale2d" as returned by \code{upscale2d}.}
   \item{\dots}{optional arguments to the \code{image.plot} function from package \pkg{fields}.}
}
\details{
Upscaling is performed via neighborhood smoothing.  Here, a boxcar kernel is convolved (using the convolution theorem with FFT's) to obtain an average over the nearest n^2 grid squares at each grid point.  This is performed on the raw forecast and verification fields.  The root mean square error (RMSE) is taken for each threshold (Yates et al., 2006; Ebert, 2008).  Further, binary fields are obtained for each smoothed field via thresholding, and frequency bias, threat score ts) and equitable threat score (ets) are calculated (Zepeda-Arce et al., 2000; Ebert, 2008).
}
\value{
\code{upscale2d} returns a list of class "upscale2d" with components: 
   \item{rmse }{numeric vector giving the root mean square error for each neighborhood size provided by object.}
   \item{bias,ts,ets }{numeric matrices giving the frequency bias, ts and ets for each neighborhood size (rows) and threshold (columns).}
}
\references{
Ebert, E. E., 2008: Fuzzy verification of high resolution gridded forecasts: A review and proposed framework.  \emph{Meteorol. Appl.}, \bold{15}, 51--64. DOI: 10.1002/met.25 Available at \url{http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf}

Yates, E., S. Anquetin, V. Ducrocq, J.-D. Creutin, D. Ricard, and K. Chancibault, 2006: Point and areal validation of forecast precipitation fields.  \emph{Meteorol. Appl.}, \bold{13}, 1--20.

Zepeda-Arce, J., E. Foufoula-Georgiou, K. K. Droegemeier, 2000: Space-time rainfall organization and its role in validating quantitative precipitation forecasts.  \emph{J. Geophys. Res.}, \bold{105}(D8), 10,129--10,146.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{hoods2d}}, \code{\link{hoods2dPrep}}, \code{\link{kernel2dsmooth}}, \code{\link{kernel2dmeitsjer}}, \code{\link{fft}}
}
\examples{
x <- matrix( 0, 50, 50)
x[ sample(1:50,10), sample(1:50,10)] <- rexp( 100, 0.25)
y <- kernel2dsmooth( x, kernel.type="disk", r=6.5)
x <- kernel2dsmooth( x, kernel.type="gauss", nx=50, ny=50, sigma=3.5)
hold <- hoods2dPrep("y", "x", thresholds=c(1e-8, 0.05, 0.1), levels=c(1, 3, 20))
look <- upscale2d( hold, thresholds=seq(0.01,1,,5))
look
plot( look)

\dontrun{plot( look)}

\dontrun{
data(geom001)
data(geom000)
data(ICPg240Locs)
hold <- hoods2dPrep( "geom001", "geom000", thresholds=c(0.01,50.01),
	levels=c(1, 3, 9, 17, 33, 65, 129, 257), loc=ICPg240Locs, units="in/100")
look <- upscale2d(hold, thresholds=c(1e-8, 0.1, 0.5, 45, 50.1), verbose=TRUE)
plot( look) # Might want to use 'pdf' to print these out so that all of them can be observed.
data(pert004)
data(pert000)
hold <- hoods2dPrep( "pert004", "pert000", thresholds=c(1,10,50),
	levels=c(1, 3, 17, 33, 65, 129, 257), loc=ICPg240Locs, units="mm/h")
look <- upscale2d( hold, thresholds=c(1e-8, 3:5, 45, 50.1), verbose=TRUE)
plot( look)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
