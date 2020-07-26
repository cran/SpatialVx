\name{pphindcast2d}
\alias{pphindcast2d}
\alias{plot.pphindcast2d}
\alias{print.pphindcast2d}
\title{
Practically Perfect Hindcast Neighborhood Verification Method
}
\description{
Function to perform the practically perfect hindcast neighborhood verification method.  Finds the
optimal threhsold, Pthresh, and calculates the desired statistic for that threshold.
}
\usage{
pphindcast2d(object, which.score = "ets", time.point = 1, obs = 1,
                 model = 1, levels = NULL, max.n = NULL, smooth.fun =
                 "hoods2dsmooth", smooth.params = NULL, rule = ">=",
                 verbose = FALSE, ...)

\method{plot}{pphindcast2d}(x, ..., mfrow = NULL,
    type = c("quilt", "line"), 
    col = heat.colors(12), horizontal = FALSE)

\method{print}{pphindcast2d}(x, ...)
}
\arguments{
  \item{object}{
A list object returned by the \code{make.SpatialVx} function.
}
  \item{which.score}{
character stating which verification score is to be used.  Must be one that is accepted by \code{vxstats}.
}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
   \item{levels}{numeric vector giving the successive values of the smoothing parameter.  For example, for the default method, these are the neighborhood lengths over which the levels^2 nearest neighbors are averaged for each point.  Values should make sense for the specific smoothing function.  For example, for the default method, these should be odd integers.}
    \item{max.n}{(optional) single numeric giving the maximum neighborhood length to use.  Only used if levels are NULL.}
    \item{smooth.fun}{character giving the name of a smoothing function to be applied.  Default is an average over the n^2 nearest neighbors, where n is taken to be each value of the \code{levels} argument.} 
    \item{smooth.params}{list object containing any optional arguments to \code{smooth.fun}.  Use NULL if none.}
   \item{rule}{character string giving the threshold rule to be applied.  See help file for \code{thresholder} function for more information.}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{x}{An object of class \dQuote{pphindcast2d} as returned by the self-same function.}
    \item{mfrow}{ mfrow parameter (see help file for \code{par}).  If NULL, then the parameter is not re-set. }
  \item{type}{character specifying whether two quilt plots (one for the score and one for Pthresh) should be made, or one line plot incorporating both the score and the Pthresh values; the latter's values being displayed on the right axis.}
  \item{col, horizontal}{arguments used in the calls by \code{image} and \code{image.plot}.}
  \item{\dots}{
\code{pphindcast2d}: optional arguments to the \code{optim} function.  May not include lower, upper or method as these are hard coded into the function.

\code{plot} method function: optional arguments to the \code{image} function.

\code{print} method function: not used.
}
}
\details{
The practically perfect hindcast method is described in Ebert (2008).  Using a similar notation as that described therein (and in the help
page for \code{hoods2d}), the method is a SO-NF approach that first compares the observed binary field (obtained from the trheshold(s) provided by \code{object}), Ix, with the smoothed binary field, <Px>s.  This smoothed binary field is thresholded by
Pthresh to obtain a new binary field.  The value of Pthresh that maximizes the verification score (provided by the which.score argument)
is then used to compare Ix with <Iy>s, the binary forecast field obtained by thresholding the smoothed binary forecast field Iy using
the value of Pthresh found above.  The verification statistic determined by which.score is calculated between Ix and <Iy>s.
}
\value{
 A list object is returned with components:
    \item{which.score }{value of which.score, same as the argument passed in.}
    \item{Pthresh }{l by q matrix giving the value of Pthresh applied at each level (rows) and threshold (columns).}
    \item{values}{l by q matrix giving the value of which.score found for each level (rows) and threshold (columns).}
}
\references{
Ebert, E. E. (2008) Fuzzy verification of high resolution gridded forecasts: A review and proposed framework.  \emph{Meteorol. Appl.}, \bold{15}, 51--64. doi:10.1002/met.25 

}
\author{
Eric Gilleland
}

\section{Warning }{The value Pthresh is optimized under the assumption that larger values of which.score are better.} 

\seealso{
\code{\link{hoods2d}}, \code{\link{kernel2dsmooth}}, \code{\link{vxstats}}, \code{\link{hoods2dPlot}}, \code{\link{optim}}
}
\examples{
x <- y <- matrix( 0, 50, 50)
x[ sample(1:50,10), sample(1:50,10)] <- rexp( 100, 0.25)
y[ sample(1:50,20), sample(1:50,20)] <- rexp( 400)

hold <- make.SpatialVx( x, y, thresholds=c(0.1, 0.5), field.type = "random")
look <- pphindcast2d(hold, levels=c(1, 3))
look
\dontrun{
data( "geom001" )
data( "geom000" )
data( "ICPg240Locs" )

hold <- make.SpatialVx( geom000, geom001, thresholds = c(0.01, 50.01),
    loc = ICPg240Locs, projection = TRUE, map = TRUE, loc.byrow = TRUE,
    data.name = "Geometric", obs.name = "geom000", model.name = "geom001",
    field.type = "Precipitation", units = "mm/h")

look <- pphindcast2d( hold, levels=c(1, 3, 65), verbose=TRUE)

plot(look, mfrow = c(1, 2) )
plot(look, mfrow = c(1, 2), type = "line")

# Alternatively:
par( mfrow = c(1, 2) )
hoods2dPlot( look$values, args = attributes( look ),
    main="Gilbert Skill Score")
}
}
\keyword{ math }
