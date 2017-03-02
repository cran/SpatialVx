\name{variographier}
\alias{variographier}
\alias{variographier.default}
\alias{variographier.SpatialVx}
\title{
Variography Score
}
\description{
Calculate the variography score between two spatial fields based on the fitted exponential variogram.
}
\usage{
variographier(x, init, zero.out = FALSE, ...)

\method{variographier}{default}( x, init, zero.out = FALSE, ..., y )

\method{variographier}{SpatialVx}( x, init, zero.out = FALSE, ...,
    obs = 1, model = 1, time.point = 1 )
}
\arguments{
  \item{x, y}{
matrices giving the fields on which to calculate the variography or a \dQuote{SpatialVx} class object (\code{x} only).
}
  \item{init}{
list with components \code{px} and \code{py} that give initial values for parameter estimates (sill + nugget and range).  If missing, default will attempt to find reasonable starting values.
}
  \item{zero.out}{
logical should the variogram be calculated over all grid points or just ones where one or both fields are non-zero?  See \code{variogram.matrix}.
}
\item{time.point}{
numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.
}
  \item{obs, model}{
numeric indicating which observation/forecast model to select for the analysis.
}

  \item{\dots}{
optional arguments to \code{vgram.matrix} or \code{variogram.matrix} (if zero.out is TRUE).  Can also have optional arguments to \code{nlminb} (but not \code{lower} or \code{upper}).
}
}
\details{
The variography score calculated here is that from Ekstr\"{o}m (2016).  So far, only the exponential variogram is allowed.

Note that in the fitting, the model g(h) = c * ( 1 - exp( -a * h ) ) is used, but the variography is calculated for theta = 3 / a.  Therefore, the values in the par component of the returned fitted variograms correspond to a, while the variography score corresponds to theta.  The score is given by:

v = 1 / sqrt( c_0^2 + c_m^2 + ( theta_0 - theta_m )^2 )

where c_0 and c_m are the sill + nugget terms for the observation and model, resp., and similarly for theta_0 and theta_m.

The parameters are *not* currently normalized, here, to give equal weight between sill + nugget and range.  If several fields are analyzed (e.g., an ensemble), then the fitted parameters could be gathered, and one could use that information to calculate the score based on a normalized version.
}
\value{
A list object of class \dQuote{variographied} is returned with components:
  \item{obs.vg, mod.vg }{Empirical variogram objects as returned by either vgram.matrix or variogram.matrix}
  \item{obs.parvg, mod.parvg }{objects returned by nlminb containing the fitted exponential variogram model parameters and some information about the optimization.  }

   \item{variography}{ single numeric giving the variography measure.}
}
\references{
Ekstr\"{o}m, M. (2016) Metrics to identify meaningful downscaling skill in WRF simulations of intense rainfall events.  \emph{Environmental Modelling and Software}, \bold{79}, 267--284, DOI: 10.1016/j.envsoft.2016.01.012.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{vgram.matrix}}, \code{\link{variogram.matrix}}
}
\examples{
data( "UKobs6" )
data( "UKfcst6" )
data( "UKloc" )

hold <- make.SpatialVx( UKobs6, UKfcst6, thresholds = c(0.01, 20.01),
    loc = UKloc, field.type = "Precipitation", units = "mm/h",
    data.name = "Nimrod", obs.name = "Observations 6", model.name = "Forecast 6",
    map = TRUE)

look <- variographier( hold )
look
plot( look )
}
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
