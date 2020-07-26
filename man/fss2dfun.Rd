\name{fss2dfun}
\alias{fss2dfun}
\alias{fuzzyjoint2dfun}
\alias{MinCvg2dfun}
\alias{multicon2dfun}
\alias{pragmatic2dfun}
\alias{upscale2dfun}
\title{Various Verification Statistics on Possibly Neighborhood-Smoothed Fields.}
\description{
Functions to calculate various verification statistics on possibly neighborhood smoothed fields.  Used by hoods2d, but can be called on their own.}
\usage{
fss2dfun(sPy, sPx, subset = NULL, verbose = FALSE)

fuzzyjoint2dfun(sPy, sPx, subset = NULL)

MinCvg2dfun(sIy, sIx, subset = NULL)

multicon2dfun(sIy, Ix, subset = NULL)

pragmatic2dfun(sPy, Ix, mIx = NULL, subset = NULL)

upscale2dfun(sYy, sYx, threshold = NULL, which.stats = c("rmse",
                 "bias", "ts", "ets"), rule = ">=", subset = NULL)

}
\arguments{
  \item{sPy}{n by m matrix giving a smoothed binary forecast field.}
  \item{sPx}{n by m matrix giving a smoothed binary observed field.}
  \item{sIy}{n by m matrix giving a binary forecast field.}
  \item{sIx}{n by m matrix giving a binary observed field (the s indicates that the binary field is obtained from a smoothed field).}
  \item{Ix}{n by m matrix giving a binary observed field.}
  \item{mIx}{(optional) single numeric giving the base rate.  If NULL, this will be calculated by the function.  Simply a computation saving step if this has already been calculated.}
  \item{sYy}{n by m matrix giving a smoothed forecast field.}
  \item{sYx}{n by m matrix giving a smoothed observed field.}
  \item{threshold}{(optional) numeric vector of length 2 giving the threshold over which to calculate the verification statistics: bias, ts and ets.  If NULL, only the rmse will be calculated.}
  \item{which.stats}{character vector naming which statistic(s) should be caluclated for \code{upscale2dfun}.}
  \item{subset}{(optional) numeric indicating over which points the summary scores should be calculated.  If NULL, all of the points are used.}
  \item{rule}{character string giving the sort of thresholding process desired.  See the help file for \code{thresholder} for more information.}
  \item{verbose}{logical, should progress information be printed to the screen?}
}
\details{
These are modular functions that calculate the neighborhood smoothing method statistics in spatial forecast verification (see, e.g., Ebert, 2008, 2009; Gilleland et al., 2009, 2010; Roberts and Lean,2008).  These functions take fields that have already had the neighborhood smoothing applied (e.g., using \code{kernele2d}) when appropriate.  They are called by \code{hoods2d}, so need not be called by the user, but they can be.
}

\value{In the case of \code{fss2dfun}, a single numeric giving the FSS value is returned.  In the other cases, list objects are returned with one or more of the following components, depending on the particular function.
  \item{fuzzy}{\code{fuzzyjoint2dfun} returns a list with this list as one component.  The list component fuzzy has the components: pod, far and ets.}
  \item{joint}{\code{fuzzyjoint2dfun} returns a list with this list as one component.  The list component joint has the components: pod, far and ets.}
  \item{pod}{numeric giving the probability of detection, or hit rate.}
  \item{far}{numeric giving the false alarm ratio.}
  \item{ets}{numeric giving the equitable threat score, or Gilbert Skill Score.}
  \item{f}{numeric giving the false alarm rate.}
  \item{hk}{numeric giving the Hanssen-Kuipers statistic.}
  \item{bs}{Brier Score}
  \item{bss}{Brier Skill Score.  The \code{pragmatic2dfun} returns the bs and bss values.  The Brier Skill Score here uses the mean square error between the base rate and the Ix field as the reference forecast.}
  \item{ts}{numeric giving the threat score.}
  \item{bias}{numeric giving the frequency bias.}
 }

\references{
Ebert, E. E. (2008) Fuzzy verification of high resolution gridded forecasts: A review and proposed framework.  \emph{Meteorol. Appl.}, \bold{15}, 51--64. doi:10.1002/met.25 

Ebert, E. E. (2009) Neighborhood verification: A strategy for rewarding close forecasts.  \emph{Wea. Forecasting}, \bold{24}, 1498--1510, doi:10.1175/2009WAF2222251.1.

Gilleland, E., Ahijevych, D., Brown, B. G., Casati, B. and Ebert, E. E. (2009) Intercomparison of Spatial Forecast Verification Methods. \emph{Wea. Forecasting}, \bold{24}, 1416--1430, doi:10.1175/2009WAF2222269.1.

Gilleland, E., Ahijevych, D. A., Brown, B. G. and Ebert, E. E. (2010) Verifying Forecasts Spatially. \emph{Bull. Amer. Meteor. Soc.}, October, 1365--1373.

Roberts, N. M. and Lean, H. W. (2008) Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events.  \emph{Mon. Wea. Rev.}, \bold{136}, 78--97. doi:10.1175/2007MWR2123.1.}
\author{Eric Gilleland}
\seealso{\code{\link{hoods2d}},\code{\link{kernel2dsmooth}},\code{\link{vxstats}}, \code{\link{thresholder}} }
\examples{
x <- y <- matrix( 0, 100, 100)
x[ sample(1:100, 10), sample(1:100, 10)] <- 1
y[ sample(1:100, 20), sample(1:100, 20)] <- 1
Px <- kernel2dsmooth( x, kernel.type="boxcar", n=9, xdim=c(100, 100))
Py <- kernel2dsmooth( y, kernel.type="boxcar", n=9, xdim=c(100, 100))
par( mfrow=c(2,2))
image( x, col=c("grey", "darkblue"), main="Simulated Observed Events")
image( y, col=c("grey", "darkblue"), main="Simulated Forecast Events")
image( Px, col=c("grey", tim.colors(256)), main="Forecast Event Frequencies (9 nearest neighbors)")
image( Py, col=c("grey", tim.colors(256)), main="Smoothed Observed Events (9 nearest neighbors)")
fss2dfun( Py, Px)

}
\keyword{ math }
