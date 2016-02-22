\name{spct}
\alias{spct}
\title{
Spatial Prediction Comparison Test
}
\description{
Spatial Prediction Comparison Test (SPCT) for spatial locations that are on a regular or irregular coordinate system.
}
\usage{
spct(d, loc, trend = 0, lon.lat = TRUE,
    dmax = NULL, vgmodel = "expvgram", vgmodel.args = NULL,
    init, alpha = 0.05, alternative = c("two.sided", "less", "greater"), mu = 0,
    verbose = FALSE, ...)
}
\arguments{
  \item{d}{
numeric vector of length n giving the (spatial) loss differential field (at a single point in time).
}
  \item{loc}{
n by 2 numeric matrix giving the spatial coordinates for each data point in \code{d}.
}
  \item{trend}{
a numeric vector of length one or n to be subtracted from d before finding the variogram and perfomring the test.
}
  \item{lon.lat}{
logical stating whether or not the values in \code{loc} are longitude/latitude coordinates or not.  If TRUE, then the \pkg{fields} function \code{rdist.earth} is used to calculate distances.  If FALSE, then the \pkg{fields} function \code{rdist} is used.
}
  \item{dmax}{
single numeric giving the maximum lag distance over which to fit the parametric variogram model.  The default uses half of the maximum lag.
}
  \item{vgmodel}{
character string naming a function defining the parametric variogram model to be used.  The default uses \code{expvgram}, the exponential variogram model.  Must have arguments \code{p} (vector of parameters), \code{h} (vector of distances) and \code{\dots}.
}
  \item{vgmodel.args}{
Optional list of other arguments to be passed to \code{vgmodel}.  Not used by the default method.
}
  \item{init}{
Initial parameter values to be used in the call to \code{nlminb} for estimating the parameters of the variogram model.  The default for the default exponential variogram is to use the square root of the first-lag of the empirical variogram for the nugget and the difference between the second and first lag variogram values (if the second lag term is positive), and the first lag term otherwise for the range parameter.
}
  \item{alpha}{
single numeric giving the desired level of significance.
}
  \item{alternative}{
character string naming which type of hypothesis test to conduct.  Default is to do a two-sided test.  Note that the SPCT is paired test.
}
  \item{mu}{
The mean loss differential value under the null hypothesis.  Usually, this will be zero (the default value).
}
  \item{verbose}{
logical, should progress information be printed to the screen?  It may also provide other useful information in the event that a problem occurs somewhere.
}
  \item{\dots}{
Optional arguments to \code{vgram} from the \pkg{fields} package.
}
}
\details{
If using a large spatial data set that occurs on a regular grid, you should probably use \code{lossdiff}, \code{empiricalVG.lossdiff}, \code{flossdiff} and \code{summary} to perform this self-same test (the SPCT), as those functions make use of special tricks for regular grids to speed things up.  Otherwise, this function should work on either type of grid.

The SPCT is a paired test introduced by Hering and Genton (2011)--and based on the time series test introduced by Diebold and Mariano (1995) of whether one of two competing forecasts is better than the other (alternative) or not (null).  Apart from being a test for spatial fields, the SPCT test fits a parametric model to the empirical variogram (instead of using the empirical one), which turns out to be more accurate.

The loss differential field is a field giving the straight difference between the two loss functions calculated for each of two forecasts.  For example, suppose Z(x,y) is an observed spatial field with (possibly irregularly spaced) locations (x, y), and Y1(x, y) and Y2(x, y) are two competing forecasts.  One might be interested in whether or not, on average, the difference in the absolute error for Y1 and Y2 is significantly different from zero.  First, g1 = abs( Y1(x, y) - Z(x, y) ) and g2 = abs( Y2(x, y) - Z(x, y) ).  Second, the loss differential field is D(x, y) = g1 - g2.  It is the average of D(x, y) that is of interest.  Because D(x, y) is likely to have a strong spatial correlation, the standard error for Dbar = mean( D(x, y) ) is calculated from the variogram.  Hering and Genton (2011) found the test to have proper size and good power, and found it to be relatively robust to contemporaneous correlation--i.e., if Y1 and Y2 are correlated (even if they are not, which is unlikely, g1 and g2 will necessarily be correlated because both involve the same field Z).

If the sample size is less than 30, a t-test is used, and a normal approximation otherwise.

See also, Gilleland (2013) for a modification of this test that accounts for location errors (coming soon).
}
\value{
A list object of class \dQuote{htest} with components:
  \item{data.name }{a character string giving the name of the loss differential field.}
  \item{loss.differential }{The original loss differential field as passed by argument d.}
  \item{nloc}{the number of spatial locations.}
  \item{trend}{Same as the argument passed in.}
  \item{optional.arguments}{list with any arguments passed into vgram.}
  \item{empirical.variogram}{the object returned by vgram giving the empirical variogram.}
  \item{parametric.vgram.fit}{the value returned by nlminb or an object of class \dQuote{try-error}.}
  \item{estimate}{the estimated mean loss differential value.}
  \item{se}{the estimated standard error estimated from the fitted variogram model.}
  \item{statistic}{the value of the statistic ( mean( d ) - mu ) / se. }
  \item{null.value}{the argument mu.}
  \item{parameter}{numeric vector giving the parameter values estimated for the variogram model.}
  \item{fitted.values}{the predicted variogram values from the fitted parametric model.}
  \item{loss.differential.detrended}{this is the loss differential field after having been de-trended. }
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{p.value}{the p-value for the test.}
  \item{conf.int}{The (1 - alpha) * 100 percent confidence interval found using the standard error based on the variogram model per hering and Genton (2011).}
  \item{method}{a character string indicating the type of test performed.}
}
\references{
Diebold, F.X. and Mariano, R.S. (1995) Comparing predictive accuracy. \emph{Journal of Business and Economic Statistics}, \bold{13}, 253--263.

Gilleland, E. (2013) Testing competing precipitation forecasts accurately and efficiently: The spatial prediction comparison test.  \emph{Mon. Wea. Rev.}, \bold{141}, (1), 340--355.

Hering, A. S. and Genton, M. G. (2011) Comparing spatial predictions.  \emph{Technometrics} \bold{53}, (4), 414--425.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{vgram}}, \code{\link{lossdiff}}, \code{\link{flossdiff}}, \code{\link{summary.lossdiff}}, \code{\link{expvgram}}
}
\examples{

\dontrun{
y1 <- predict( Tps( fields::ozone$x, fields::ozone$y ) )
y2 <- predict( Krig( fields::ozone$x, fields::ozone$y, theta = 20 ) )

y <- fields::ozone$y

spct( abs( y1 - y ) - abs( y2 - y ), loc = fields::ozone$x )

spct( abs( y1 - y ) - abs( runif( 20, 1, 5 ) - y ), loc = fields::ozone$x )
}
}
\keyword{ htest }
\keyword{ models }
