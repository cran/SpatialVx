\name{lossdiff}
\alias{lossdiff}
\alias{lossdiff.default}
\alias{lossdiff.SpatialVx}
\alias{empiricalVG.lossdiff}
\alias{flossdiff}
\alias{summary.lossdiff}
\alias{plot.lossdiff}
\alias{print.lossdiff}
\title{
Test for Equal Predictive Ability on Average Over a Regularly Gridded Space
}
\description{
Test for equal predictive ability (for two forecast models) on average over a regularly gridded space using the method of Hering and Genton (2011).
}
\usage{
lossdiff(x, ...)

\method{lossdiff}{default}(x, ..., xhat1, xhat2, threshold = NULL,
    lossfun = "corrskill", loc = NULL, zero.out = FALSE)

\method{lossdiff}{SpatialVx}(x, ..., time.point = 1, obs = 1, model = c(1, 2),
    threshold = NULL, lossfun = "corrskill", zero.out = FALSE)

empiricalVG.lossdiff( x, trend = 0, maxrad, dx = 1, dy = 1 )

flossdiff(object, vgmodel = "expvg", ...)

\method{summary}{lossdiff}(object, ...)

\method{plot}{lossdiff}(x, ..., icol = c("gray", tim.colors(64)))

\method{print}{lossdiff}(x, ...)
}
\arguments{
  \item{x,xhat1, xhat2}{
\code{lossdiff}: m by n matrices defining the (gridded) verification set where \code{xhat1} and \code{xhat2} are the two forecast models being compared.  \code{plot.lossdiff}: \code{x} is a list returned by \code{lossdiff}.
}
  \item{object}{
 \code{flossdiff} this is the output returned by \code{lossdiff}.  \code{summary.lossdiff}: list object returned by \code{lossdiff} or \code{flossdiff}.
}
  \item{threshold}{numeric vector of length one, two or three giving a threshold under which (non-inclusive) all values will be set to zero.  If length is one, the same threshold is used for all fields (observed, and both models).  If length is two, the same threshold will be used for both models (the second value of \code{threshold}).  Otherwise, the first entry is used for the observed field, the second for the first model and the third for the second model.}
  \item{lossfun}{
character anming a loss function to use in finding the loss differential for the fields.  Default is to use correlation as the loss function.  Must have arguments \code{x} and \code{y}, and may have any additional arguments.
}
  \item{trend}{a matrix (of appropriate dimension) or single numeric (if constant trend) giving the value of the spatial trend.  the value is simply subtracted from the loss differential field before finding the empirical variogram.  If \code{zero.out} is TRUE, then wherever the original three fields all had zero-valued grid points are returned back to zero before continuing (hence ignored in the computation of the variogram).}
  \item{loc}{(optional) mn by 2 matrix giving location coordinates for each grid point.  If NULL, they are taken to be the grid expansion of the dimension of \code{x} (i.e., cbind(rep(1:dim(x)[1],dim(x)[2]), rep(1:dim(x)[2],each=dim(x)[1]))).  This argument is not used by \code{lossdiff}, but may be used subsequently by the \code{plot} method function.}
  \item{maxrad}{
numeric giving the maximum radius for finding variogram differences per the \code{R} argument of \code{vgram.matrix}.
}
  \item{dx, dy}{
\code{dx} and \code{dy} of \code{vgram.matrix}.
}
  \item{zero.out}{logical, should the variogram be computed only over non-zero values of the process?  If TRUE, a modified version of \code{vgram.matrix} is used (\code{variogram.matrix}).}
  \item{vgmodel}{character string naming a variogram model function to use.  Default is the exponential variogram, \code{expvg}. }
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
    \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
    \item{icol}{(optional) color scheme.}
  \item{\dots}{
\code{lossdiff}: optional additional arguments to \code{lossfun}.  Not used by the summary or plot functions.
}
}
\details{
Hering and Genton (2011) introduce a test procedure for comparing spatial fields, which is based on a time series test introduced by Diebold and Mariano (1995).  First, a loss function, g(x,y), is calculated, which can be any appropriate loss function.  This is calculated for each of two forecast fields.  The loss differential field is then given by:

		D(s) = g(x(s),y1(s)) - g(x(s),y2(s)), where s are the spatial locations, x is the verification field, and y1 and y2 are the two forecast fields.

It is assumed that D(s) = phi(s) + psi(s), where phi(s) is the mean trend and psi(s) is a mean zero stationary process with unknown covariance function C(h) = cov(psi(s),psi(s+h)).  In particular, the argument trend represents phi(s), and the default is that the mean is equal (and zero) over the entire domain.  If it is believed that this is not the case, then it should be removed before finding the covariance.

To estimate the trend, see e.g. Hering and Genton (2011) and references therein.

A test is constructed to test the null hypothesis of equal predictive ability on average.  That is,

		H_0: 1/|D| int_D E[D(s)]ds = 0, where |D| is the area of the domain, 

The test statistic is given by

		S_V = mean(D(s))/sqrt(mean(C(h))),

where C(h) = gamma(infinity|p) - gamma(h|p) is a fitted covariance function for the loss differential field.  The test statistic is assumed to be N(0,1) so that if the p-value is smaller than the desired level of significance, the null hypothesis is not accepted.

For 'flossdiff', an exponential variogram is used. Specifically,

		gamma(h | theta=(s,r)) = s^2*(1 - exp(-h/r)),

where s is sqrt(sill) and r is the range (nugget effects are not accounted for here).  If \code{flossdiff} should fail, and the empirical variogram appears to be reasonable (e.g., use the \code{plot} method function on \code{lossdiff} output to check that the empirical variogram is concave), then try giving alternative starting values for the \code{nls} function by using the \code{start.list} argument.  The default is to use the variogram value for the shortest separation distance as an initial estimate for s, and \code{maxrad} as the initial estimate for r.

Currently, it is not possible to fit other variogram models with this function.  Such flexibility may possibly be added in a future release.  In the meantime, use \code{flossdiff} as a template to make your own similar function; just be sure to return an object of class \dQuote{nls}, and it should work seamlessly with the \code{plot} and \code{summary} method functions for a \dQuote{lossdiff} object.  For example, if it is desired to include the nugget or an extra factor (e.g., 3 as used in Hering and Genton, 2011), then a new similar function would need to be created.

Also, although the testing procedure can be applied to irregularly spaced locations (non-gridded), this function is set up only for gridded fields in order to take advantage of computational efficiencies (i.e., use of vgram.matrix), as these are the types of verification sets in mind for this package.  For irregularly spaced grids, the function \code{spct} can be used.

The above test assumes constant spatial trend.  It is possible to remove any spatial trend in D(s) before applying the test.

The procedure requires four steps (hence four functions).  The first is to calculate the loss differential field using \code{lossdiff}.  Next, calculate the empirical variogram of the loss differential field using \code{empiricalVG.lossdiff}.  This second step was originally included within the first step in \code{lossdiff}, but that setup presented a problem for determining if a spatial trend exists or not.  It is important to determine if a trend exists, and if so, to (with care) estimate the trend, and remove it.  If a trend is detected (and estimated), it can be removed before calling \code{empiricalVG.lossdiff} (then use the default \code{trend} = 0), or it can be passed in via the \code{trend} argument; the advantage (or disadvantage) of which is that the trend term will be included in the output object.  The third step is to fit a parametric variogram model to the empirical one using \code{flossdiff}.  The final, fourth step, is to conduct the test, which is performed by the \code{summary} function.

In each step, different aspects of the model assumptions can be checked.  For example, isotropy can be checked by the plot in the lower right panel of the result of the \code{plot} method function after having called \code{empiricalVG.lossdiff}.  The function \code{nlminb} is used to fit the variogram model.

For application to precipitation fields, and introduction to the image warp (coming soon) and distance map loss functions, see Gilleland (2013).

}
\value{
  A list object is returned with possible components:
  \item{data.name}{character vector naming the fields under comparison}
  \item{lossfun,lossfun.args,vgram.args }{same as the arguments input to the lossdiff function.}
  \item{d}{m by n matrix giving the loss differential field, D(s).}
  \item{trend.fit}{An OLS trend fitting the locations to the field via lm.}
  \item{loc}{the self-same value as the argument passed in, or if NULL, it is the expanded grid coordinates.}

empiricalVG.lossdiff returns all of the above (carried over) along with 
  \item{lossdiff.vgram}{list object as returned by vgram.matrix}
  \item{trend}{it is the self-same as the value passed in.}

flossdiff returns all of the above plus:
  \item{vgmodel}{list object as returned by nls containing the fitted exponential variogram model where s is the estimate of sqrt(sill), and r of the range parameter (assuming 'flossdiff' was used to fit the variogram model).}

summary.lossdiff invisibly returns the same list object as above with additional components:
  \item{Dbar}{the estimated mean loss differential (over the entire field).}
  \item{test.statistic}{the test statistic.}
  \item{p.value}{list object with components: two.sided--the two-sided alternative hypothesis--, less--the one-sided alternative hypothesis that the true value mu(D) < 0--and greater--the one-sided alternative hypothesis that mu(D) > 0--, p-values under the assumption of standard normality of the test statistic.}

}

\references{

Diebold, F. X. and Mariano, R. S. (1995) Comparing predictive accuracy. \emph{Journal of Business and Economic Statistics}, \bold{13}, 253--263.

Gilleland, E. (2013) Testing competing precipitation forecasts accurately and efficiently: The spatial prediction comparison test.  \emph{Mon. Wea. Rev.}, \bold{141}, (1), 340--355.

Hering, A. S. and Genton, M. G. (2011) Comparing spatial predictions.  \emph{Technometrics} \bold{53}, (4), 414--425.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{vgram.matrix}}, \code{\link{nls}}, \code{\link{corrskill}}, \code{\link{abserrloss}}, \code{\link{sqerrloss}}, \code{\link{distmaploss}}
}
\examples{
grid<- list( x = seq( 0, 5,, 25), y = seq(0,5,,25) )
obj<-Exp.image.cov( grid = grid, theta = .5, setup = TRUE)

look<- sim.rf( obj )
look[ look < 0 ] <- 0
look <- zapsmall( look )
     
look2 <- sim.rf( obj ) * .25
look2[ look2 < 0 ] <- 0
look2 <- zapsmall( look2 )

look3 <- sim.rf( obj) * 2 + 5
look3[ look3 < 0 ] <- 0 
look3 <- zapsmall( look3 )

res <- lossdiff( x = look, xhat1 = look2, xhat2 = look3, lossfun = "abserrloss" )
res <- empiricalVG.lossdiff( res, maxrad = 8 )
res <- flossdiff( res )
res <- summary( res )

plot( res )

}
\keyword{ htest }
