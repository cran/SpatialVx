\name{spatMLD}
\alias{spatMLD}
\alias{spatMLD.default}
\alias{spatMLD.SpatialVx}
\alias{fit.spatMLD}
\alias{summary.spatMLD}
\alias{plot.spatMLD}
\alias{print.spatMLD}
\title{
Test for Equal Predictive Ability on Average Over a Regularly Gridded Space
}
\description{
Test for equal predictive ability (for two forecast models) on average over a regularly gridded space using the method of Hering and Genton (2011).
}
\usage{
spatMLD(x, ...)

\method{spatMLD}{default}(x, ..., xhat1, xhat2, lossfun = "corrskill",
    trend = "ols", loc = NULL, maxrad = 20, dx = 1, dy = 1,
    zero.out = FALSE)

\method{spatMLD}{SpatialVx}(x, ..., time.point = 1, model = c(1, 2), lossfun =
                 "corrskill", trend = "ols", maxrad, dx = 1, dy = 1,
                 zero.out = FALSE)

fit.spatMLD(object, vgmodel = "expvg", ...)

\method{summary}{spatMLD}(object, ...)

\method{plot}{spatMLD}(x, ..., icol = c("gray", tim.colors(64)))

\method{print}{spatMLD}(x, ...)
}
\arguments{
  \item{x,xhat1, xhat2}{
\code{spatMLD}: m by n matrices defining the (gridded) verification set where \code{xhat1} and \code{xhat2} are the two forecast models being compared.  \code{plot.spatMLD}: \code{x} is a list returned by \code{spatMLD}.
}
  \item{object}{
 \code{fit.spatMLD} this is the output returned by \code{spatMLD}.  \code{summary.spatMLD}: list object returned by \code{spatMLD} or \code{fit.spatMLD}.
}
  \item{lossfun}{
character anming a loss function to use in finding the loss differential for the fields.  Default is to use correlation as the loss function.  Must have arguments \code{x} and \code{y}, and may have any additional arguments.
}
  \item{trend}{character saying \dQuote{ols} if it is desired to have the spatial trend of the loss differential estimated by ordinary least squares, or a matrix (of appropriate dimension) or single numeric giving the value of the spatial trend, or anything else (e.g., 0, NA or NULL) if trend is assumed to be constant.}
  \item{loc}{(optional) mn by 2 matrix giving location coordinates for each grid point.  Only used if trend is \dQuote{ols}.  If NULL, they are taken to be the grid expansion of the dimension of \code{x} (i.e., cbind(rep(1:dim(x)[1],dim(x)[2]), rep(1:dim(x)[2],each=dim(x)[1]))). }
  \item{maxrad}{
numeric giving the maximum radius for finding variogram differences per the \code{R} argument of \code{vgram.matrix}.
}
  \item{dx, dy}{
\code{dx} and \code{dy} of \code{vgram.matrix}.
}
  \item{zero.out}{logical, should the variogram be computed only over non-zero values of the process?  If TRUE, a modified version of \code{vgram.matrix} is used (\code{variogram.matrix}).}
  \item{vgmodel}{character string naming a variogram model function to use.  Default is the exponential variogram, \code{expvg}. }
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
    \item{model}{numeric indicating which forecast model to select for the analysis.}
    \item{icol}{(optional) color scheme.}
  \item{\dots}{
\code{spatMLD}: optional additional arguments to \code{lossfun}.  Not used by the summary or plot functions.
}
}
\details{
Hering and Genton (2011) introduce a test procedure for comparing spatial fields.  First, a loss function, g(x,y), is calculated, which can be any appropriate loss function.  This is calculated for each of two forecast fields.  The loss differential field is then given by:

		D(s) = g(x(s),y1(s)) - g(x(s),y2(s)), where s are the spatial locations, x is the verification field, and y1 and y2 are the two forecast fields.

It is assumed that D(s) = phi(s) + psi(s), where phi(s) is the mean trend and psi(s) is a mean zero stationary process with unknown covariance function C(h) = cov(psi(s),psi(s+h)).  In particular, the argument trend represents phi(s), and the default is that the mean is equal (and zero) over the entire domain.  If it is believed that this is not the case, then it should be removed before finding the covariance.  Currently, trend estimation is performed via lm, but it is also allowed to remove the trend using some other trend in numeric form of appropriate dimension (e.g., 1 or m by n, or something else that is allowed for M - N, where M is an m by n matrix).  To estimate the trend in another way, see e.g. Hering and Genton (2011) and references therein.

A test is constructed to test the null hypothesis of equal predictive ability on average.  That is,

		H_0: 1/|D| int_D E[D(s)]ds = 0, where |D| is the area of the domain, 

The test statistic is given by

		S_V = mean(D(s))/sqrt(mean(C(h))),

where C(h) = gamma(infinity|p) - gamma(h|p) is a fitted covariance function for the loss differential field.  The test statistic is assumed to be N(0,1) so that if the p-value is smaller than the desired level of significance, the null hypothesis is not accepted.

For 'fit.spatMLD', an exponential variogram is used. Specifically,

		gamma(h | theta=(s,r)) = s^2*(1 - exp(-h/r)),

where s is sqrt(sill) and r is the range (nugget effects are not accounted for here).  If \code{fit.spatMLD} should fail, and the empirical variogram appears to be reasonable (e.g., use the \code{plot} method function on \code{spatMLD} output to check that the empirical variogram is concave), then try giving alternative starting values for the \code{nls} function by using the \code{start.list} argument.  The default is to use the variogram value for the shortest separation distance as an initial estimate for s, and \code{maxrad} as the initial estimate for r.

Currently, it is not possible to fit other variogram models with this function.  Such flexibility may possibly be added in a future release.  In the meantime, use \code{fit.spatMLD} as a template to make your own similar function; just be sure to return an object of class \dQuote{nls}, and it should work seamlessly with the \code{plot} and \code{summary} method functions for a \dQuote{spatMLD} object.  For example, if it is desired to include the nugget or an extra factor (e.g., 3 as used in Hering and Genton, 2011), then a new similar function would need to be created.

Also, although the testing procedure can be applied to irregularly spaced locations (non-gridded), this function is set up only for gridded fields in order to take advantage of computational efficiencies (i.e., use of vgram.matrix), as these are the types of verification sets in mind for this package.  Eventually, a similar function will be added for non-gridded fields.

The above test assumes constant spatial trend.  It is possible to remove any spatial trend in D(s) before applying the test.

The actual test statistic is computed by the \code{summary} function.  Isotropy can be checked by the plot in the lower right panel of the result of the \code{plot} method function.  The function \code{nls} is used to fit the variogram model.

For application to precipitation fields, and introduction to the image warp (coming soon) and distance map loss functions, see Gilleland (2013).

}
\value{
  A list object is returned with possible components:
  \item{data.name}{character vector naming the fields under comparison}
  \item{lossfun,lossfun.args,vgram.args }{same as the arguments input to the spatMLD function.}
  \item{d}{m by n matrix giving the loss differential field, D(s).}
  \item{trend}{If trend is \dQuote{ols}, a list object of class \dQuote{lm}.  If trend is a numeric, it is the same as the value passed in.}
  \item{loc}{If trend is \dQuote{ols} or a numeric, then this is the same as the argument passed in, or if NULL, it is the expanded grid coordinates.}
  \item{lossdiff.vgram}{list object as returned by vgram.matrix}
  \item{vgmodel}{list object as returned by nls containing the fitted exponential variogram model where s is the estimate of sqrt(sill), and r of the range parameter (assuming 'fit.spatMLD' was used to fit the variogram model).}

summary.spatMLD invisibly returns the same list object as above with additional components:
  \item{Dbar}{the estimated mean loss differential (over the entire field).}
  \item{test.statistic}{the test statistic.}
  \item{p.value}{list object with components: two.sided--the two-sided alternative hypothesis--, less--the one-sided alternative hypothesis that the true value mu(D) < 0--and greater--the one-sided alternative hypothesis that mu(D) > 0--, p-values under the assumption of standard normality of the test statistic.}
}
\references{
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
grid<- list( x= seq( 0,5,,25), y= seq(0,5,,25))
obj<-Exp.image.cov( grid=grid, theta=.5, setup=TRUE)
look<- sim.rf( obj)
look[ look < 0] <- 0
look <- zapsmall( look)
     
look2 <- sim.rf( obj)*.25
look2[ look2 < 0] <- 0
look2 <- zapsmall( look2)

look3 <- sim.rf( obj)*2 + 5
look3[ look3 < 0] <- 0 
look3 <- zapsmall( look3)

res <- spatMLD(x=look, xhat1=look2, xhat2=look3, lossfun="abserrloss", maxrad=8)
res <- fit.spatMLD(res)
res <- summary(res)
plot(res)

\dontrun{
data(pert000)
data(pert004)
data(pert006)
data(ICPg240Locs)

hold <- make.SpatialVx(pert000, list(pert004, pert006), loc=ICPg240Locs,
    projection=TRUE, map=TRUE, loc.byrow = TRUE,
    field.type="Precipitation", units="mm/h",
    data.name=c("ICP Perturbed Cases", "pert000", "pert004", "pert006"))

look <- spatMLD(hold, lossfun="abserrloss", maxrad=8)
look <- fit.spatMLD(look)
plot(look)
summary(look)
}
}
\keyword{ htest }
