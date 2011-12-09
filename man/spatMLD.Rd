\name{spatMLD}
\alias{spatMLD}
\alias{fit.spatMLD}
\alias{summary.spatMLD}
\alias{plot.spatMLD}
\title{
Test for equal predictive ability on average over a regularly gridded space
}
\description{
Test for equal predictive ability (for two forecast models) on average over a regularly gridded space using the method of Hering and Genton (2011).
}
\usage{
spatMLD(x, y1, y2, lossfun = "corrskill", trend = "ols", loc = NULL, maxrad = 20, dx = 1, dy = 1, zero.out = FALSE, ...)
fit.spatMLD(object)
\method{summary}{spatMLD}(object, ...)
\method{plot}{spatMLD}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y1,y2}{
\code{spatMLD}: m X n matrices defining the (gridded) verification set where \code{y1} and \code{y2} are the two forecast models being compared.  \code{plot.spatMLD}: \code{x} is a list returned by \code{spatMLD}.
}
  \item{object}{
 \code{fit.spatMLD} this is the output returned by \code{spatMLD}.  \code{summary.spatMLD}: list object returned by \code{spatMLD} or \code{fit.spatMLD}.
}
  \item{lossfun}{
character anming a loss function to use in finding the loss differential for the fields.  Default is to use correlation as the loss function.  Must have arguments \code{x} and \code{y}, and may have any additional arguments.
}
  \item{trend}{character saying "ols" if it is desired to have the spatial trend of the loss differential estimated by ordinary least squares, or a matrix (of appropriate dimension) or single numeric giving the value of the spatial trend, or anything else (e.g., 0, NA or NULL) if trend is assumed to be constant.}
  \item{loc}{(optional) mn X 2 matrix giving location coordinates for each grid point.  Only used if trend is "ols".  If NULL, they are taken to be the grid expansion of the dimension of \code{x} (i.e., cbind(rep(1:dim(x)[1],dim(x)[2]), rep(1:dim(x)[2],each=dim(x)[1]))). }
  \item{maxrad}{
numeric giving the maximum radius for finding variogram differences per the \code{R} argument of \code{vgram.matrix}.
}
  \item{dx, dy}{
\code{dx} and \code{dy} of \code{vgram.matrix}.
}
  \item{zero.out}{logical, should the variogram be computed only over non-zero values of the process?  If TRUE, a modified version of \code{vgram.matrix} is used (\code{variogram.matrix}).}
  \item{\dots}{
\code{spatMLD}: optional additional arguments to \code{lossfun}.  Not used by the summary or plot functions.
}
}
\details{
Hering and Genton (2011) introduce a test procedure for comparing spatial fields.  First, a loss function, g(x,y), is calculated, which can be any appropriate loss function.  This is calculated for each of two forecast fields.  The loss differential field is then given by:

		D(s) = g(x(s),y1(s)) - g(x(s),y2(s)), where s are the spatial locations, x is the verification field, and y1 and y2 are the two forecast fields.

It is assumed that D(s) = phi(s) + psi(s), where phi(s) is the mean trend and psi(s) is a mean zero stationary process with unknown covariance function C(h) = cov(psi(s),psi(s+h)).  In particular, the argument trend represents phi(s), and the default is that the mean is equal (and zero) over the entire domain.  If it is believed that this is not the case, then it should be removed before finding the covariance.  Currently, trend estimation is performed via lm, but it is also allowed to remove the trend using some other trend in numeric form of appropriate dimension (e.g., 1 or m X n, or something else that is allowed for M - N, where M is an m X n matrix).  To estimate the trend in another way, see e.g. Hering and Genton (2011) and references therein.

A test is constructed to test the null hypothesis of equal predictive ability on average.  That is,

		H_0: 1/|D| int_D E[D(s)]ds = 0, where |D| is the area of the domain, 

The test statistic is given by

		S_V = mean(D(s))/sqrt(mean(C(h))),

where C(h) = gamma(infinity|p) - gamma(h|p) is a fitted covariance function for the loss differential field.  The test statistic is assumed to be N(0,1) so that if the p-value is smaller than the desired level of significance, the null hypothesis is not accepted.

For this function, an exponential variogram is used (at some point, this function will be amended to allow general variogram models).  Specifically,

		gamma(h|sigma,theta) = sigma^2*(1 - exp(-h/theta))

Also, although the testing procedure can be applied to irregularly spaced locations (non-gridded), this function is set up only for gridded fields in order to take advantage of computational efficiencies (i.e., use of vgram.matrix), as these are the types of verification sets in mind for this package.  Eventually, a similar function will be added for non-gridded fields.

The above test assumes constant spatial trend.  It is possible to remove any spatial trend in D(s) before applying the test.

The actual test statistic is computed by the summary function.  Isotropy can be checked by the plot in the lower right panel of the result of the plot method function.  The function nls is used to fit the variogram model.

For application to precipitation fields, and introduction to the image warp (coming soon) and distance map loss functions, see Gilleland (2012).
}
\value{
  A list object is returned with possible components:
  \item{Vx.name,Mod1.name,Mod2.name }{characters naming the fields under comparison}
  \item{lossfun,lossfun.args,vgram.args }{same as the arguments input to the spatMLD function.}
  \item{d}{m X n matrix giving the loss differential field, D(s).}
  \item{trend}{If trend is "ols", a list object of class "lm".  If trend is a numeric, it is the same as the value passed in.}
  \item{loc}{If trend is "ols" or a numeric, then this is the same as the argument passed in, or if NULL, it is the expanded grid coordinates.}
  \item{lossdiff.vgram}{list object as returned by vgram.matrix}
  \item{vgmodel}{list object as returned by nls containing the fitted exponential variogram model where s is the estimate of sigma, and r of theta (i.e., the range).}

summary.spatMLD invisibly returns the same list object as above with additional components:
  \item{Dbar}{the estimated mean loss differential (over the entire field).}
  \item{test.statistic}{the test statistic.}
  \item{p.value}{the (two-sided) p-value under the assumption of standard normality of the test statistic.}
}
\references{
Gilleland, E., 2012: Testing competing precipitation forecasts accurately and efficiently: The spatial prediction comparison test.  Submitted to \emph{Mon. Wea. Rev.}

Hering, A. S. and M. G. Genton, 2011: Comparing spatial predictions.  \emph{Technometrics} \bold{53}, (4), 414 - 425.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{vgram.matrix}}, \code{\link{nls}}, \code{\link{corrskill}}, \code{\link{abserrloss}}, \code{\link{sqerrloss}}, \code{\link{distmaploss}}
}
\examples{
data(pert000)
data(pert004)
data(pert006)
look <- spatMLD(x=pert000, y1=pert004, y2=pert006, lossfun="abserrloss", maxrad=8)
look <- fit.spatMLD(look)
plot(look)
summary(look)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
