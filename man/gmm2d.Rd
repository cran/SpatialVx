\name{gmm2d}
\alias{gmm2d}
\alias{gmm2d.default}
\alias{gmm2d.SpatialVx}
\alias{plot.gmm2d}
\alias{predict.gmm2d}
\alias{print.gmm2d}
\alias{summary.gmm2d}
\title{
2-d Gaussian Mixture Models Verification
}
\description{
Use 2-d Gaussian Mixture Models (GMM) to assess forecast performance.
}
\usage{
gmm2d(x, ...)

\method{gmm2d}{default}(x, ..., xhat, K = 3, gamma = 1, threshold = NULL,
    initFUN = "initGMM", verbose = FALSE)

\method{gmm2d}{SpatialVx}(x, ..., time.point = 1, obs = 1, model = 1, K = 3, gamma = 1,
    threshold = NULL, initFUN = "initGMM", verbose = FALSE)

\method{plot}{gmm2d}(x, ..., col = c("gray", tim.colors(64)),
    zlim = c(0, 1), horizontal = TRUE)

\method{predict}{gmm2d}(object, ..., x)

\method{print}{gmm2d}(x, ...)

\method{summary}{gmm2d}(object, ...)
}
\arguments{
  \item{x,xhat}{
    Default: m by n numeric matrices giving the verification and forecast fields, resp.

    \code{gmm2d.SpatialVx}: object of class \dQuote{SpatialVx}.

    \code{plot} and \code{print}: object returned by \code{gmm2d}.

    \code{predict}: k by 2 matrix of lon/lat coordinates on which to predict the model.
}
  \item{object}{output from \code{gmm2d}.}
  \item{K}{
	single numeric giving the number of mixture components to use.
}
  \item{gamma}{
Value of the gamma parameter from Eq (11) of Lakshmanan and Kain (2010).  This affects the number of times a location is repeated.
}
  \item{threshold}{
numeric giving a threshold over which (and including) the GMM is to be fit (zero-valued grid points are not included in the estimation here for speed).  If NULL, no thresholding is applied.
}
  \item{initFUN}{
character naming a function to provide initial estimates for the GMM.  Must take an m by n matrix as input, and return a dataframe a component called \code{ind} that is a vector indicating the order of the rows for which the first \code{K} will be used, a third column giving the x-coordinates of the initial estimate of the mean for the x direction, fourth column giving the initial estimate for the mean of the  y-direction, and fifth and sixth columns giving initial estimates for the standard deviations of the x- and y-directions.  The default identifies all connected components using the \code{disjointer} function, then uses their centroids as the initial estimates of the means, and their axes as initial estimates for the standard deviations.  The \code{ind} component gives the order of the object areas from largest to smallest so that the \code{K} largest objects are used to provide initial estimates.  Note that this differs from the initial estimates in Lakshmanan and Kain (2010) where they break the field into different areas first.
}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
    \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
	\item{col, zlim, horizontal}{optional arguments to \pkg{fields} function(s) \code{poly.image}, \code{image.plot}.}
  \item{\dots}{
In the case of \code{gmm2d}: optional arguments to \code{initFUN}.  In the case of \code{plot}: not used.  In the case of \code{predict}: N by 2 matrix of grid point locations on which to predict the probability from the 2-d GMM model.  In the case of \code{summary}: this can include the arguments: \code{silent}, logical stating whether to print summaries to the screen (FALSE) or not (TURE), \code{e1}, \code{e2}, ..., \code{e5}, giving alternative weights in calculating the overall error (Eq 15 in Lakshmanan and Kain, 2010, but see details section below).
}
}
\details{
   These functions carry out the spatial verification approach described in Lakshmanan and Kain (2010), which fits a 2-d Gaussian Mixture Model (GMM) to the locations for each field in the verification set, and makes comparisons using the estimated parameters.  In fitting the GMMs, first an initial estimate is provided by using the initFUN argument, which is a function.  The default function is relatively fast (it might seem slow, but for what it does, it is very fast!), but is typically the slowest part of the process.  Although the EM algorithm is a fairly computationally intensive procedure, acceleration algorithms are employed (via the turboem function of the turboEM package) so that once initial estimates are found, the procedure is very fast.

Because the fit is to the locations only, Lakshmanan and Kain (2010) suggest two ways to incorporate intensity information.  The first is to repeat points with higher intensities, and the second is to multiply the results by the total intensities over the fields.  The points are repeated M times according to the formula (Eq 11 in Lakshmanan and Kain, 2010):

	M = 1 + gamma * round( CFD(I_{xy})/frequency(I_MODE)),

where CFD is the cumulative *frequency* distribution (here estimated from the histogram using the \sQuote{hist} function), I_{xy} is intensity at grid point (x,y), I_MODE is the mode of intensity values, and gamma is a user-supplied parameter controlling how much to repeat points where higher numbers will result in larger repetitions of high intensity values.

The function \code{gmm2d} fits the 2-d GMM to both fields, \code{plot.gmm2d} first uses \code{predict.gmm2d} to obtain probabilities for each grid point, and then makes a plot similar to those in Lakshmanan and Kain (2010) Figs. 3, 4 and 5, but giving the probabilities instead of the probabilities times A.  Note that \code{predict.gmm2d} can be very slow to compute so that \code{plot.gmm2d} can also be very slow.  Less effort was put into speeding these functions up because they are not necessary for obtaining results via the parameters.  However, they can give the user an idea of how good the fit is.

The 2-d GMM is given by

G(x,y) = A*sum(lambda*f(x,y))

where lambda and f(x,y) are numeric vectors of length K, lambda components describe the mixing, and f(x,y) is the bivariate normal distribution with mean (mu.x, mu.y) and covariance function.  \sQuote{A} is the total sum of intensities over the field.

Comparisons between forecast and observed fields are carried out finally by the summary method function.  In particular, the translation error

	e.tr = sqrt((mu.xf - mu.xo)^2 + (mu.yf - mu.yo)^2),

where f means forecast and o verification fields, resp., and mu .x is the mean in the x- direction, and mu.y in the y- direction.  The rotation error is given by

	e.rot = (180/pi)*acos(theta),

where theta is the dot product between the first eigenvectors of the covariance matrices for the verification and forecast fields.  The scaling error is given by

	e.sc = Af*lambda.f/Ao*lambda.o,

where lambda is the mixture component and Af/Ao is the forecast/observed total intensity.

The overall error (Eq 15 of Lakshmanana and Kain, 2010) is given by

	e.overall = e1 * min(e.tr/e2, 1) + e3*min(e.rot,180 - e.rot)/e4 + e5*(max(e.sc,1/e.sc)-1),

where e1 to e5 can be supplied by the user, but the defaults are those given by Lakshmanan and Kain (2010).  Namely, e1 = 0.3, e2 = 100, e3=0.2, e4 = 90, and e5=0.5.
}
\value{
For gmm2d, a list object of class \dQuote{gmm2d} is returned with components:
\item{fitX,fitY}{list objects returned by the \code{turboem} function from the \pkg{turboEM} package that describe the EM estimates of the 2-d GMM parameters for the verification and forecast fields, resp.}
\item{initX,initY }{numeric vectors giving the initial estimates used in the EM algorithm for the verification and forecast fields, resp.  The first 2*K values are the initial mean estiamtes for the x- and y- directions, resp.  The next 4*K values are the initial estiamtes of the covariances (note that the cross-covariance terms are zero regardless of initialization function employed (maybe this will be improved in the future).  The final K values are the initial estimates for lambda.}
\item{sX,sY}{ N by 2 matrix giving the repeated coordinates calculated per M as described in the details section for the verification and forecast fields, resp.}
\item{k}{single numeric giving the value of K}
\item{Ax,Ay}{single numerics giving the value of A (the total sum of intensities over the field) for the verifiaction and forecast fields, resp.}

For plot.gmm2d no value is returned.  A plot is created.

For predict.gmm2d, a list is returned with components:
\item{predX,predY}{numeric vectors giving the GMM predicted values for the verification and forecast fields, resp.}

For summary.gmm2d, a list is returned invisibly (if silent is FALSE, information is printed to the screen) with components:
\item{meanX,meanY}{Estimated mean vectors for each GMM component for the verification and forecast fields, resp.}
\item{covX,covY}{Estimated covariances for each GMM component for the verification and forecast fields, resp.}
\item{lambdasX,lambdasY}{Estimated mixture components for each GMM component for the verification and forecast fields, resp.}
\item{e.tr,e.rot,e.sc,e.overall}{K by K matrices giving the errors between each GMM component in the verification field (rows) to each GMM component in the forecast field (columns).  The errors are: translation (e.tr), rotation (e.rot), scaling (e.sc), and overall (e.overall).}
}
\references{
Lakshmanan, V. and Kain, J. S. (2010) A Gaussian Mixture Model Approach to Forecast Verification. \emph{Wea. Forecasting}, \bold{25} (3), 908--920.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{turboem}}, \code{\link{disjointer}}, \code{\link{connected}}
}
\examples{
\dontrun{
data( "ExampleSpatialVxSet" )

x <- ExampleSpatialVxSet$vx
xhat <- ExampleSpatialVxSet$fcst

u <- min(quantile(c(x[x > 0]), probs = 0.75),
    quantile(c(xhat[xhat > 0]), probs = 0.75))

look <- gmm2d(x, xhat=xhat, threshold=u, verbose=TRUE)
summary(look)
plot(look)
}
\dontrun{
# Alternative method to skin the cat.
hold <- make.SpatialVx( x, xhat, field.type = "MV Gaussian w/ Exp. Cov.",
    units = "units", data.name = "Example", obs.name = "x",
    model.name = "xhat" )

look2 <- gmm2d( hold, threshold = u, verbose = TRUE)
summary(look2)
plot(look2)

}
}
\keyword{math }
\keyword{ models }
\keyword{ univar}
