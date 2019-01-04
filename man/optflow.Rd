\name{optflow}
\alias{optflow}
\title{
Optical Flow
}
\description{
Estimate the optical flow from one gridded field (image) to another.
}
\usage{
optflow(initial, final, grads.diff = 1, mean.field = NULL, ...)
}
\arguments{
  \item{initial,final}{m by n matrices where the optical flow is determined from initial (forecast) to final (observation).}
  \item{grads.diff}{either 1 or 2, where 1 calculates first derivatives with first differences and 2 first derivatives with second differences.}
  \item{mean.field}{Should they first be centered?  If so, give the value for the centering here (usually the mean of initial).}
  \item{\dots}{optional arguments to the \code{optim} function (cannot be \code{par}, \code{fn}, \code{gr} or \code{method}).}
}
\details{
This function estimates the optical flow from the initial field (image) to the final one as described in Marzban and Sandgathe (2010).  Letting I_o(x,y) and I_f(x,y) represent the intensities of each field at coordinate (x,y), the collection of pairs (dx, dy) is the optical flow field, where:

I_o(x,y) ~ I_f(x,y) + [partial(I_f) wrt x]*dx + [partial(I_f) wrt y]*dy.

The procedure follows that proposed by Lucas and Kanade (1981) whereby for some window, W, it is assumed that all dx (dy) are assumed constant, and least squares estimation is used to estimate dx and dy (see Marzban and Sandgathe, 2010 for more on this implementation).  It is assumed that the fields (initial and final) include only the window around the point of interest (i.e., this function finds the optical flow estimate for a single window).  See the function \code{OF}, which iteratively calls this function, for performing optical flow over the entire field.

The above formulation is linear in the parameters.  Marzban and Sandgathe (2010) also introduce an additive error component, which leads to a nonlinear version of the above.  Namely,

I_o(x,y) ~ I_f(x,y) + [partial(I_f) wrt x]*dx + [partial(I_f) wrt y]*dy + A(x,y).

See Marzban and Sandgathe for more details.
}
\value{
numeric vector whose first three components are the optimized estimates (returned by  the par component of optim) for the regression I_o(x,y) - I_f(x,y) = a0 + a1*[partial(I_f) wrt x] + a2*[partial(I_f) wrt y] (i.e., a1 and a2 are the estimates for dx and dy, resp.) and the latter three values are the initial estimates to optim as determined by linear regression (i.e., returned from the lm function).
}
\references{
Lucas, B D. and Kanade, T. (1981)  An iterative image registration technique with an application to stereo vision.  \emph{Proc. Imaging Understanding Workshop}, DARPA, 121--130.

Marzban, C. and Sandgathe, S. (2010)  Optical flow for verification.  \emph{Wea. Forecasting}, \bold{25}, 1479--1494, doi:10.1175/2010WAF2222351.1.
}
\author{
Caren Marzban, marzban \dQuote{at} u.washington.edu, and modified by Eric Gilleland
}

\seealso{
\code{\link{OF}}, \code{\link{optim}}, \code{\link{lm}}
}
\examples{
x <- y <- matrix(0, 10, 10)
x[1:2,3:4] <- 1
y[3:4,5:6] <- 2

optflow(x,y)

\dontrun{
initial <- hump$initial
final <- hump$final
look <- OF(final, initial, W=9, verbose=TRUE)
plot(look) # Compare with Fig. 1 in Marzban and Sandgathe (2010).
hist(look) # 2-d histogram.
plot(look, full=TRUE) # More plots.
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ regression }% __ONLY ONE__ keyword per line
