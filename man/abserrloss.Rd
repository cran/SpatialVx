\name{abserrloss}
\alias{abserrloss}
\alias{corrskill}
\alias{sqerrloss}
\alias{distmaploss}
\title{
Loss functions for the spatial prediction comparison test (SPCT)
}
\description{
Loss functions for applying the spatial prediction comparison test (SPCT) for competing forecasts.
}
\usage{
abserrloss(x, y, ...)
corrskill(x, y, ...)
sqerrloss(x, y, ...)
distmaploss(x, y, threshold = 0, const = Inf, ...)
}
\arguments{
  \item{x,y}{
	m by n numeric matrices against which to calculate the loss (or skill) functions.
}
  \item{threshold}{numeric giving the threshold over which (and including) binary fields are created from \code{x} and \code{y} in order to make a distance map.}
  \item{const}{numeric giving the constant beyond which the differences in distance maps between \code{x} and \code{y} are set to zero.  If \code{Inf} (default), then no cut-off is taken.  The SPCT is probably not powerful for large values of \code{const}.}
  \item{\dots}{
   Not used by \code{abserrloss} or \code{sqerrloss} (there for consistency only, and in order to work with \code{lossdiff}).  For \code{corrskill}, these are optional arguments to \code{sd}.  For \code{distmaploss}, these are optional arguments to the \code{distmap} function from pacakge \pkg{spatstat}.
}
}
\details{
These are simple loss functions that can be used in conjunction with \code{lossdiff} to carry out the spatial prediction comparison test (SPCT) as introduced in Hering and Genton (2011); see also Gilleland (2013) in particular for details about the distance map loss function.

The distance map loss function does not zero-out well as the other loss functions do.  Therefore, \code{zero.out} should be \code{FALSE} in the call to \code{lossdiff}.  Further, as pointed out in Gilleland (2013), the distance map loss function can easily be hedged by having a lot of correct negatives.  The image warp loss function is probably better for this purpose if, e.g., there are numerous zero-valued grid points in all fields.
}
\value{
numeric m by n matrices containing the value of the loss (or skill) function at each location i of the original set of locations (or grid of points).
}
\references{
Gilleland, E. (2013) Testing competing precipitation forecasts accurately and efficiently: The spatial prediction comparison test.  \emph{Mon. Wea. Rev.}, \bold{141}, (1), 340--355.

Hering, A. S. and Genton, M. G. (2011) Comparing spatial predictions.  \emph{Technometrics} \bold{53}, (4), 414--425.

}
\author{
Eric Gilleland
}
\seealso{
\code{\link{lossdiff}}, \code{\link{vgram.matrix}}, \code{\link{vgram}}
}
\examples{
# See help file for lossdiff for examples.
}
\keyword{ math }
