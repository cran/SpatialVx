\name{abserrloss}
\alias{abserrloss}
\alias{corrskill}
\alias{sqerrloss}
\alias{distmaploss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Loss functions for applying the spatial prediction comparison test (SPCT)
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y}{
	m X n numeric matrices against which to calculate the loss (or skill) functions.
}
  \item{threshold}{numeric giving the threshold over which (and including) binary fields are created from \code{x} and \code{y} in order to make a distance map.}
  \item{const}{numeric giving the constant beyond which the differences in distance maps between \code{x} and \code{y} are set to zero.  If 'Inf' (default), then no cut-off is taken.  The SPCT is probably not powerful for large values of \code{const}.}
  \item{\dots}{
   Not used by \code{abserrloss} or \code{sqerrloss} (there for consistency only, and in order to work with \code{spatMLD}).  For \code{corrskill}, these are optional arguments to \code{sd}.  For \code{distmaploss}, these are optional arguments to the \code{distmap} function from pacakge \pkg{spatstat}.
}
}
\details{
These are simple loss functions that can be used in conjunction with spatMLD to carry out the spatial prediction comparison test (SPCT) as introduced in Hering and Genton (2011); see also Gilleland (2012) in particular for details about the distance map loss function.
}
\value{
numeric m X n matrices containing the value of the loss (or skill) function at each location i of the original set of locations (or grid of points).
}
\references{
Gilleland, E., 2012: Testing competing precipitation forecasts accurately and efficiently: The spatial prediction comparison test.  Submitted to \emph{Mon. Wea. Rev.}

Hering, A. S. and M. G. Genton, 2011: Comparing spatial predictions.  \emph{Technometrics} \bold{53}, (4), 414 - 425.

}
\author{
Eric Gilleland
}
\seealso{
\code{\link{spatMLD}}, \code{\link{vgram.matrix}}, \code{\link{vgram}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, ...) 
return(abs(x - y))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
