\name{expvg}
\alias{expvg}
\alias{predict.flossdiff.expvg}
\alias{print.flossdiff.expvg}
\title{
Exponential Variogram
}
\description{
Compute the exponential variogram.
}
\usage{
expvg(p, vg, ...)

\method{predict}{flossdiff.expvg}(object, newdata, ...)

\method{print}{flossdiff.expvg}(x, ...)
}
\arguments{
  \item{p}{
numeric vector of length two.  Each component should be positively valued.  The first component is the nugget and the second is the range parameter.
}
  \item{vg}{
A list object with component \code{d} giving a numeric vector of distances over which the variogram is to be calculated.
}
  \item{object, x}{A list object returned by \code{flossdiff} using \code{expvg} as the variogram model.}
  \item{newdata}{Numeric giving the distances over which to use the fitted exponential variogram model to make predictions.  The default is to go from zero to the maximum lag distance for a given data set, which is not the usual convention for the generic \code{predict}, which usually defaults to operate on the lags used in performing the fit.}
  \item{\dots}{
Not used.
}
}
\details{
A very simple function used mainly internally by \code{flossdiff} when fitting the exponential variogram to the empirical one, and by the \code{predict}, \code{print} and \code{summary} method functions for \code{lossdiff} objects.  For those wishing to use a different variogram model than the exponential, use this function and its method functions as a template.  Be sure to create \code{predict} and \code{print} method functions to operate on objects of class \dQuote{flossdiff.XXX} where \dQuote{XXX} is the name of the variogram function you write (so, \dQuote{expvg} in the current example).
}
\value{
Numeric vector of length equal to that of the \code{d} component of \code{vg} giving the corresponding exponential variogram values with nugget and range defined by \code{p}.
}
\references{
Cressie, N. A. (2015) Statistics for Spatial Data.  Wiley-Interscience; Revised Edition edition (July 27, 2015), ISBN-10: 1119114616, ISBN-13: 978-1119114611, 928 pp.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{lossdiff}}, \code{\link{flossdiff}}
}
\examples{
##
## For examples, see lossdiff and flossdiff
##
}
\keyword{ manip }
\keyword{ spatial }
