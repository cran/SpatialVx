\name{expvgram}
\alias{expvgram}
\title{
Exponential Variogram
}
\description{
Calculates the empirical variogram for use with function spct.
}
\usage{
expvgram(p, h, ...)
}
\arguments{
  \item{p}{ numeric vector of length two giving the nugget and range parameter values, resp.
}
  \item{h}{
numeric vector of separation distances.
}
  \item{\dots}{
Not used.
}
}
\details{
Simple function to work with \code{spct} to calculate the exponential variogram for given parameters and separation distances.  The exponential variogram employed here is parameterized by

gamma(h) = sigma * ( 1 - exp( - h * theta ) )

where \code{p} is \code{c( sigma, theta )}.
}
\value{
A numeric vector of variogram values for each separation distance in \code{h}.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{spct}}
}
\examples{
# See help file for spct for examples.
}
\keyword{ array }
