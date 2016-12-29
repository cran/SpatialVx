\name{MCdof}
\alias{MCdof}
\alias{sig.cor.t}
\alias{sig.cor.Z}
\alias{fisherz}
\title{
Monte Carlo Degrees of Freedom
}
\description{
Estimate the distribution of the proportion of spatial locations that contain significant correlations with randomly generated data along the lines of Livezey and Chen (1983).
}
\usage{
MCdof(x, ntrials = 5000, field.sig = 0.05, zfun = "rnorm", zfun.args = NULL,
      which.test = c("t", "Z", "cor.test"), verbose = FALSE, ...)

sig.cor.t(r, len = 40, ...)

sig.cor.Z(r, len = 40, H0 = 0)

fisherz(r)
}
\arguments{
  \item{x}{
n by m numeric matrix whose rows represent temporal points, and whose columns are spatial locations.
}
  \item{ntrials}{
numeric/integer giving the number of times to generate random samples of size n, and correlate them with the columns of \code{x}.
}
  \item{field.sig}{
numeric between 0 and 1 giving the desired fields significance level.
}
  \item{zfun}{
character naming a random number generator that takes \code{n} (the size of the sample to be drawn) as an argument, and any other arguments necessary.
}
  \item{zfun.args}{
list object giving the values for additional arguments to the function named by \code{zfun}.
}
  \item{which.test}{
character naming which type of test to do (default, \dQuote{t}, is a t-test, calls \code{sig.cor.t}). \dQuote{Z} does Fisher's Z transform (calls \code{sig.cor.Z}).   \dQuote{cor.test} calls \code{cor.test} giving more options, but is also considerably slower than \dQuote{t} or \dQuote{Z}.
}
  \item{r}{numeric giving the correlation value(s).}
  \item{len}{numeric giving the size of the data for the test.}
  \item{H0}{numeric giving the null hypothesis value (not used by \code{MCdof}).}
  \item{verbose}{
logical, should progress information (including total run time) be printed to the screen?
}
  \item{\dots}{
optional arguments to \code{sig.cor.t} (not used), \code{sig.cor.Z}, or \code{cor.test} depending on argument \code{which.test}.
}
}
\details{
This function does the Livezey and Chen (1983) Monte Carlo step 2 (a) from Elmore et al. (2006).  It generates a random sample of size n, and finds the p-values of a correlation test with this random sample and each column of \code{x}.  From this, it estimates the proportion of spatial locations that could contain significant bias purely by chance.
}
\value{
MCdof returns a list object with components:
  \item{MCprops }{numeric vector of length ntrials giving the proportion of locations with significant bias found by chance for each repition of the experiment.}
  \item{minsigcov }{ single numeric giving the 1 - field.sig quantile of the resulting proportions given by MCprops.}

sig.cor.t and sig.cor.Z return umeric vectors of p-values, and fisherz returns a numeric vector of test statistics.
}
\references{
Elmore, K. L., Baldwin, M. E. and Schultz, D. M. (2006) Field significance revisited: Spatial bias errors in forecasts as applied to the Eta model.  \emph{Mon. Wea. Rev.}, \bold{134}, 519--531.

Livezey, R. E. and Chen, W. Y. (1983) Statistical field significance and its determination by Monte Carlo techniques.  \emph{Mon. Wea. Rev.}, \bold{111}, 46--59.
}
\author{
Kimberly L. Elmore, Kim.Elmore \dQuote{at} noaa.gov, and Eric Gilleland
}

\seealso{
\code{\link{spatbiasFS}}, \code{\link{LocSig}}, \code{\link{cor.test}}, \code{\link{rnorm}}, \code{\link{runif}}, \code{\link{rexp}}, \code{\link{rgamma}}
}
\examples{
data( "GFSNAMfcstEx" )
data( "GFSNAMobsEx" )
data( "GFSNAMlocEx" )

id <- GFSNAMlocEx[,"Lon"] >=-90 & GFSNAMlocEx[,"Lon"] <= -75 & GFSNAMlocEx[,"Lat"] <= 40

look <- MCdof(GFSNAMfcstEx[,id] - GFSNAMobsEx[,id], ntrials=500)

stats(look$MCprops)
look$minsigcov

fisherz( abs(cor(rnorm(10),rexp(10), use="pairwise.complete.obs")))

}
\keyword{ distribution }
\keyword{ univar }
\keyword{ htest }
\keyword{ nonparametric }

