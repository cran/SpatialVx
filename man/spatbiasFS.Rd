\name{spatbiasFS}
\alias{spatbiasFS}
\alias{plot.spatbiasFS}
\alias{summary.spatbiasFS}
\title{
Field Significance Method of Elmore et al. (2006)
}
\description{
Apply field significance method of Elmore et al. (2006).
}
\usage{
spatbiasFS(X, Y, loc = NULL, block.length = NULL, alpha.boot = 0.05, field.sig = 0.05,
    bootR = 1000, ntrials = 1000, verbose = FALSE)

\method{summary}{spatbiasFS}(object, ...)

\method{plot}{spatbiasFS}(x, ...)
}
\arguments{
  \item{X,Y}{
m by n matrices giving the verification and forecast fields, resp., for each of m time points (rows) and n locations (columns).
}
  \item{x,object}{list object as returned by \code{spatbiasFS}.}
  \item{loc}{
optional (for subsequent plotting) n by 2 matrix giving the lon/lat coordinates for the locations.
}
  \item{block.length}{
numeric giving the block length to be used n the block bootstrap algorithm.  If NULL, floor(sqrt(n)) is used.
}
  \item{alpha.boot}{
numeric between 0 and 1 giving the confidence level desired for the bootstrap algorithm.
}
  \item{field.sig}{
numeric between 0 and 1 giving the desired field significance level.
}
  \item{bootR}{
numeric integer giving the number of bootstrap replications to use.
}
  \item{ntrials}{
numeric integer giving the number of Monte Carol iterations to use.
}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
\item{\dots}{not used.}
}
\details{
See Elmore et al. (2006) for details.
}
\value{
A list object with components:
\item{data.name}{character vector giving the name of the verification and forecast spatio-temporal fields used, and the associated location object (if not NULL).}
\item{block.boot.results }{object of class LocSig}
\item{sig.results}{list object containing information about the significance of the results.}
\item{field.significance,alpha.boot}{field significance level and bootstrap CI level as input by field.sig alpha.boot arguments.}
\item{bootR,ntrials}{same as arguments above.}
}
\references{
Elmore, K. L., Baldwin, M. E. and Schultz, D. M. (2006) Field significance revisited: Spatial bias errors in forecasts as applied to the Eta model.  \emph{Mon. Wea. Rev.}, \bold{134}, 519--531.
}
\author{
Eric Gilleland and Kimberly L. Elmore
}

\seealso{
\code{\link{MCdof}}, \code{\link{LocSig}}, \code{\link{tsboot}}
}
\examples{
data(GFSNAMfcstEx)
data(GFSNAMobsEx)
data(GFSNAMlocEx)
id <- GFSNAMlocEx[,"Lon"] >=-95 & GFSNAMlocEx[,"Lon"] <= -75 & GFSNAMlocEx[,"Lat"] <= 32
loc <- GFSNAMlocEx[id,]
GFSobsSub <- GFSNAMobsEx[,id]
GFSfcstSub <- GFSNAMfcstEx[,id]
look <- spatbiasFS(GFSobsSub, GFSfcstSub, loc=loc, bootR=500, ntrials=500)
plot(look)
summary(look)
}
\keyword{ math }
