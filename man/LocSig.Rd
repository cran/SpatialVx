\name{LocSig}
\alias{LocSig}
\alias{plot.LocSig}
\title{
Temporal Block Bootstrap Keeping Locations in Space Constant
}
\description{
Temporal block bootstrap for data at spatial locations (holding locations constant at each iteration).  This is a wrapper function to the tsboot or boot functions for use with the field significance approach of Elmore et al. (2006).
}
\usage{
LocSig(Z, numrep = 1000, block.length = NULL, bootfun = "mean",
    alpha = 0.05, bca = FALSE, ...)

\method{plot}{LocSig}(x, loc = NULL, nx = NULL, ny = NULL, ...)
}
\arguments{
  \item{Z}{
n by m numeric matrix whose rows represent contiguous time points, and whose columns represent spatial locations.
}
  \item{numrep}{
numeric/integer giving the number of bootstrap replications to use.
}
  \item{block.length}{
positive numeric/integer giving the desired block lengths.  If NULL, \code{floor(sqrt(n))} is used.  If 1, then the IID bootstrap is performed, and the BCa method may be used to find CI's, if \code{bca} is TRUE.
}
  \item{bootfun}{
character naming an R function to be applied to each replicate sample.  Must return a single number, but is otherwise the \code{statistic} argument for function \code{tsboot} (or \code{boot} if \code{block.length} = 1).
}
  \item{alpha}{
numeric giving the value of \code{alpha} to obtain (1-\code{alpha})*100 percent CI's for \code{bootfun}.
}
  \item{bca}{
logical, should bias-corrected and adjusted (BCa) CI's be calculated?  Only used if \code{block.length} = 1.  Will give a warning if this argument is TRUE, and \code{block.length} > 1, and will use the percentile method.
}
  \item{x}{ data frame of class \dQuote{LocSig} as returned by \code{LocSig}.}
  \item{loc}{ m by 2 matrix of location coordinates.}
  \item{nx,ny}{ If \code{loc} is NULL, then \code{nx} and \code{ny} must be supplied.  These give the number of rows and columns of a grid to make an image (using \code{as.image}) for plotting.  If these are used, the data \code{Z} must be from a regular grid of points.}
  \item{\dots}{
\code{LocSig}: optional additional arguments to the \code{tsboot} (or \code{boot} if \code{block.length}=1) function.
\code{plot.LocSig}: optional additional arguments to \code{image.plot}.
}
}
\details{
This function performs the circular block bootstrap algorithm over time at each of m locations (columns of \code{x}).  So, at each bootstrap iteration, entire blocks of rows of x are resampled with replacement.  If \code{Z} represents forecast errors at grid points, and \code{bootfun}=\dQuote{mean}, then this finds the grid-point CI's in steps 1 (a) to 1 (c) of Elmore et al. (2006).
}
\value{
LocSig: A data frame with class attribute \dQuote{LocSig} with components:
  \item{Estimate }{numeric giving the estimated values of bootfun (the statistic for which CI's are computed).}
  \item{Lower, Upper }{numeric giving the estimated lower (upper) (1-alpha)*100 percent CI's.}

plot.LocSig: invisibly returns a list containing the estimate as returned by LocSig, and the confidence range.
}
\references{
Elmore, K. L., Baldwin, M. E. and Schultz, D. M. (2006) Field significance revisited: Spatial bias errors in forecasts as applied to the Eta model.  \emph{Mon. Wea. Rev.}, \bold{134}, 519--531.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{spatbiasFS}}, \code{\link{tsboot}}, \code{\link{boot}}, \code{\link{boot.ci}}, \code{\link{MCdof}}, \code{\link{sig.cor.t}}, \code{\link{sig.cor.Z}}, \code{\link{cor.test}}, \code{\link{image.plot}}, \code{\link{as.image}}
}
\examples{
\dontrun{
data( "GFSNAMfcstEx" )
data( "GFSNAMobsEx" )
data( "GFSNAMlocEx" )

id <- GFSNAMlocEx[,"Lon"] >=-90 & GFSNAMlocEx[,"Lon"] <= -75 & GFSNAMlocEx[,"Lat"] <= 40

look <- LocSig(GFSNAMfcstEx[,id] - GFSNAMobsEx[,id], numrep=500)

stats(look)

plot(look, loc = GFSNAMlocEx[ id, ] )
}
}
\keyword{ distribution }
\keyword{ univar }
\keyword{ htest }
\keyword{ nonparametric }
