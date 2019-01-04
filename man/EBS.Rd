\name{EBS}
\alias{EBS}
\alias{plot.EBS}
\title{
Elmore, Baldwin and Schultz Method for Field Significance for Spatial Bias Errors
}
\description{
Apply the method of Elmore, Baldwin and Schultz (2006) for calculating field significance of spatial bias errors.
}
\usage{
EBS(object, model = 1, block.length = NULL, alpha.boot = 0.05,
    field.sig = 0.05, bootR = 1000, ntrials = 1000,
    verbose = FALSE)

\method{plot}{EBS}(x, ..., mfrow = c(1, 2), col, horizontal)
}
\arguments{
  \item{object}{ list object of class \dQuote{SpatialVx}.
}
    \item{x}{object of class \dQuote{EBS} as returned by \code{EBS}.}
  \item{model}{
number or character describing which model (if more than one in the \dQuote{SpatialVx} object) to compare.
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
    \item{mfrow}{mfrow parameter (see help file for \code{par}).  If NULL, then the parameter is not re-set.  }
    \item{col, horizontal}{optional arguments to \code{image.plot} from \pkg{fields}.}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
    \item{\dots}{optional arguments to \code{image.plot} from \pkg{fields}.}
}
\details{
this is a wrapper function for the \code{spatbiasFS} function utilizing the \dQuote{SpatialVx} object class to simplify the arguments.
}
\value{
A list object of class \dQuote{EBS} with the same attributes as the input object and additional attribute (called \dQuote{arguments})that is a named vector giving information provided by the user.  Components of the list include:
  \item{block.boot.results }{object of class \dQuote{LocSig}.}
  \item{sig.results }{list object containing information about the significance of the results.}
}
\references{
Elmore, K. L., Baldwin, M. E. and Schultz, D. M. (2006) Field significance revisited: Spatial bias errors in forecasts as applied to the Eta model.  \emph{Mon. Wea. Rev.}, \bold{134}, 519--531.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{boot}}, \code{\link{tsboot}}, \code{\link{spatbiasFS}}, \code{\link{LocSig}}, \code{\link{poly.image}}, \code{\link{image.plot}}, \code{\link{make.SpatialVx}}
}
\examples{
data( "GFSNAMfcstEx" )
data( "GFSNAMobsEx" )
data( "GFSNAMlocEx" )

id <- GFSNAMlocEx[,"Lon"] >=-95
id <- id & GFSNAMlocEx[,"Lon"] <= -75
id <- id & GFSNAMlocEx[,"Lat"] <= 32

##
## This next step is a bit awkward, but these data
## are not in the format of the SpatialVx class.
## These are being set up with arbitrarily chosen
## dimensions (49 X 48) for the spatial part.  It
## won't matter to the analyses or plots.
##
Vx <- GFSNAMobsEx
Fcst <- GFSNAMfcstEx
Ref <- array(t(Vx), dim=c(49, 48, 361))
Mod <- array(t(Fcst), dim=c(49, 48, 361)) 

hold <- make.SpatialVx(Ref, Mod, loc=GFSNAMlocEx,
    projection=TRUE, map=TRUE, loc.byrow = TRUE, subset=id,
    field.type="Precipitation", units="mm",
    data.name = "GFS/NAM", obs.name = "Reference", model.name = "Model" )

look <- EBS(hold, bootR=500, ntrials=500, verbose=TRUE)
plot( look )

\dontrun{
# Same as above, but now we'll do it for all points.
# A little slower, but not terribly bad.

hold <- make.SpatialVx(Ref, Mod, loc = GFSNAMlocEx,
    projection = TRUE, map = TRUE, loc.byrow = TRUE,
    field.type = "Precipitation", reg.grid = FALSE, units = "mm",
    data.name = "GFS/NAM", obs.name = "Reference", model.name = "Model" )

look <- EBS(hold, bootR=500, ntrials=500, verbose=TRUE)
plot( look )
}
}
\keyword{ htest }
\keyword{ nonparametric }
\keyword{ spatial }
\keyword{ ts }
