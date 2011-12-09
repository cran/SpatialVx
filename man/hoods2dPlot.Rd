\name{hoods2dPlot}
\alias{hoods2dPlot}
\title{Quilt Plot and a Matrix Plot.}
\description{Function to make a quilt plot and a matrix plot for a matrix whose rows represent neighborhood lengths, and whose columns represent different threshold choices.}
\usage{
hoods2dPlot(x, args, ...)
}
\arguments{
  \item{x}{l by q numeric matrix.}
  \item{args}{list object with components: threshold (numeric vector giving the threshold values), qs (optional numeric vector giving the quantiles used if the thresholds represent quantiles rather than hard values), levels (numeric giving the neighborhood lengths (in grid squares) used, units (optional character giving the units for the thresholds)}
  \item{\dots}{optional arguments to \code{image} and \code{image.plot} functions.  May not include xaxt, yaxt, lab, lab, col, or legend.only}
}
\details{Used by \code{plot.hoods2d}, but can be useful for other functions.  Generally, however, this is an internal function that should not be called by the user.  However, it might be called instead of \code{plot.hoods2d} in order to make a subset of the available plots.}
\value{No value is returned.  A plot is created.}
\references{
Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}.  Wadsworth & Brooks/Cole.
}
\author{Eric Gilleland}

\seealso{
 \code{\link{matplot}}, \code{\link{image}}, \code{\link{image.plot}}, \code{\link{hoods2dPrep}}, \code{\link{plot.hoods2d}}, \code{\link{hoods2d}}
}
\examples{
x <- y <- matrix( 0, 50, 50)
x[ sample(1:50,10), sample(1:50,10)] <- rexp( 100, 0.25)
y[ sample(1:50,20), sample(1:50,20)] <- rexp( 400)
hold <- hoods2dPrep("y", "x", thresholds=c(0.1, 0.5), levels=c(1, 3, 20))
look <- hoods2d( hold, which.methods=c("multi.event", "fss"))
hoods2dPlot( look$multi.event$hk, args=hold, main="Hanssen Kuipers Score (Multi-Event Cont. Table)")
\dontrun{
data(geom001)
data(geom000)
data(ICPg240Locs)
hold <- hoods2dPrep( "geom001", "geom000", thresholds=c(0.01,50.01),
    levels=c(1, 3, 5, 17, 33, 65),
    loc=ICPg240Locs, units="in/100")
look <- hoods2d(hold, verbose=TRUE)
hoods2dPlot(look$pragmatic$bss, args=hold)
}
}
\keyword{ hplot }
