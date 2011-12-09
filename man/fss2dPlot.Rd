\name{fss2dPlot}
\alias{fss2dPlot}
\alias{upscale2dPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create several graphics for list objects returned from hoods2d
}
\description{
Creates several graphics for list objects returned from hoods2d.  Mostly quilt and matrix plots for displaying results of smoothing fields over different neighborhood lengths and thresholds.
}
\usage{
fss2dPlot(x, ...)
upscale2dPlot(object, args, ...)
}
\arguments{
  \item{x}{list object with components fss, fss.random and fss.uniform.  Effectively, it does the same thing as \code{hoods2dPlot}, but adds the fss.random and fss.uniform horizontal lines to the matrix plot.}
  \item{object}{list object with named components: rmse (numeric vector), ets, ts and bias all matrices whose rows represent neighborhood lengths, and whose columns represent thresholds.}
  \item{args}{list object passed to \code{hoods2dPlot}, see its help file for more details.}
  \item{\dots}{Optional arguments to \code{image} and \code{image.plot} for \code{fss2dPlot}, and optional arguments to \code{hoods2dPlot} for \code{upscale2dPlot}}
}
\details{makes quilt and matrix plots for output from \code{hoods2d}.}
\value{
No value is returned.  A series of plots are created.  It may be useful to use this function in conjunction with \code{pdf} in order to view all of the plots.  See the help file for \code{hoods2dPlot} to plot individual results.
}

\author{Eric Gilleland}

\seealso{
\code{\link{hoods2dPlot}}, \code{\link{matplot}}, \code{\link{image}}, \code{\link{image.plot}}, \code{\link{hoods2d}}, \code{\link{pdf}}
}
\examples{
##
## This is effectively an internal function, so the example is commented out
## in order for R's check to run faster.
##
\dontrun{data(geom001)
data(geom000)
data(ICPg240Locs)
hold <- hoods2dPrep( "geom001", "geom000", thresholds=c(0.01,50.01), levels=c(1, 3, 5, 33, 65), loc=ICPg240Locs, units="in/100")
look <- hoods2d(hold, verbose=TRUE)
plot( look)}
}
\keyword{ hplot }
