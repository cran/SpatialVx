\name{fss2dPlot}
\alias{fss2dPlot}
\alias{upscale2dPlot}
\title{
Create Several Graphics for List Objects Returned from hoods2d
}
\description{
Creates several graphics for list objects returned from hoods2d.  Mostly quilt and matrix plots for displaying results of smoothing fields over different neighborhood lengths and thresholds.
}
\usage{
fss2dPlot(x, ..., matplotcol = 1:6, mfrow = c(1, 2), add.text = FALSE)

upscale2dPlot(object, args, ..., type = c("all",
                 "gss", "ts", "bias", "rmse"))
}
\arguments{
  \item{x}{list object with components fss, fss.random and fss.uniform.  Effectively, it does the same thing as \code{hoods2dPlot}, but adds the fss.random and fss.uniform horizontal lines to the matrix plot.}
  \item{object}{list object with named components: rmse (numeric vector), ets, ts and bias all matrices whose rows represent neighborhood lengths, and whose columns represent thresholds.}
  \item{args}{list object passed to \code{hoods2dPlot}, see its help file for more details.}
    \item{mfrow}{mfrow parameter (see help file for \code{par}).  If NULL, then the parameter is not re-set.  }
    \item{add.text}{logical, if TRUE, FSS values will be added to the quilt plot as text (in addition to the color).}
    \item{type}{character string stating which plots to make (default is \dQuote{all}).}
  \item{\dots}{Optional arguments to \code{image} and \code{image.plot} for \code{fss2dPlot}, and optional arguments to \code{hoods2dPlot} for \code{upscale2dPlot}}
  \item{matplotcol}{col argument to function \code{matplot}.}
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
\dontrun{
data( "geom001" )
data( "geom000" )
data( "ICPg240Locs" )

hold <- make.SpatialVx( geom000, geom001, thresholds = c(0.01,50.01),
    loc = ICPg240Locs, map = TRUE, projection = TRUE, loc.byrow = TRUE,
    units = "mm/h", data.name = "Geometric", obs.name = "observation",
    model.name = "case 1" )

look <- hoods2d(hold, levels=c(1, 3, 5, 33, 65),
    verbose=TRUE)
plot( look)
}
}
\keyword{ hplot }
