\name{griddedVgram}
\alias{griddedVgram}
\alias{plot.griddedVgram}
\title{
Variograms for a Gridded Verification Set
}
\description{
Find (and plot) variograms for each field in a gridded verification set.
}
\usage{
griddedVgram(object, zero.in = TRUE, zero.out = TRUE, time.point = 1, 
    obs = 1, model = 1, ...)

\method{plot}{griddedVgram}(x, ... )
}
\arguments{
  \item{object}{
list object of class \dQuote{SpatialVx} containing information on the verification set.
}
   \item{zero.in,zero.out}{logical, should the variogram be calculated over the entire field (zero.in), and/or over only the non-zero values (zero.out)?}
   \item{x}{list object as returned by \code{griddedVgram}.}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
  \item{\dots}{
In the case of \code{griddedVgram}, these are optional arguments to the \code{vgram.matrix} function from package \code{fields}.  In the case of \code{plot.griddedVgram}, these are optional arguments to \code{plot.vgram.matrix}, which in turn are optional arguments to \code{image.plot}.
}
}
\details{
Here, the terms semi-variogram and variogram are used interchangeably.

This is a simple wrapper function to \code{vgram.matrix} (entire field) from \pkg{fields} and/or \code{variogram.matrix} (non-zero grid points only) for finding the variogram between two gridded fields.  It calls this function for each of two fields in a verification set.  This function allows one to do the diagnostic analysis proposed in Marzban and Sangathe (2009).
}
\value{
A list object containing the entire list passed in by the object argument, and components:
  \item{Vx.cgram.matrix, Fcst.vgram.matrix }{list objects as returned by vgram.matrix containing the variogram information for each field.}

No value is returned by plot.griddedVgram, plots are created showing the empirical variogram (circles), along with directional empirical variograms (dots), and the variogram by direction (image plot).
}
\references{
Marzban, C. and Sandgathe, S. (2009) Verification with variograms.  \emph{Wea. Forecasting}, \bold{24} (4), 1102--1120, doi:10.1175/2009WAF2222122.1.

}
\author{
Eric Gilleland
}

\seealso{
\code{\link{vgram.matrix}}, \code{\link{make.SpatialVx}}}

\examples{
data( "ExampleSpatialVxSet" )

x <- ExampleSpatialVxSet$vx
xhat <- ExampleSpatialVxSet$fcst

hold <- make.SpatialVx( x, xhat, field.type = "contrived",
    units="none", data.name = "Example", obs.name = "x",
    model.name = "xhat" )

res <- griddedVgram( hold, R = 8 )
plot( res )

}
\keyword{ math }
