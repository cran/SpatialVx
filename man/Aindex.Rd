\name{Aindex}
\alias{Aindex}
\alias{Aindex.default}
\alias{Aindex.SpatialVx}
\title{
Area Index
}
\description{
Calculate Area index described in AghaKouchak et al. (2011).
}
\usage{
Aindex(x, thresh = NULL, dx = 1, dy = 1, ...)

\method{Aindex}{default}(x, thresh = NULL, dx = 1, dy = 1, ...)

\method{Aindex}{SpatialVx}(x, thresh = NULL, dx = 1, dy = 1, ...,
    time.point=1, obs = 1, model=1)
}
\arguments{
    \item{x}{
Default: m by n numeric matrix giving the field for which the area index is to be calculated.

\code{Aindex.SpatialVx}: list object of class \dQuote{SpatialVx}.
}
    \item{thresh}{
Values under this threshold are set to zero.  If NULL, it will be set to 1e-8 (a very small value).
}
    \item{dx,dy}{
numeric giving the grid point size in each direction if it is desired to apply such a correction.  However, the values are simply canceled out in the index, so these arguments are probably not necessary.  If it is desired to only get the area of the non-zero values in the field, or the convex hull, then these make sense.
}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
    \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
    \item{\dots}{Not used.}
}
\details{
The area index introduced in AghaKouchak et al. (2011) is given by

	Aindex = A/Aconvex,

where A is the area of the pattern, and Aconvex the area of its convex hull (area.owin from package spatstat is used to calculate this latter area, and the functions as.im and solutionset from spatstat are also used by this function).  Values are between 0 and 1.  Values closer to unity indicate a more structured pattern, and values closer to zero indicate higher dispersiveness of the pattern, but note that two highly structured patterns far away from each other may also give a low value (see examples below).  Because of this property, this measure is perhaps best applied to individual features in a field.
}
\value{
numeric vector (or two-row matrix in the case of \code{Aindex.SpatialVx}) with named components (columns):
\item{Aindex }{numeric giving the area index.}
\item{A,Aconvex }{numeric giving the area of th epattern and the convex hull, resp.}
\item{dx,dy}{the values of dx and dy as input to the function.}
}
\references{
AghaKouchak, A., Nasrollahi, N., Li, J., Imam, B. and Sorooshian, S. (2011) Geometrical characterization of precipitation patterns.  \emph{J. Hydrometeorology}, \bold{12}, 274--285, doi:10.1175/2010JHM1298.1.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{as.im}}, \code{\link{solutionset}}, \code{\link{convexhull}}, \code{\link{Cindex}}, \code{\link{Sindex}}
}
\examples{
# Gemetric shape that is highly structured.
data( "geom000" )

Aindex( geom000 )
Aindex( geom000, dx = 4, dy = 4 )

}
\keyword{ math }
