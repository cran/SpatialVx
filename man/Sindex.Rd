\name{Sindex}
\alias{Sindex}
\alias{Sindex.default}
\alias{Sindex.SpatialVx}
\title{
Shape Index
}
\description{
Calculate the shape index (Sindex) as described in AghaKouchak et al. (2011)
}
\usage{
Sindex(x, thresh = NULL, ...)

\method{Sindex}{default}(x, thresh = NULL, ...,
    loc = NULL)

\method{Sindex}{SpatialVx}(x, thresh = NULL, ...,
    time.point = 1, obs = 1, model = 1)
}
\arguments{
    \item{x}{
Default: m by n numeric matrix giving the field for which the shape index is to be calculated.

\code{Sindex.SpatialVx}: list object of class \dQuote{SpatialVx}.
}
    \item{thresh}{
numeric giving a threshold under which (and including, i.e., <=) all values are set to zero, and the shape index is calculated for the non-zero (positive-valued) grid-points.
}
    \item{loc}{
(optional) mn by 2 numeric matrix giving the grid point locations.  If NULL, the expanded grid with x=1:m and y=1:n is used.
}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
    \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
    \item{\dots}{Not used.}
}
\details{
The shape index introduced in AghaKouchak et al. (2011) is defined as

	Sindex = Pmin/P,

where for n = the number of positive-valued grid points, Pmin = 4*sqrt(n) if floor(sqrt(n)) = sqrt(n), and Pmin = 2 * floor(2*sqrt(n)+1) otherwise.  P is the permieter of the non-zero grid points.  Range is 0 to 1.  Values closer to 1 indicate shapes that are closer to circular.
}
\value{
numeric with named components:

    \item{Sindex}{the shape index}
    \item{Pmin,P}{the numerator and denominator (perimeter) that make the Sindex.}

For \dQuote{SpatialVx} objects, the routine is applied to both the verification and forecast objects so that a two-row matrix is returned containing the above vectors for each field.
}
\references{
AghaKouchak, A., Nasrohllahi, N., Li, J., Imam, B. and Sorooshian, S. (2011) Geometrical characterization of precipitation patterns.  \emph{J. Hyrdometeorology}, \bold{12}, 274--285, doi:10.1175/2010JHM1298.1.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{Cindex}}, \code{\link{Aindex}}
}
\examples{
# Re-create Fig. 7a from AghaKouchak et al. (2011).
tmp <- matrix(0, 8, 8)
tmp[3,2:4] <- 1
tmp[5,4:6] <- 1
tmp[7,6:7] <- 1
Sindex(tmp)


}
\keyword{ math }
