\name{Sindex}
\alias{Sindex}
\title{
Shape index
}
\description{
Calculate the shape index (Sindex) as described in AghaKouchak et al. (2011)
}
\usage{
Sindex(x, thresh = NULL, loc = NULL)
}
\arguments{
  \item{x}{
m X n numeric matrix giving the field for which the shape index is to be calculated.
}
  \item{thresh}{
numeric giving a threshold under which (and including, i.e., <=) all values are set to zero, and the shape index is calculated for the non-zero (positive-valued) grid-points.
}
  \item{loc}{
(optional) mn X 2 numeric matrix giving the grid point locations.  If NULL, the expanded grid with x=1:m and y=1:n is used.
}
}
\details{
The shape index introduced in AghaKouchak et al. (2011) is defined as

	Sindex = Pmin/P,

where for n = the number of positive-valued grid points, Pmin = 4*sqrt(n) if floor(sqrt(n)) = sqrt(n), and Pmin = 2*(floor(2*sqrt(n)+1) otherwise.  P is the permieter of the non-zero grid points.  Range is 0 to 1.  Values closer to 1 indicate shapes that are closer to circular.
}
\value{
numeric
}
\references{
AghaKouchak, A., N. Nasrohllahi, J. Li, B. Imam, and S. Sorooshian, 2011: Geometrical characterization of precipitation patterns.  \emph{J. Hyrdometeorology}, \bold{12}, 274--285, DOI: 10.1175/2010JHM1298.1
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
