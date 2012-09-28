\name{Aindex}
\alias{Aindex}
\title{
Area index
}
\description{
Calculate Area index described in AghaKouchak et al. (2011).
}
\usage{
Aindex(x, thresh = NULL, dx = 1, dy = 1)
}
\arguments{
  \item{x}{
m X n numeric matrix giving the field for which the area index is to be calculated.
}
  \item{thresh}{
Values under this threshold are set to zero.  If NULL, it will be set to 1e-8 (a very small value).
}
  \item{dx,dy}{
numeric giving the grid point size in each direction if it is desired to apply such a correction.  However, the values are simply canceled out in the index, so these arguments are probably not necessary.  If it is desired to only get the area of the non-zero values in the field, or the convex hull, then these make sense.
}
}
\details{
The area index introduced in AghaKouchak et al. (2011) is given by

	Aindex = A/Aconvex,

where A is the area of the pattern, and Aconvex the area of its convex hull (area.owin from package spatstat is used to calculate this latter area, and the functions as.im and solutionset from spatstat are also used by this function).  Values are between 0 and 1.  Values closer to unity indicate a more structured pattern, and values closer to zero indicate higher dispersiveness of the pattern, but note that two highly structured patterns far away from each other may also give a low value (see examples below).  Because of this, this measure is perhaps best applied to individual features in a field.
}
\value{
A list with components:
\item{Aindex }{numeric giving the area index.}
\item{A,Aconvex }{numeric giving the area of th epattern and the convex hull, resp.}
\item{dx,dy}{the values of dx and dy as input to the function.}
}
\references{
AghaKouchak, A., N. Nasrollahi, J. Li, B. Imam, and S. Sorooshian, 2011: Geometrical characterization of precipitation patterns.  \emph{J. Hydrometeorology}, \bold{12}, 274--285, DOI: 10.1175/2010JHM1298.1
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{as.im}}, \code{\link{solutionset}}, \code{\link{convexhull}}, \code{\link{Cindex}}, \code{\link{Sindex}}
}
\examples{
# Gemetric shape that is highly structured.
data(geom000)
Aindex(geom000)
Aindex(geom000,dx=4,dy=4)

\dontrun{
# Two separate areas with highly structured shapes, but far away from each other.
data(pert000)
Aindex(pert000, thresh=20, dx=4, dy=4)
image(pert000>20)
}
}
\keyword{ math }
