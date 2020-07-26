\name{Mij}
\alias{Mij}
\title{
Raw Image Moments.
}
\description{
Calculate the raw Hu image moment Mij.
}
\usage{
Mij(x, s, i = 0, j = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix.
}
  \item{s}{
A two-column matrix giving the location coordinates.  May be missing in which case they are assumed to be integers giving the row and column numbers.
}
  \item{i,j}{
Integer giving the moment order for each coordinate x and y, resp.
}
}
\details{
The raw moment M_{ij} (Hu 1962) is calculated by

M_{ij} = sum(x^i * y^j * Im[i, j])

where x and y are the pixel coordinates and Im is the (image) matrix.  Various useful properties of an image may be gleaned from certain moments.  For example, the image area is given by M_{00}, and the image centroid is (M_{10} / M_{00}, M_{01} / M_{00}).  The image orientation angle can also be derived.
}
\value{
A single numeric giving the desired moment is returned.
}
\references{
Hu, M. K. (1962) Visual Pattern Recognition by Moment Invariants. \emph{IRE Trans. Info. Theory}, \bold{IT-8}, 179--187.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{imomenter}}
}
\examples{
data( "geom000" )

Mij( geom000 ) # area
}
\keyword{arith}
