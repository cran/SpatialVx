\name{Cindex}
\alias{Cindex}
\title{
Connectivity Index
}
\description{
Calculate the connectivity index of an image.
}
\usage{
Cindex(x, thresh = NULL, connect.method = "C")
}
\arguments{
  \item{x}{
m by n numeric matrix for which the connectivity index is to be calculated.
}
  \item{thresh}{
Set values under (strictly less than) this threshold to zero, and calculate the connectivity index for the resulting image.  If NULL, no threshold is applied.
}
  \item{connect.method}{
character string giving the \code{method} argument for the \code{connected} function of package \pkg{spatstat}.  This must be one of \dQuote{C} or \dQuote{interpreted}.  See the help file for \code{connected} for more details.
}
}
\details{
The connectivity index is introduced in AghaKouchak et al. (2011), and is designed to automaticlaly determine how connected an image is.  It is defined by

Cindex = 1 - (NC - 1)/(sqrt(NP) + NC),

where 0 <= Cindex <= 1 is the connectivity index (values close to zero are less connected, and values close to 1 are more connected), NP is the number of nonzero pixels, and NC is the number of isolated clusters.  

The function \code{connected} from package \pkg{spatstat} is used to identify the number of isolated clusters.
}
\value{
single numeric giving the connectivity index.
}
\references{
AghaKouchak, A., Nasrollahi, N., Li, J., Imam, B. and Sorooshian, S. (2011) Geometrical characterization of precipitation patterns.  \emph{J. Hydrometerology}, \bold{12}, 274--285, doi:10.1175/2010JHM1298.1.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{connected}}, \code{\link{as.im}}, \code{\link{Sindex}}, \code{\link{Aindex}}
}
\examples{
data(geom000)
Cindex(geom000)

data(pert000)
Cindex(pert000)

}
\keyword{ math }
