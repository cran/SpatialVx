\name{Cindex}
\alias{Cindex}
\alias{Cindex.default}
\alias{Cindex.SpatialVx}
\title{
Connectivity Index
}
\description{
Calculate the connectivity index of an image.
}
\usage{
Cindex(x, thresh = NULL, connect.method = "C", ...)

\method{Cindex}{default}(x, thresh = NULL, connect.method = "C", ...)

\method{Cindex}{SpatialVx}(x, thresh = NULL, connect.method = "C", ...,
    time.point = 1, obs = 1, model = 1)
}
\arguments{
    \item{x}{
Default: m by n numeric matrix giving the field for which the shape index is to be calculated.

\code{Sindex.SpatialVx}: list object of class \dQuote{SpatialVx}.
}
    \item{thresh}{
Set values under (strictly less than) this threshold to zero, and calculate the connectivity index for the resulting image.  If NULL, no threshold is applied.
}
    \item{connect.method}{
character string giving the \code{method} argument for the \code{connected} function of package \pkg{spatstat}.  This must be one of \dQuote{C} or \dQuote{interpreted}.  See the help file for \code{connected} for more details.
}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
    \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
    \item{\dots}{Not used.}
}
\details{
The connectivity index is introduced in AghaKouchak et al. (2011), and is designed to automatically determine how connected an image is.  It is defined by

Cindex = 1 - (NC - 1)/(sqrt(NP) + NC),

where 0 <= Cindex <= 1 is the connectivity index (values close to zero are less connected, and values close to 1 are more connected), NP is the number of nonzero pixels, and NC is the number of isolated clusters.  

The function \code{connected} from package \pkg{spatstat} is used to identify the number of isolated clusters.
}
\value{
numeric giving the connectivity index.
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
data( "geom000" )
Cindex( geom000 )

}
\keyword{ math }
