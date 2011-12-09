\name{centdist}
\alias{centdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the centroid distance between two owin objects
}
\description{
Find the centroid distance between two owin objects.
}
\usage{
centdist(x, y)
}
\arguments{
  \item{x,y}{objects of class "owin" (package \pkg{spatstat}) containing binary images of features of interest.
}
}
\details{
This is a simple wrapper function that calls centroid.owin from package spatstat for each of x and y (to get their centroids), and then finds the (Euclidean) distance between them.  If (xc1, xc2) is the centroid position for x, and (yc1, yc2) is the centroid position for y, then sqrt((xc1 - yc1)^2 + (xc2 - yc2)^2) is returned.
}
\value{
numeric giving the centroid (Euclidean) distance.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{centroid.owin}}, \code{\link{as.im}}, \code{\link{solutionset}}, \code{\link{FeatureMatchAnalyzer}}, \code{\link{FeatureComps}}
}
\examples{
x <- y <- matrix(0, 10, 12)
x[2:3,c(3:6, 8:10)] <- 1
y[c(4:7, 9:10),c(7:9, 11:12)] <- 1

x <- as.im(x)
x <- solutionset(x>0)
y <- as.im(y)
y <- solutionset(y>0)
centdist(x,y)

}
\keyword{ math }
