\name{centdist}
\alias{centdist}
\title{
Centroid Distance Between Two Identified Objects
}
\description{
Find the centroid distance between two identified objects/features.
}
\usage{

centdist(x, y, distfun = "rdist", loc = NULL, ...)

}

\arguments{

  \item{x,y}{objects of class \dQuote{owin} (package \pkg{spatstat}) containing binary images of features of interest.}
  \item{distfun}{character string naming a distance function that should take arguments \code{x1} and \code{x2} as 1 by 2 matrices, and return a single numeric value.  Default uses the \pkg{fields} function, \code{rdist}, where the fields function \code{rdist.earth} is an obvious alternative.}
  \item{loc}{two-column matrix giving the location values for which to calculate the centroids.  If NULL, indices according to the dimension of the fields are used.}
  \item{\dots}{optional arguments to \code{distfun}.}

}

\details{
This is a simple function that calculates the centroid for each of x and y (to get their centroids), and then finds the distance between them according to \code{distfun}.  The centroids are calculated using \code{FeatureProps}.
}
\value{
numeric giving the centroid distance.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{FeatureProps}}, \code{\link{as.im}}, \code{\link{solutionset}}, \code{\link{FeatureMatchAnalyzer}}, \code{\link{FeatureComps}}
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
