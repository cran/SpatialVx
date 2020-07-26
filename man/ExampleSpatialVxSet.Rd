\name{ExampleSpatialVxSet}
\alias{ExampleSpatialVxSet}
\docType{data}
\title{
Simulated Spatial Verification Set
}
\description{
A simulated spatial verification set for use by various examples for this package.
}
\usage{data(ExampleSpatialVxSet)}
\format{
  The format is:
List of 2
 $ vx  : num [1:50, 1:50] 0 0 0 0 0 0 0 0 0 0 ...
 $ fcst: num [1:50, 1:50] 0.0141 0 0 0 0 ...
}
\details{
The data here were generated using the \code{sim.rf} function from \pkg{fields} (Furrer et al., 2012):

x <- y <- matrix(0, 10, 12)
x[2:3,c(3:6, 8:10)] <- 1
y[c(1:2, 9:10),c(3:6)] <- 1

grid <- list(x=seq(0,5,,50), y=seq(0,5,,50))
obj <- Exp.image.cov(grid=grid, theta=0.5, setup=TRUE)
x <- sim.rf(obj)
x[x < 0] <- 0
x <- zapsmall(x)

y <- sim.rf(obj)
y[y < 0] <- 0
y <- zapsmall(y)

}
\references{
Reinhard Furrer, Douglas Nychka and Stephen Sain (2012). fields: Tools for spatial data. R package version 6.6.3.  http://CRAN.R-project.org/package=fields
}
\examples{
data( "ExampleSpatialVxSet" )
x <- ExampleSpatialVxSet$vx
xhat <- ExampleSpatialVxSet$fcst
par(mfrow=c(1,2))
image.plot(x, col=c("gray",tim.colors(64)))
image.plot(xhat, col=c("gray",tim.colors(64)))
}
\keyword{datasets}
