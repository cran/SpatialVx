\name{saller}
\alias{saller}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to analyze features in a field (image)
}
\description{
Functions to analyze features in a field (image) in conjunction, e.g., with the FeatureSuite function.
}
\usage{
saller(x, object, y = NULL, matches = NULL, d = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
For \code{saller} \code{x} is a list object returned by \code{convthresh} or other feature identification function that returns a list with components \code{X.feats}, \code{Y.feats} (themselves lists with owin class objects defining separate features in the verification and forecast fields, resp.), and \code{X.labeled}, \code{Y.labeled} (fields with the numbers from 0 to the number of features also defining the separate feature locations (e.g., as returned by the \code{connected} function of package \pkg{spatstat}.
}
  \item{object}{
A list object of class "FeatureSuitePrep" as returned by \code{FeatureSuitePrep} providing information about the original fields.
}
  \item{y}{
(optional) if not NULL, it is used instead of \code{x}, and otherwise has the same stucture thereof.  When used with \code{FeatureSuite}, one might want to merge objects in one or both fields before calling the analysis function.  If so, \code{y} is the means by which \code{FeatureSuite} can ensure that this happens.
}
  \item{matches}{
Not used by \code{saller} function.  If matching of features between fields is required, this will be a two-column matrix describing which features from the forecast (first column) match to those of the verificaiton field (second column).
}
  \item{d}{
(optional) the SAL (\code{saller}) method requires division by the longest distance between two border points.  If NULL, this is taken to be simply the length of the longest side.
}
}
\details{
saller: Computes S, A, and L of the SAL method introduced by Wernli et al. (2008).
}
\value{
saller returns a list with components:
  \item{A }{numeric giving the amplitude component.}
  \item{L }{numeric giving the lcoation component.}
  \item{S}{numeric giving the structure component.}
}
\references{
Wernli, H., M. Paulat, M. Hagen, and C. Frei, 2008: SAL--A novel quality measure for the verification of quantitative precipitation forecasts.  \emph{Mon. Wea. Rev.}, \bold{136}, 4470--4487, doi: 10.1175/2008MWR2415.1.
}
\author{
Eric Gilleland
}
\note{
There are several ways to identify features, and some are provided by this package, but only a few.  For example, the method for identifying features in the SAL method as introduced by Wernli et al. (2008) utilizes information from a contour field of a particular variable, and is therefore not currently included in this package.  Users are encouraged to write their own such functions, and should feel free to contribute them to this package by contacing the maintainer.

The SAL method typically looks at a small domain, and it is up to the user to set this up before calling these functions, as they are not designed to handle such a situation.
}

\seealso{
\code{\link{FeatureSuite}}, \code{\link{centroid.owin}}, \code{\link{connected}}, \code{\link{tiles}}, \code{\link{tess}}, \code{\link{deltamm}}
}
\examples{
grid<- list( x= seq( 0,5,,100), y= seq(0,5,,100))
obj <- stationary.image.cov( grid=grid, V=cbind(c(0.25,0.15),c(0.45,0.75)), setup=TRUE)
look<- sim.rf(obj)
look2 <- sim.rf(obj)
look[look<0] <- 0
look2[look2<0] <- 0
q <- mean( c(c(look[look>0]),c(look2[look2>0])), na.rm=TRUE)

hold <- FeatureSuitePrep("look2", "look")
hold2 <- convthresh(hold, smoothpar=5, thresh=q)

if(!is.null(hold2$X.feats) & !is.null(hold2$Y.feats)) {
   N <- max(c(length(hold2$X.feats), length(hold2$Y.feats)), na.rm=TRUE)
   par(mfrow=c(1,2))
   image.plot(hold2$X.labeled, col=c("grey", rainbow(N)), zlim=c(0,N))
   image.plot(hold2$Y.labeled, col=c("grey", rainbow(N)), zlim=c(0,N))

   saller(hold2, hold)   
} # else try a lower threshold.

\dontrun{
data(pert000)
data(pert004)
hold <- FeatureSuitePrep("pert004", "pert000")
look <- convthresh( hold, smoothpar=10.5)
par( mfrow=c(1,2))
zl <- range(c(c(look$X.labeled),c(look$Y.labeled)),finite=TRUE)
image.plot(look$X.labeled, zlim=zl)
image.plot(look$Y.labeled, zlim=zl)

look2 <- deltamm( look, verbose=TRUE)
image.plot(look2$X.labeled)
image.plot(look2$Y.labeled)

   }

}
\keyword{ math }
