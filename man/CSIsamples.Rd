\name{CSIsamples}
\alias{CSIsamples}
\alias{summary.CSIsamples}
\alias{plot.CSIsamples}
\alias{plot.summary.CSIsamples}
\title{
Forecast Verification with Cluster Analysis: The Variation
}
\description{
A variation on cluster analysis for forecast verification as proposed by Marzban and Sandgathe (2008).
}
\usage{
CSIsamples(X, Y, nbr.csi.samples = 100, threshold = 20, k = 100, width = 25, stand = TRUE,
	z.mult = 0, hit.threshold = 0.1, max.csi.clust = 100, diss.metric = "euclidean",
	linkage.method = "average", verbose = FALSE)
\method{summary}{CSIsamples}(object, ...)
\method{plot}{CSIsamples}(x, ...)
\method{plot}{summary.CSIsamples}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X,Y}{matrices giving the verification and forecast fiedls, resp.}
  \item{object,x}{list objects of class "CSIsamples" or "summary.CSIsamples".}
  \item{nbr.csi.samples}{integer giving the number of samples to take at each level of the CA.}
  \item{threshold}{numeric giving a value over which is to be considered an event.}
  \item{k}{numeric giving the value for \code{centers} in the call to \code{kmeans}.}
  \item{width}{numeric giving the size of the samples for each cluster sample.}
  \item{stand}{logical, should the data first be standardized before applying CA?}
  \item{z.mult}{numeric giving a value by which to multiply the z- component.  If zero, then the CA is performed on locations only.  Can be used to give more or less weight to the actual values at these locations.}
  \item{hit.threshold}{numeric between zero and one giving the threshold for the proportion of a cluster that is from the verification field vs the forecast field used for determining whether the cluster consitutes a hit (vs false alarm or miss depending).}
  \item{max.csi.clust}{integer giving the maximum number of clusters allowed.}
  \item{diss.metric}{character giving which \code{method} to use in the call to \code{dist} (which dissimilarity metric should be used?).}
  \item{linkage.method}{character giving the name of a linkage method acceptable to the \code{method} argument from the \code{hclust} function of package \code{fastcluster}.}
  \item{verbose}{logical, should progress information be printed to the screen?}
  \item{\dots}{In the \code{summary} method function, the argument \code{silent} may be specified, which is a logical stating whether to print the information to the screen (FALSE) or not (TRUE).  If not given, the summary information will be printed to the screen.  Not used by the plot method functions.}
}
\details{
This function carries out the procedure described in Marzban and Sandgathe (2008) for verifying forecasts.  Effectively, it combines the verification and forecast fields (keeping track of which values belong to which field) and applies CA to the combined field.  Clusters identified with a proportion of values belonging to the verification field within a certain range (defined by the hit.threshold argument) are determined to be hits, misses or false alarms.  From this information, the CSI (at each number of clusters; scale) is calculated.  A sampling scheme is used to speed up the process.

The plot and summary functions all give the same information, but in different formats: i.e., CSI by number of clusters (scale).
}
\value{
A list is returned by CSIsamples with components:
\item{Vx.name, Fcst.name }{characters giving the names of the verification and forecast fields analyzed, resp.}
\item{call }{an object of class "call" giving the function call.}
\item{results}{max.csi.clust X nbr.csi.samples matrix giving the caluclated CSI for each sample and iteration of CA.}

The summary method function invisibly returns the same list, but with the additional component:
\item{csi}{vector of length max.csi.clust giving the sample average CSI for each iteration of CA.}

The plot method functions do not return anything.  Plots are created.
}
\references{
Marzban, C., S. Sandgathe, 2008: Cluster Analysis for Object-Oriented Verification of Fields: A Variation. \emph{Mon. Wea. Rev.}, \bold{136}, (3), 1013 - 1025.
}
\author{
Hillary Lyons (modified by Eric Gilleland)
}
\note{
Special thanks to Caren Marzban for making the CSIsamples (originally called csi.samples) function available for use with this package.
}

\seealso{
\code{\link[fastcluster]{hclust}}, \code{\link[stats]{hclust}}, \code{\link{kmeans}}, \code{\link{clusterer}}
}
\examples{
grid<- list( x= seq( 0,5,,100), y= seq(0,5,,100))
obj<-Exp.image.cov( grid=grid, theta=.5, setup=TRUE)
look<- sim.rf( obj)
look2 <- sim.rf( obj)

res <- CSIsamples(X=look, Y=look2, 10, threshold=0, k=100,
		width=2, z.mult=0, hit.threshold=0.25, max.csi.clust=75)
plot(res)
y <- summary(res)
plot(y)

\dontrun{
data(UKfcst6)
data(UKobs6)
res <- CSIsamples(UKobs6, UKfcst6, threshold=0, k=200, z.mult=0.3, hit.threshold=0.2,
		max.csi.clust=150, verbose=TRUE)
plot(res)
summary(res)
y <- summary(res)
plot(y)
}
}
\keyword{ cluster }
\keyword{ hplot }
