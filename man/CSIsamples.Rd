\name{CSIsamples}
\alias{CSIsamples}
\alias{CSIsamples.default}
\alias{CSIsamples.SpatialVx}
\alias{summary.CSIsamples}
\alias{plot.CSIsamples}
\alias{plot.summary.CSIsamples}
\alias{print.CSIsamples}
\title{
Forecast Verification with Cluster Analysis: The Variation
}
\description{
A variation on cluster analysis for forecast verification as proposed by Marzban and Sandgathe (2008).
}
\usage{
CSIsamples(x, ...)

\method{CSIsamples}{default}(x, ..., xhat, nbr.csi.samples = 100, threshold = 20, 
    k = 100, width = 25, stand = TRUE, z.mult = 0, hit.threshold = 0.1, 
    max.csi.clust = 100, diss.metric = "euclidean", linkage.method = "average", 
    verbose = FALSE)

\method{CSIsamples}{SpatialVx}(x, ..., time.point = 1, obs = 1, model = 1, nbr.csi.samples = 100, 
    threshold = 20, k = 100, width = 25, stand = TRUE, z.mult = 0, 
    hit.threshold = 0.1, max.csi.clust = 100, diss.metric = "euclidean", 
    linkage.method = "average", verbose = FALSE)

\method{summary}{CSIsamples}(object, ...)

\method{plot}{CSIsamples}(x, ...)

\method{plot}{summary.CSIsamples}(x, ...)

\method{print}{CSIsamples}(x, ...)
}
\arguments{
  \item{x,xhat}{default method: matrices giving the verification and forecast fields, resp.

    \dQuote{SpatialVx} method: \code{x} is an object of class \dQuote{SpatialVx}.

    \code{plot}, \code{print} methods:  list object of class \dQuote{CSIsamples} or \dQuote{summary.CSIsamples} (in the case of \code{plot}).
}
  \item{object}{list object of class \dQuote{CSIsamples}.}
  \item{nbr.csi.samples}{integer giving the number of samples to take at each level of the CA.}
  \item{threshold}{numeric giving a value over which is to be considered an event.}
  \item{k}{numeric giving the value for \code{centers} in the call to \code{kmeans}.}
  \item{width}{numeric giving the size of the samples for each cluster sample.}
  \item{stand}{logical, should the data first be standardized before applying CA?}
  \item{z.mult}{numeric giving a value by which to multiply the z- component.  If zero, then the CA is performed on locations only.  Can be used to give more or less weight to the actual values at these locations.}
  \item{hit.threshold}{numeric between zero and one giving the threshold for the proportion of a cluster that is from the verification field vs the forecast field used for determining whether the cluster consitutes a hit (vs false alarm or miss depending).}
  \item{max.csi.clust}{integer giving the maximum number of clusters allowed.}
  \item{diss.metric}{character giving which \code{method} to use in the call to \code{dist} (which dissimilarity metric should be used?).}
  \item{linkage.method}{character giving the name of a linkage method acceptable to the \code{method} argument from the \code{hclust} function of package \pkg{fastcluster}.}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
  \item{verbose}{logical, should progress information be printed to the screen?}
  \item{\dots}{
Not used by \code{CSIsamples} method functions.

\code{summary} method function: the argument \code{silent} may be specified, which is a logical stating whether to print the information to the screen (FALSE) or not (TRUE).  If not given, the summary information will be printed to the screen.

Not used by the \code{plot} method function.}
}
\details{
This function carries out the procedure described in Marzban and Sandgathe (2008) for verifying forecasts.  Effectively, it combines the verification and forecast fields (keeping track of which values belong to which field) and applies CA to the combined field.  Clusters identified with a proportion of values belonging to the verification field within a certain range (defined by the hit.threshold argument) are determined to be hits, misses or false alarms.  From this information, the CSI (at each number of clusters; scale) is calculated.  A sampling scheme is used to speed up the process.

The \code{plot} and \code{summary} functions all give the same information, but in different formats: i.e., CSI by number of clusters (scale).
}
\value{
A list is returned by CSIsamples with components:
\item{data.name}{character vector giving the names of the verification and forecast fields analyzed, resp.}
\item{call }{an object of class \dQuote{call} giving the function call.}
\item{results}{max.csi.clust by nbr.csi.samples matrix giving the caluclated CSI for each sample and iteration of CA.}

The summary method function invisibly returns the same list, but with the additional component:
\item{csi}{vector of length max.csi.clust giving the sample average CSI for each iteration of CA.}

The plot method functions do not return anything.  Plots are created.
}
\references{
Marzban, C., Sandgathe, S. (2008) Cluster Analysis for Object-Oriented Verification of Fields: A Variation. \emph{Mon. Wea. Rev.}, \bold{136}, (3), 1013--1025.
}
\author{
Hillary Lyons, h.lyons \dQuote{at} comcast.net, and modified by Eric Gilleland
}
\note{
Special thanks to Caren Marzban, marzban \dQuote{at} u.washington.edu, for making the CSIsamples (originally called csi.samples) function available for use with this package.
}

\seealso{
\code{\link[fastcluster]{hclust}}, \code{\link[stats]{hclust}}, \code{\link{kmeans}}, \code{\link{clusterer}}
}
\examples{
\dontrun{
grid<- list( x= seq( 0,5,,100), y= seq(0,5,,100))
obj<-Exp.image.cov( grid=grid, theta=.5, setup=TRUE)
look<- sim.rf( obj)
look2 <- sim.rf( obj)

res <- CSIsamples(x=look, xhat=look2, 10, threshold=0, k=100,
                  width=2, z.mult=0, hit.threshold=0.25, max.csi.clust=75)
plot(res)
y <- summary(res)
plot(y)
}
\dontrun{
data( "UKfcst6" )
data( "UKobs6" )
data( "UKloc" )

hold <- make.SpatialVx(UKobs6, UKfcst6, thresholds=0,
    loc=UKloc, map=TRUE, field.type="Rainfall", units="mm/h",
    data.name = "Nimrod", obs.name = "obs 6", model.name = "fcst 6" )

res <- CSIsamples( hold, threshold = 0, k = 200, z.mult = 0.3, hit.threshold = 0.2,
                  max.csi.clust = 150, verbose = TRUE)
plot( res )

summary( res )

y <- summary( res )

plot( y )
}
}
\keyword{ cluster }
\keyword{ hplot }
