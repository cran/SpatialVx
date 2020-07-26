\name{compositer}
\alias{compositer}
\alias{compositer.features}
\alias{compositer.combined}
\alias{compositer.matched}
\alias{plot.composited}
\title{
Create Composite Features
}
\description{

After identifying features in a verification set, re-grid them so that their centroids are all the same, and the new grid is as small as possible to completely contain all of the features in the verification set.

}
\usage{
compositer(x, level = 0, verbose = FALSE, ...)

\method{compositer}{features}(x, level = 0, verbose = FALSE, ...)

\method{compositer}{matched}(x, level = 0, verbose = FALSE, ...)

\method{compositer}{combined}(x, level = 0, verbose = FALSE, ...) 

\method{plot}{composited}(x, ..., type = c("all", "X", "Xhat", "X|Xhat", "Xhat|X"),
    dist.crit = 100, FUN = "mean", col = c("gray", tim.colors(64)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{compositer}: an object of class \dQuote{features}, \dQuote{matched}, or \dQuote{combined}.

\code{plot}: an object of class \dQuote{composited}.
}
  \item{type}{character, stating which composite features should be plotted.  Default makes a two by two panel of plots with all of the choices.}
  \item{dist.crit}{maximum value beyond which any minimum centroid distances are considered too far for features to be \dQuote{present} in the area.}
  \item{FUN}{name of a function to be applied to the composites in order to give the distributional summary.}
  \item{col}{color pallette to be used.}
  \item{level}{numeric used in shrinking the grid to a smaller size.}
  \item{verbose}{
logical, should information be printed to the screen?
}
  \item{\dots}{
Not used by \code{compositer}.  

\code{plot}: Optional arguments to \code{image.plot} (only for the scale legend), such as \code{legend.only}, \code{legend.lab}, \code{legend.mar}, \code{horizontal}, etc.}


}
\details{

This is functionality for performing an analysis similar to the composite verification method of Nachamkin (2004).  See also Nachamkin et al. (2005) and Nachamkin (2009).  The main difference is that this function centers all features to the same point, then re-sizes the grid to the smallest possible size to contain all features.  The "existence" of a feature at the same time point is determined by the centroid distance (because, here, the compositing is done for a large field rather than a small area), but it does not allow for having half of the feature in the domain in order to be considered.

\code{compositer} takes an object of class \dQuote{features} or \dQuote{matched} and centers all of the identified features onto the same point so that all features have the same centroid.  It also then re-grids the composited features so that they are contained on the smallest possible domain that includes all of the features in the verification set.

Generally, because the composite approach is distributional in nature, it makes sense to look at features across multiple time points.  The function \code{combiner} allows for combining features from more than one object of class \dQuote{features} or \dQuote{matched} in order to subsequently run with \code{compositer}.

\code{plot} takes the composite features and adds them together creating a density of the composite features, then, Depending on the \code{type} argument, the verification (\code{type} = \dQuote{X}), model (\code{type} = \dQuote{Xhat}), verification conditioned on the model (\code{type} = \dQuote{X|Xhat}), or the model conditioned on the verification composite features are plotted.  In the case of \code{type} = \dQuote{all}, then a panel of four plots are made with all of these choices.  In the case of the conditional plots, the sum of composites for one field are masked out so that only the density of the other field is plotted where composited features from the first field exist.

}
\value{
A list object of class \dQuote{composited} is returned with all of the same components and attributes as the x argument, but with additional components:
  \item{distances}{List with components X and Xhat giving the minimum centroid distances from each feature in X (Xhat) to a feature in the other field (used for determining the conditional distributions; i.e., a feature is present if its centroid distance is less than some pre-specified amount).}
  \item{Xcentered, Ycentered }{list of \dQuote{owin} objects containing each feature similar to X.feats and Y.feats, but centered on the same spot and re-gridded}
}
\references{
Nachamkin, J. E. (2004) Mesoscale verification using meteorological composites. \emph{Mon. Wea. Rev.}, \bold{132}, 941--955.

Nachamkin, J. E. (2009) Application of the Composite Method to the Spatial Forecast Verification Methods Intercomparison Dataset. \emph{Wea. Forecasting}, \bold{24} (5), 1390--1400, DOI: 10.1175/2009WAF2222225.1.

Nachamkin, J. E., Chen, S. and Schmidt, J. S. (2005) Evaluation of heavy precipitation forecasts using composite-based methods: A distributions-oriented approach. \emph{Mon. Wea. Rev.}, \bold{133}, 2163--2177.
}
\author{
Eric Gilleland
}
\note{
centroids are rounded to the nearest whole number so that interpolation is not necessary.  This may introduce a slight bias in results, but it should not be a major issue.
}


\seealso{
Identifying features: \code{\link{FeatureFinder}}
}
\examples{
x <- y <- matrix(0, 100, 100)
x[2:3,c(3:6, 8:10)] <- 1
y[c(4:7, 9:10),c(7:9, 11:12)] <- 1

x[30:50,45:65] <- 1
y[c(22:24, 99:100),c(50:52, 99:100)] <- 1

hold <- make.SpatialVx( x, y, field.type = "contrived", units = "none",
    data.name = "Example", obs.name = "x", model.name = "y" )

look <- FeatureFinder(hold, smoothpar=0.5)

look2 <- compositer(look)
plot(look2, horizontal = TRUE)

}
\keyword{ manip }
\keyword{ hplot }
