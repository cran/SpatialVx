\name{craer}
\alias{craer}
\alias{print.craered}
\title{
Contiguous Rain Area
}
\description{
Calculate the contiguous rain area (CRA) for matched features.
}
\usage{
craer(x, type = c("regular", "fast"), rotate = FALSE, loss, loss.args = NULL,
    interp = "bicubic", method = "BFGS", stages = TRUE, verbose = FALSE, ...)

\method{print}{craered}(x, ...)
}
\arguments{
  \item{x}{
\code{craer}: A list object of class \dQuote{matched}.

\code{print}: Object returned by \code{craer}.
}
  \item{type}{
character string that must be one of \dQuote{regular} or \dQuote{fast}.  See details section for more information.
}
  \item{rotate}{
logical, should the features also be rotated into alignement (in addition to their being translated)?
}
  \item{loss}{
character string naming a function to calculate the loss function when using \code{type} \dQuote{regular}.  Default uses \code{QlossRigid}, which is based on an assumption of iid normally distributed errors (only important if you try to obtain CIs for the rigid transformations, which is currently not readily available), but is similar to squared error loss.  See \code{QcorrRigid} for another option based on correlations.  This argument is ignored if \code{type} is \dQuote{fast}.
}
  \item{loss.args}{
list object with optional arguments for \code{loss}.  Only used if \code{type} is \dQuote{fast}.
}
  \item{interp}{
character string naming the interpolation method used in calls to \code{Fint2d}.  Must be one of \dQuote{round}, \dQuote{bilinear} or \dQuote{bicubic}.  Only used if \code{type} is \dQuote{fast}.
}
  \item{method}{
character string naming the numerical optimization method to use (see \code{optim} for choices).  Only used if \code{type} is \dQuote{fast}.
}
  \item{stages}{
logical, if \code{rotate} is TRUE, then should the optimal translation be found first, and then further try to find the optimal translation and rotation using the results from the initial optimization?  The idea is to help the optimization, but not clear that it matters.
}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{\dots}{
Optional arguments to \code{optim} besides \code{method}.  Not used by \code{print}.
}
}
\details{
The contiguous rain area (CRA) spatial verification method was first introduced by Ebert and McBride (2000).  After identifying features and matching them across fields using your favorite feature identification and feature matching functions, the method attempts to find an optimal rigid transformation between the matched features, and then applies the mean square error (MSE) breakdown given by:

MSE( total ) = MSE( displacement ) + MSE( volume ) + MSE( pattern ),

where MSE( displacement ) = MSE( total ) - MSE( shifted ), and MSE( shifted ) is the MSE between the rigidly transformed forecast field and the observed one. MSE( volume )  is the squared bias (differences in means) between the rigidly transformed forecast and the observed field.  MSE( pattern ) is the difference between MSE( shift ) and MSE( volume ).  In the case of rotations with \code{stages} = FALSE, these values are the same, but the displacement also incorporates the rotation.  It is more complicated if this argument is TRUE.  In this case, the values are as with its being FALSE, but some additional MSE breakdowns are given for the translation only situation, but MSE( volume ) and MSE( pattern ) are based on the translation and rotation together.

Finding the optimal rigid transformation between two features can be difficult.  Two options are supplied here.  The first (\code{type} = \dQuote{regular}) attempts to find the optimal rigid transformation according to a given loss function via numerical optimization of the loss function.  The second (\code{type} = \dQuote{fast}) uses the image moments to find the distance between feature centroids (as the optimal translation) and difference between feature major axis angles (as the optimal rotation).

In general, Ebert and McBride (2000) only considered translations (not rotations), and it is unclear whether or not the addition of rotations is worthwhile. 
}
\value{
Returns a matrix of class \dQuote{craered} with named columns specifying the results for each feature match in the order they were supplied via the \dQuote{matched} list object.
}
\references{
Ebert, E. E. and McBride, J. L. (2000) Verification of precipitation in weather systems: determination of systematic errors. \emph{J. Hydrology}, \bold{239}, 179--202.
}
\author{
Eric Gilleland
}
\note{
The rigid transformation is performed by the \code{rigider} function, which is experimental at this stage in the game.
}

\seealso{
\code{\link{optim}}, \code{\link{imomenter}}

Feature identification function: \code{\link{FeatureFinder}}

For some feature matching functions, see:

\code{\link{centmatch}}, \code{\link{deltamm}}, \code{\link{minboundmatch}}

For finding the optimal rigid transformations:

\code{\link{rigider}}, \code{\link{imomenter}}

For rigidly transforming a field: \code{\link{rigidTransform}}

Two-dimensional interpolation (useful because the transformations typically do not map exactly to a grid square):

\code{\link{Fint2d}}
}
\examples{
x <- y <- matrix(0, 100, 100)
x[2:3,c(3:6, 8:10)] <- 1
y[c(4:7, 9:10),c(7:9, 11:12)] <- 1

x[30:50,45:65] <- 1
y[c(22:24, 99:100),c(50:52, 99:100)] <- 1

hold <- make.SpatialVx(x, y, field.type="contrived", units="none",
    data.name=c("Example", "x", "y"))

look <- FeatureFinder(hold, smoothpar=0.5)

look2 <- minboundmatch( look, type = "multiple", mindist = 20 )
look2 <- MergeForce( look2 )

craer( look2, type = "fast", verbose = TRUE)

\dontrun{
data(pert000)
data(pert004)
data(ICPg240Locs)

hold <- make.SpatialVx(pert000, pert004,
    loc=ICPg240Locs, projection=TRUE, map=TRUE, loc.byrow = TRUE,
    field.type="Precipitation", units="mm/h",
    data.name=c("ICP Perturbed Cases", "pert000", "pert004"))

look <- FeatureFinder(hold, smoothpar=10.5, thresh = 5)
plot(look)

look2 <- minboundmatch(look, verbose = TRUE)
plot(look2)

craer( look2 )
   }

}
\keyword{ manip }
\keyword{ math }
