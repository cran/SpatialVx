\name{deltamm}
\alias{deltamm}
\alias{plot.deltamm}
\alias{centmatch}
\title{
Merge and/or match identified features within two fields
}
\description{
Merge and/or match identified features within two fields using the delta metric method described in Gilleland et al. (2008), or the mathcing method of Davis et al. (2006a).
}
\usage{
deltamm(x, y = NULL, object = NULL, max.delta = Inf, verbose = FALSE, ...)
\method{plot}{deltamm}(x, ...)
centmatch(x, y = NULL, object = NULL, criteria = 1, const = 14, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y}{
For \code{deltamm}, a list object with components \code{X.feats} and \code{Y.feats}, each of which are list objects containing numbered components within which are objects of class "owin" containing logical matrices that define objects for the forecast (Y) and verification (X) fields, resp.  For example, as returned by the \code{convthresh} function.  For \code{plot.deltamm} a list object as returned by \code{deltamm}.  Argument \code{y} is used if it is not NULL, otherwise argument \code{x} is used (but only one of \code{x} or \code{y} is used).
}
  \item{object}{
Not used.  Provided for compatibility with the \code{FeatureSutie} function.
}
  \item{max.delta}{single numeric giving a cut-off value for delta that disallows two objects to be merged or matched if the delta between them is larger than this value.}
   \item{criteria}{1, 2 or 3 telling which criteria for determining a match based on centroid distance, D, to use.  The first (1) is a match if D is less than the sum of the sizes of the two features in question (size is the square root of the area of the feature).  The second is a match if D is less than the average size of the two features in question.  The third is a match if D is less than a constant given by the argument \code{const}.}
   \item{const}{numeric giving the number of grid squares whereby if the centroid distance (D) is less than this value, a match is declared (only used if \code{criteria} is 3.}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{\dots}{
For \code{deltamm}: additional optional arguments to the \code{deltametric} function from package \pkg{spatstat}.  For \code{plot.deltamm}, additional arguments to the \code{plot} function.
}
}
\details{
Gilleland et al. (2008) describe a method for automatically merging, and simultaneously, matching identified objects within two fields (a verification set).  The method was proposed with the general method for spatial forecast verification introduced by Davis et al. (2006 a,b) in mind.  It relies heavily on use of a binary image metric introduced by Baddeley (1992a,b) for comparing binary images; henceforth referred to as the delta metric, or just delta.

The procedure is as follows.  Suppose there are m identified forecast features and n identified verification features.

1. Compute delta for each feature identified in the forecast field against each feature identified in the verification field.  Store these values in an m X n matrix, Upsilon.

2. For each of the m rows of Upsilon, rank the values of delta to identify the features, j_1, ..., j_n that provide the lowest (best) to highest (worst) value, and do the same for each of the n columns to find the forecast features i1, ...,i_m that yield the lowest to highest values for each verification feature.

3. Create a new m X n matrix, Psi, whose columns contain delta computed between each of the individual features in the forecast and (first column) the corresponding j_1 feature from the verification field, and each successive column, k, has delta between the i-th forecast feature and the union of j_1, j_2, ..., j_k.

4. Create a similar m X n matrix, Ksi, that has delta computed between each individual feature in the verification field and the successively bigger unions i_1, ..., i_l for the l-th column.

5. Let Q=[Upsilon, Psi, Ksi], and merge and match objects based on the rankings of delta in Q.  That is, find the smallest delta in Q, and determine which mergings (if any) and matchings correspond to this value.  Remove the appropriate row(s) and column(s) of Q corresponding to the already determined matchings and/or mergings.  Repeat this until all features in at least one field have been exhausted.

The above algorithm suffers from two deficiencies.  First, features that are merged in one field cannot be matched to merged features in another field.  One possible remedy for this is to run this algorithm twice, though this is not a universally good solution.  Second, features can be merged and/or matched to features that are very different from each other.  A possible remedy for this is to use the cut-off argument, max.delta, to disallow mergings or matchings between features whose delta value is not <= this cut-off.  In practice, these two deficiencies are not likely very problematic.
}
\value{
A list object is returned containing several components.
\item{X.feats, Y.feats }{list objects with one component for each possibly merged feature for the each of the forecast (Y.feats) and verification (X.feats) fields.  Each component of these lists contains an object of class "owin" containing binary fields indicating the whereabouts of the possibly merged features. }
\item{X.labeled,Y.labeled}{matrices of dimension equal to the original fields identifying the newly merged features of each field.  If l features are matched between fields, then the first l features are labeled in 1 to l such that features labeled r <= l in one field match the corresponding feature labeled r in the other field. }
\item{mm.old.labels}{an l X 2 matrix indicating which features from one field match to features of the other based on the feature labels that were passed into the function (i.e., these labels have now changed).}
\item{mm.new.labels}{list object with components mm (an l X 2 matrix showing which features from one field match to which features of the other; moot because they are labeled identically), and unmatched (list object with components fcst and vx, which identify the newly labeled features in each field that were not matched to features in the other field.}
\item{Q}{an array of dimension n X m X 3 giving all of the delta values that were computed in determining the mergings and matchings.}
}
\references{
Baddeley, A., 1992a.  An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, W. Forstner and S. Ruwiedel, Eds., Wichmann, 59--78.

Baddeley, A., 1992b.  Errors in binary images and an Lp version of the Hausdorff metric.  \emph{Nieuw Arch. Wiskunde}, \bold{10}, 157--183.

Davis CA, BG Brown, and RG Bullock, 2006a. Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.

Davis CA, BG Brown, and RG Bullock, 2006b. Object-based verification of precipitation forecasts, Part II: Application to convective rain systems. \emph{Mon. Wea. Rev.}, \bold{134}, 1785--1795.

Gilleland, E., T. C. M. Lee, J. Halley Gotway, R. G. Bullock, and B. G. Brown, 2008: Computationally efficient spatial forecast verification using Baddeley's delta image metric.  \emph{Mon. Wea. Rev.}, \bold{136}, 1747--1757.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{convthresh}}, \code{\link{disjointer}}, \code{\link{deltametric}}, \code{\link{FeatureSuite}}, \code{\link{owin}}, \code{\link{tess}}, \code{\link{tiles}}, \code{\link{connected}}
}
\examples{
x <- y <- matrix(0, 100, 100)
x[2:3,c(3:6, 8:10)] <- 1
y[c(4:7, 9:10),c(7:9, 11:12)] <- 1

x[30:50,45:65] <- 1
y[c(22:24, 99:100),c(50:52, 99:100)] <- 1

hold <- FeatureSuitePrep("y", "x")
look <- convthresh( hold, smoothpar=0.5) 
par( mfrow=c(1,2))
image.plot( look$X.labeled)
image.plot( look$Y.labeled)

look2 <- deltamm( look)

look3 <- centmatch(look)

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ manip }
\keyword{ graphs}