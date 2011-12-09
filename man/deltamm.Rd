\name{deltamm}
\alias{deltamm}
\alias{plot.matched}
\alias{centmatch}
\title{
Merge and/or Match Identified Features Within Two Fields
}
\description{
Merge and/or match identified features within two fields using the delta metric method described in Gilleland et al. (2008), or the matching only method of Davis et al. (2006a).
}
\usage{
deltamm(object, x, y = NULL, max.delta = Inf, time.point = 1, 
    model = 1, verbose = FALSE, ...)

\method{plot}{matched}(x, ..., set.pw = FALSE,
    horizontal = TRUE, loc.byrow = TRUE)

centmatch(object, x, y = NULL, criteria = 1, const = 14, time.point = 1, 
    model = 1, verbose = FALSE)
}
\arguments{
    \item{object}{
A list object of class \dQuote{SpatialVx}.  Provided mainly for compatibility with the \code{FeatureSutie} function, but does utilize information from its attributes.  Can be missing.
}
  \item{x,y}{
For \code{deltamm}, a list object with components \code{X.feats} and \code{Y.feats}, each of which are list objects containing numbered components within which are objects of class \dQuote{owin} containing logical matrices that define objects for the forecast (Y) and verification (X) fields, resp.  For example, as returned by the \code{convthresh} function.  For \code{plot.matched} a list object as returned by \code{deltamm}.  Argument \code{y} is used if it is not NULL, otherwise argument \code{x} is used (but only one of \code{x} or \code{y} is used).
}
  \item{max.delta}{single numeric giving a cut-off value for delta that disallows two objects to be merged or matched if the delta between them is larger than this value.}
   \item{criteria}{1, 2 or 3 telling which criteria for determining a match based on centroid distance, D, to use.  The first (1) is a match if D is less than the sum of the sizes of the two features in question (size is the square root of the area of the feature).  The second is a match if D is less than the average size of the two features in question.  The third is a match if D is less than a constant given by the argument \code{const}.}
   \item{const}{numeric giving the number of grid squares whereby if the centroid distance (D) is less than this value, a match is declared (only used if \code{criteria} is 3.}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{model}{numeric indicating which forecast model to select for the analysis.}
    \item{set.pw}{logical, should a panel of plots be determined and set by the function.}
    \item{horizontal}{optional argument to \code{image.plot} from \pkg{fields}.}
   \item{loc.byrow}{logical, only used if field is a projection, this determines how the locations should be put into matrices.}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{\dots}{
For \code{deltamm}: additional optional arguments to the \code{deltametric} function from package \pkg{spatstat}.  For \code{plot.matched}, additional arguments to the \code{plot} function.  May use \code{no.label} (logical, should a title be put on the plot or not) and \code{plot.set} (anything other than NULL will keep this function from re-setting the par option \code{mfrow}; otherwise, it will set it to have one row and two columns).
}
}
\details{
\code{deltamm}:

Gilleland et al. (2008) describe a method for automatically merging, and simultaneously, matching identified objects within two fields (a verification set).  The method was proposed with the general method for spatial forecast verification introduced by Davis et al. (2006 a,b) in mind.  It relies heavily on use of a binary image metric introduced by Baddeley (1992a,b) for comparing binary images; henceforth referred to as the delta metric, or just delta.

The procedure is as follows.  Suppose there are m identified forecast features and n identified verification features.

1. Compute delta for each feature identified in the forecast field against each feature identified in the verification field.  Store these values in an m by n matrix, Upsilon.

2. For each of the m rows of Upsilon, rank the values of delta to identify the features, j_1, ..., j_n that provide the lowest (best) to highest (worst) value, and do the same for each of the n columns to find the forecast features i1, ...,i_m that yield the lowest to highest values for each verification feature.

3. Create a new m by n matrix, Psi, whose columns contain delta computed between each of the individual features in the forecast and (first column) the corresponding j_1 feature from the verification field, and each successive column, k, has delta between the i-th forecast feature and the union of j_1, j_2, ..., j_k.

4. Create a similar m by n matrix, Ksi, that has delta computed between each individual feature in the verification field and the successively bigger unions i_1, ..., i_l for the l-th column.

5. Let Q=[Upsilon, Psi, Ksi], and merge and match objects based on the rankings of delta in Q.  That is, find the smallest delta in Q, and determine which mergings (if any) and matchings correspond to this value.  Remove the appropriate row(s) and column(s) of Q corresponding to the already determined matchings and/or mergings.  Repeat this until all features in at least one field have been exhausted.

The above algorithm suffers from two deficiencies.  First, features that are merged in one field cannot be matched to merged features in another field.  One possible remedy for this is to run this algorithm twice, though this is not a universally good solution.  Second, features can be merged and/or matched to features that are very different from each other.  A possible remedy for this is to use the cut-off argument, max.delta, to disallow mergings or matchings between features whose delta value is not <= this cut-off.  In practice, these two deficiencies are not likely very problematic.

\code{centmatch}:

This function works similarly as \code{deltamm} in that it is designed to work with \code{FeatureSuite} with the same end of matching (but not merging) objects between fields.  It is based on the method proposed by Davis et al. (2006a).
}
\value{
A list object of class \dQuote{matched} is returned by both centmatch and deltamm containing several components.
\item{X.feats, Y.feats }{list objects with one component for each possibly merged feature for the each of the forecast (Y.feats) and verification (X.feats) fields.  Each component of these lists contains an object of class \dQuote{owin} containing binary fields indicating the whereabouts of the possibly merged features. }
\item{X.labeled,Y.labeled}{matrices of dimension equal to the original fields identifying the newly merged features of each field.  If l features are matched between fields, then the first l features are labeled in 1 to l such that features labeled r <= l in one field match the corresponding feature labeled r in the other field. }
\item{mm.old.labels}{an l by 2 matrix indicating which features from one field match to features of the other based on the feature labels that were passed into the function (i.e., these labels have now changed).}
\item{mm.new.labels}{list object with components mm (an l by 2 matrix showing which features from one field match to which features of the other; moot because they are labeled identically), and unmatched (list object with components fcst and vx, which identify the newly labeled features in each field that were not matched to features in the other field.}
\item{Q}{(deltamm only) an array of dimension n by m by 3 giving all of the delta values that were computed in determining the mergings and matchings.}

If the argument \sQuote{object} is passed in, then the list object will also contain nearly the same attributes, with the data.name attribute possibly changed to reflect the specific model used.  It will also contain a time.point and model attribute.

}
\references{
Baddeley, A. (1992a)  An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, W. Forstner and S. Ruwiedel, Eds., Wichmann, 59--78.

Baddeley, A. (1992b)  Errors in binary images and an Lp version of the Hausdorff metric.  \emph{Nieuw Arch. Wiskunde}, \bold{10}, 157--183.

Davis, C. A., Brown, B. G. and Bullock, R. G. (2006a) Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.

Davis, C. A., Brown, B. G. and Bullock, R. G. (2006b) Object-based verification of precipitation forecasts, Part II: Application to convective rain systems. \emph{Mon. Wea. Rev.}, \bold{134}, 1785--1795.

Gilleland, E., Lee, T. C. M.,  Halley Gotway, J., Bullock, R. G. and Brown, B. G. (2008) Computationally efficient spatial forecast verification using Baddeley's delta image metric.  \emph{Mon. Wea. Rev.}, \bold{136}, 1747--1757.
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

hold <- make.SpatialVx(x, y, field.type="contrived", units="none",
    data.name=c("Example", "x", "y"))

look <- convthresh(hold, smoothpar=0.5) 
par( mfrow=c(1,2))
image.plot(look$X.labeled)
image.plot(look$Y.labeled)

look2 <- deltamm(hold, x=look)

look3 <- centmatch(hold, x=look)

\dontrun{
data(pert000)
data(pert004)
data(ICPg240Locs)

hold <- make.SpatialVx(pert000, pert004,
    loc=ICPg240Locs, projection=TRUE, map=TRUE,
    field.type="Precipitation", units="mm/h",
    data.name=c("ICP Perturbed Cases", "pert000", "pert004"))

look <- convthresh(hold, smoothpar=10.5)
plot(look, set.pw=TRUE)

look2 <- deltamm(hold, x=look, verbose=TRUE)
plot(look2, set.pw=TRUE)
   }

}
\keyword{ math }
\keyword{ manip }
\keyword{ graphs}
