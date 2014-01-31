\name{deltamm}
\alias{deltamm}
\alias{plot.matched}
\alias{print.matched}
\alias{centmatch}
\title{
Merge and/or Match Identified Features Within Two Fields
}
\description{
Merge and/or match identified features within two fields using the delta metric method described in Gilleland et al. (2008), or the matching only method of Davis et al. (2006a).
}
\usage{
deltamm(x, p = 2, max.delta = Inf, const = Inf, verbose = FALSE, ...)

\method{plot}{matched}(x, ..., type = c("both", "obs", "model"))

\method{print}{matched}(x, ...)

centmatch(x, criteria = 1, const = 14, distfun = "rdist", areafac = 1,
    verbose = FALSE, ...)
}
\arguments{
  \item{x}{
For \code{deltamm} and \code{centmatch}, an object of class \dQuote{features} with components \code{X.labeled}, \code{Y.labeled} (matrices with numbered features), as well as, \code{X.feats} and \code{Y.feats}, each of which are list objects containing numbered components within which are objects of class \dQuote{owin} containing logical matrices that define features for the forecast (Y) and verification (X) fields, resp.  For example, as returned by the \code{convthresh} function.  For \code{plot} an object of class \dQuote{matched}.  Argument \code{y} is used if it is not NULL, otherwise argument \code{x} is used (but only one of \code{x} or \code{y} is used).  If \code{x} and \code{y} are missing, but not \code{object}, then \code{object} will be used, in which case it must be of class \dQuote{features}.
}
  \item{p}{Baddeley delta metric parameter.  A value of 1 gives arithmetic averages, Inf gives the Hausdorff metric and -Inf gives a minimum.  The default of 2 is most common.}
  \item{max.delta}{single numeric giving a cut-off value for delta that disallows two features to be merged or matched if the delta between them is larger than this value.}
  \item{const}{\code{deltamm}: a constant value over which the shortest distances in a distance map are set.  See Gilleland (2011) and Schwedler and Baldwin (2011) for more information about this parameter.
  \code{centmatch}: numeric giving the number of grid squares whereby if the centroid distance (D) is less than this value, a match is declared (only used if \code{criteria} is 3.
}
   \item{criteria}{1, 2 or 3 telling which criteria for determining a match based on centroid distance, D, to use.  The first (1) is a match if D is less than the sum of the sizes of the two features in question (size is the square root of the area of the feature).  The second is a match if D is less than the average size of the two features in question.  The third is a match if D is less than a constant given by the argument \code{const}.}
    \item{distfun}{character string naming a distance function.  Default uses \code{rdist} from the \pkg{fields} package.}
    \item{areafac}{single numeric used to multiply by grid-space based area in order to at least approximate the correct distance (e.g., using the ICP test cases, 4 would make the areas approximately square km instead of grid points).  This should not be used unless \code{distfun} is \dQuote{rdist.earth} in which case it will use the locations given in the call to \code{make.SpatialVx}, which are assumed to be lat/lon coordinates.}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{type}{character string naming which plot should be created.}
  \item{\dots}{
For \code{deltamm}: additional optional arguments to the \code{distmap} function from package \pkg{spatstat}.  For \code{centmatch}: optional arguments to \code{distfun}.  For \code{plot.matched}, additional arguments to \code{image.plot} concerning the color legend bar only.
}
}
\details{
\code{deltamm}:

Gilleland et al. (2008) describe a method for automatically merging, and simultaneously, matching identified features within two fields (a verification set).  The method was proposed with the general method for spatial forecast verification introduced by Davis et al. (2006 a,b) in mind.  It relies heavily on use of a binary image metric introduced by Baddeley (1992a,b) for comparing binary images; henceforth referred to as the delta metric, or just delta.

The procedure is as follows.  Suppose there are m identified forecast features and n identified verification features.

1. Compute delta for each feature identified in the forecast field against each feature identified in the verification field.  Store these values in an m by n matrix, Upsilon.

2. For each of the m rows of Upsilon, rank the values of delta to identify the features, j_1, ..., j_n that provide the lowest (best) to highest (worst) value, and do the same for each of the n columns to find the forecast features i1, ...,i_m that yield the lowest to highest values for each verification feature.

3. Create a new m by n matrix, Psi, whose columns contain delta computed between each of the individual features in the forecast and (first column) the corresponding j_1 feature from the verification field, and each successive column, k, has delta between the i-th forecast feature and the union of j_1, j_2, ..., j_k.

4. Create a similar m by n matrix, Ksi, that has delta computed between each individual feature in the verification field and the successively bigger unions i_1, ..., i_l for the l-th column.

5. Let Q=[Upsilon, Psi, Ksi], and merge and match features based on the rankings of delta in Q.  That is, find the smallest delta in Q, and determine which mergings (if any) and matchings correspond to this value.  Remove the appropriate row(s) and column(s) of Q corresponding to the already determined matchings and/or mergings.  Repeat this until all features in at least one field have been exhausted.

The above algorithm suffers from two deficiencies.  First, features that are merged in one field cannot be matched to merged features in another field.  One possible remedy for this is to run this algorithm twice, though this is not a universally good solution.  Second, features can be merged and/or matched to features that are very different from each other.  A possible remedy for this is to use the cut-off argument, max.delta, to disallow mergings or matchings between features whose delta value is not <= this cut-off.  In practice, these two deficiencies are not likely very problematic.

\code{centmatch}:

This function works similarly as \code{deltamm}, though it does not merge features.  It is based on the method proposed by Davis et al. (2006a).  It is possible for more than one object to be matched to the same object in another field.  As a result, when plotting, it might appear that features have been merged, but they have not been.  For informational purposes, the criteria, appelled \code{criteria.values} (as determined by the \code{criteria} argument), along with the centroid distance matrix, appelled \code{centroid.distances}, are returned.

	\code{plot}: The plot method function for matched features plots matched features across fields in the same color using \code{rainbow}.  Unmatched features in either field are all colored gray.  Zero values are colored white.  In the case of features created by \code{centmatch}, multiple features from one field could be matched to the same feature in the other so that, in such cases, it is not possible to decipher which feature in one field is mached to which feature in the other exactly.  The \code{print} method function will tell you which features matched between fields, so one can plot the originally derived features (e.g., from \code{convthresh}) to identify matched features.  On the other hand, one can subjectively determine if same colored features should be merged or not.

}

\value{

A list object of class \dQuote{matched} is returned by both centmatch and deltamm containing several components added to the value of x or y passed in, and possibly with attributes inhereted from object.

    \item{match.message}{A character string stating how features were matched.}

    \item{match.type}{character string naming the matching funciton used.}

    \item{matches}{two-column matrix with forecast object numbers in the first column and corresponding matched observed features in the second column.  If no matches, this will have value integer(0) for each column giving a matrix with dimension 0 by 2.}

    \item{unmatched}{list with components X and Xhat giving the unmatched object numbers, if any, from the observed and forecast fields, resp.  If none, the value will be integer(0).}

    \item{Q}{(deltamm only) an array of dimension n by m by 3 giving all of the delta values that were computed in determining the mergings and matchings.}

    \item{criteria}{(centmatch only) 1, 2, or 3 as given by the criteria argument.}

    \item{criteria.values, centroid.distances}{(centmatch only) matrices giving the forecast by observed object criteria and centroid distances.}

    \item{implicit.merges}{(centmatch only) list displaying multiple matches for each field (this could define potential merges).  Each component of the list is a unified set of matched features in the form of two-column matrices analogous to the matches component.  If there are no implicit mergings or no matched features, this component will be named, but also NULL.  Note: such implicit mergings may or may not make physical sense, and are not considered to be merged generally, but will show up as having been merged/clustered when plotted.}

If the argument \sQuote{object} is passed in, then the list object will also contain nearly the same attributes, with the data.name attribute possibly changed to reflect the specific model used.  It will also contain a time.point and model attribute.

}

\references{

Baddeley, A. (1992a)  An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, W. Forstner and S. Ruwiedel, Eds., Wichmann, 59--78.

Baddeley, A. (1992b)  Errors in binary images and an Lp version of the Hausdorff metric.  \emph{Nieuw Arch. Wiskunde}, \bold{10}, 157--183.

Davis, C. A., Brown, B. G. and Bullock, R. G. (2006a) Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.

Davis, C. A., Brown, B. G. and Bullock, R. G. (2006b) Object-based verification of precipitation forecasts, Part II: Application to convective rain systems. \emph{Mon. Wea. Rev.}, \bold{134}, 1785--1795.

Gilleland, E. (2011) Spatial Forecast Verification: Baddeley's Delta Metric Applied to the ICP Test Cases. \emph{Wea. Forecasting}, \bold{26} (3), 409--415.

Gilleland, E., Lee, T. C. M.,  Halley Gotway, J., Bullock, R. G. and Brown, B. G. (2008) Computationally efficient spatial forecast verification using Baddeley's delta image metric.  \emph{Mon. Wea. Rev.}, \bold{136}, 1747--1757.

Schwedler, B. R. J. and Baldwin, M. E. (2011) Diagnosing the sensitivity of binary image measures to bias, location, and event frequency within a forecast verification framework. \emph{Wea. Forecasting}, \bold{26}, 1032--1044.
}

\author{
Eric Gilleland
}

\seealso{
To identify features: \code{\link{convthresh}}, \code{\link{threshfac}}, \code{\link{threshsizer}}

To force merges (implicit or otherwise): \code{\link{MergeForce}}, 

Other functions used to identify features within the above mentioned functions (all from \pkg{spatstat}):

\code{\link{disjointer}}, \code{\link{deltametric}}, \code{\link{owin}}, \code{\link{tess}}, \code{\link{tiles}}, \code{\link{connected}}
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

look2 <- deltamm(look)

FeatureTable(look2)

look3 <- centmatch(look)

FeatureTable(look3)

\dontrun{
data(pert000)
data(pert004)
data(ICPg240Locs)

hold <- make.SpatialVx(pert000, pert004,
    loc=ICPg240Locs, projection=TRUE, map=TRUE, loc.byrow = TRUE,
    field.type="Precipitation", units="mm/h",
    data.name=c("ICP Perturbed Cases", "pert000", "pert004"))

look <- convthresh(hold, smoothpar=10.5, thresh = 5)
plot(look)

look2 <- deltamm(look, verbose = TRUE)
plot(look2)
   }

}
\keyword{ math }
\keyword{ manip }
\keyword{ graphs}
