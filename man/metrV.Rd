\name{metrV}
\alias{metrV}
\alias{metrV.default}
\alias{metrV.SpatialVx}
\alias{print.metrV}
\title{
Binary Location Metric Proposed in Zhu et al. (2011)
}
\description{
Calculate the metric metrV proposed in Zhu et al (2011), which is a linear combination
of the square root of the sum of squared error between two binary fields, and the
mean error distance (Peli and Malah, 1982); or the difference in mean error distances
between two forecast fields and the verification field, if the comparison is performed
between two forecast models against the same verification field.
}
\usage{
metrV(x, ...)

\method{metrV}{default}(x, xhat, xhat2 = NULL, thresholds, lam1 = 0.5, lam2 = 0.5, 
    distfun = "distmapfun", a = NULL, verbose = FALSE, ...)

\method{metrV}{SpatialVx}(x, time.point = 1, obs = 1, model = 1, lam1 = 0.5, lam2 = 0.5, 
    distfun = "distmapfun", verbose = FALSE, ...)

\method{print}{metrV}(x, ...)
}
\arguments{
  \item{x}{
Either a list object as returned by \code{make.SpatialVx} or a matrix representing a verificaiton grid.  For the \code{print} method, this is an object returned by \code{metrV}.
}
  \item{xhat, xhat2}{
(xhat2 is optional) matrix representing a forecast grid.
}
    \item{thresholds}{q X 2 or q X 3 (if \code{xhat} is not NULL) matrix giving the thresholds to apply to the verification field (first column) and each forecast field.}
  \item{lam1}{
numeric giving the weight to be applied to the square root of the sum of squared errors of binary fields term in metrV.
}
  \item{lam2}{
numeric giving the weight to be applied to the mean error distance term in metrV.
}
  \item{distfun}{character naming a function with which to calculate the shortest distances between each point x in the grid and the set of events.  Default is the Euclidean distance metric (see the help file for \code{locperf} for more information). }
    \item{a}{list object giving certain information about the verification set.  These are the attributes of the \dQuote{SpatailVx} object.  May be used here to include information (as attributes of the returned object) that would otherwise not be available to the \code{print} method function.  In particular, the components, \code{msg}, \code{data.name} and \code{qs} are printed if available.}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.  May have length one or two.  If it has length two, the second value is taken to be the second forecast model (i.e., \code{xhat2} in the call to \code{metrV.default}).}
  \item{verbose}{
logical, should progress information be printed ot the screen.
}
  \item{\dots}{
Optional arguments to the \code{distfun} function.
}
}
\details{
The binary location metric proposed in Zhu et al. (2011) is a linear combination of two measures: the amount of overlap between events in two fields, given by \code{distOV} (simply the square root of sum of squared errors between two binary fields), and (if there are events in both fields) the mean error distance described in Peli and Malah (1982); see also Baddeley (1992).  The metric can be computed between a forecast field, M1, and the verificaiton field, V, or it can be compared between two foreast models M1 and M2 with reference to V.  That is,

metrV(M1,M2) = lam1*distOV(I.M1,I.M2) + lam2*distDV(I.M1,I.M2),

where I.M1 (I.M2) is the binary field determined by M1 >= threshold (M2 >= threshold), distOV(I.M1,I.M2) = sqrt( sum( (I.M1 - I.M2)^2)), distDV(I.M1,I.M2) = abs(distob(I.V,I.M1) - distob(I.V,I.M2)), where distob(A,B) is the mean error distance between A and B, given by:

e(A,B) = 1/(N(A))*sqrt( sum( d(x,B)), where the summation is over all the points x corresponding to events in A, and d(x,B) is the minimum of the shortest distance from the point x to each point in B.  e(A,B) is calculated by using the distance transform as calculated by the \code{distmap} function from package \code{spatstat} for computational efficiency.

Note that if there are no events in both fields, then by definition, the term distob(A,B) = 0, and if there are no events in one and only one of the two fields, then a large constant (here, the maximum dimension of the field), is returned.  In this way, distob differs from the mean error distance described in Peli and Malah (1982).

If comparing between the verification field and one forecast model, then the distDV term simplifies to just distob(I.V,I.M1).

One final note is that Eq (6) that defines \code{distOV} in Zhu et al. (2011) is correct (or rather, what is used in the paper).  It is not, as is stated below Eq (6) in Zhu et al. (2011) the root *mean* square error, but rather the root square error.  This function computes Eq (6) as written.
}
\value{
list object of class \dQuote{metrV} with components:

\item{OvsM1}{k by 3 matrix whose rows represent thresholds and columns give the component distOV, distob and metrV between the verification field and the forecast model 1.}
\item{OvsM2}{If object2 supplied, k by 3 matrix whose rows represent thresholds and columns give the component distOV, distob and metrV between the verification field and the forecast model 2.}
\item{M1vsM2}{If object2 supplied, k by 3 matrix whose rows represent thresholds and columns give the component distOV, distob and metrV between model 1 and model 2.}

May also contain attributes as passed by either the a argument or the \dQuote{SpatialVx} object.
}

\references{
Baddeley, A. J. (1992) An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, W. Forstner and S. Ruwiedel, Eds., Wichmann, 59--78.

Peli, T. and Malah, D. (1982) A study on edge detection algorithms.  \emph{Computer Graphics and Image Processing}, \bold{20}, 1--21.

Zhu, M., Lakshmanan, V. Zhang, P. Hong, Y. Cheng, K. and Chen, S. (2011) Spatial verification using a true metric.  \emph{Atmos. Res.}, \bold{102}, 408--419, doi:10.1016/j.atmosres.2011.09.004.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{distob}}, \code{\link{distmap}}, \code{\link{im}}, \code{\link{solutionset}}, \code{\link{deltametric}}, \code{\link{locmeasures2d}}, \code{\link{make.SpatialVx}}
}
\examples{

A <- B <- B2 <- matrix( 0, 10, 12)
A[2,3] <- 3
B[4,7] <- 400
B2[10,12] <- 17
hold <- make.SpatialVx( A, list(B, B2), thresholds = c(0.1, 3.1, 500),
    field.type = "contrived", units = "none",
    data.name = "Example", obs.name = "A",
    model.name = c("B", "B2") )

metrV(hold)

metrV(hold, model = c(1,2) )

\dontrun{

data( "geom000" )
data( "geom001" )

testobj <- make.SpatialVx( geom000, geom001, thresholds = 0,
    projection = TRUE, map = TRUE, loc = ICPg240Locs, loc.byrow = TRUE,
    field.type = "Precipitation", units = "mm/h",
    data.name = "ICP Geometric Cases", obs.name = "geom000",
    model.name = "geom001" )

metrV(testobj)

# compare above to results in Fig. 2 (top right panel)
# of Zhu et al. (2011).  Note that they differ wildly.
# Perhaps because an actual elliptical area is taken in
# the paper instead of finding the values from the fields
# themselves?
}

}
\keyword{ math }
