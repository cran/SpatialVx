\name{metrV}
\alias{metrV}
\title{
Calculate the binary location metric proposed in Zhu et al. (2011)
}
\description{
Calculate the metric metrV proposed in Zhu et al (2011), which is a linear combination
of the square root of the sum of squared error between two binary fields, and the
mean error distance (Peli and Malah, 1982); or the difference in mean error distances
between two forecast fields and the verificaiton field, if the comparison is performed
between two forecast models against the same verification field.
}
\usage{
metrV(object1, object2 = NULL, lam1 = 0.5, lam2 = 0.5, distfun = "distmapfun", verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object1}{
list object as returned from the \code{locmeasures2dPrep} function.
}
  \item{object2}{
(optional) list object returned from the \code{locmeasures2dPrep} function containing information for a second forecast model to be compared against the first forecast model (object1) against the same verification field, taken from object1.
}
  \item{lam1}{
numeric giving the weight to be applied to the square root of the sum of squared errors of binary fields term in metrV.
}
  \item{lam2}{
numeric giving the weight to be applied to the mean error distance term in metrV.
}
  \item{distfun}{character naming a function with which to calculate the shortest distances between each point x in the grid and the set of events.  Default is the Euclidean distance metric (see the help file for \code{locperf} for more information). }
  \item{verbose}{
logical, should progress information be printed ot the screen.
}
  \item{\dots}{
Optional arguments to the \code{distfun} function.
}
}
\details{
The binary location metric proposed in Zhu et al. (2011) is a linear combination of two measures: the amount of overlap between events in two fields, given by distOV (simply the square root of sum of squared errors between two binary fields), and (if there are events in both fields) the mean error distance described in Peli and Malah (1982); see also Baddeley (1992).  The metric can be computed between a forecast field, M1, and the verificaiton field, V, or it can be compared between two foreast models M1 and M2 with reference to V.  That is,

metrV(M1,M2) = lam1*distOV(I.M1,I.M2) + lam2*distDV(I.M1,I.M2),

where I.M1 (I.M2) is the binary field determined by M1 >= threshold (M2 >= threshold), distOV(I.M1,I.M2) = sqrt( sum( (I.M1 - I.M2)^2)), distDV(I.M1,I.M2) = abs(distob(I.V,I.M1) - distob(I.V,I.M2)), where distob(A,B) is the mean error distance between A and B, given by:

e(A,B) = 1/(N(A))*sqrt( sum( d(x,B)), where the summation is over all the points x corresponding to events in A, and d(x,B) is the minimum of the shortest distance from the point x to each point in B.  e(A,B) is calculated by using the distance transform as calculated by the \code{distmap} function from package \code{spatstat} for computational efficiency.

Note that if there are no events in both fields, then by definition, the term distob(A,B) = 0, and if there are no events in one and only one of the two fields, then a large constant (here, the maximum dimension of the field), is returned.  In this way, distob differs from the mean error distance described in Peli and Malah (1982).

If comparing between the verification field and one forecast model, then the distDV term simplifies to just distob(I.V,I.M1).

One final note is that Eq (6) that defines distOV in Zhu et al. (2011) is correct (or rather, what is used in the paper).  It is not, as is stated below Eq (6) in Zhu et al. (2011) the root *mean* square error, but rather the root square error.  This is function computes Eq (6) as written.
}
\value{
list object of class "metrV" with components:

\item{prep.obj1}{character giving the name of the \code{locmeasures2dPrep} object as in the argument object1.}
\item{prep.obj2}{If supplied, character giving the name of the \code{locmeasures2dPrep} object as in the argument object2.}
\item{OvsM1}{k X 3 matrix whose rows represent thresholds and columns give the component distOV, distob and metrV between the verification field and the forecast model 1.}
\item{OvsM2}{If object2 supplied, k X 3 matrix whose rows represent thresholds and columns give the component distOV, distob and metrV between the verification field and the forecast model 2.}
\item{M1vsM2}{If object2 supplied, k X 3 matrix whose rows represent thresholds and columns give the component distOV, distob and metrV between model 1 and model 2.}
}

\references{
Baddeley, A. J., 1992: An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, W. Forstner and S. Ruwiedel, Eds., Wichmann, 59--78.  Available at: http://www.google.com/url?sa=t&rct=j&q=an%20error%20metric%20for%20binary%20images&source=web&cd=1&ved=0CB4QFjAA&url=http%3A%2F%2Fciteseerx.ist.psu.edu%2Fviewdoc%2Fdownload%3Fdoi%3D10.1.1.52.3879%26rep%3Drep1%26type%3Dpdf&ei=GdXoTrSRJaSIsQLh1tiYCQ&usg=AFQjCNE2Ngehh_J2dXIsonPCMAvAL-fz_g&cad=rja

Peli, T. and D. Malah, 1982: A study on edge detection algorithms.  \emph{Computer Graphics and Image Processing}, \bold{20}, 1--21.

Zhu, M., V. Lakshmanan, P. Zhang, Y. Hong, K. Cheng, and S. Chen, 2011: Spatial verification using a true metric.  \emph{Atmos. Res.}, \bold{102}, 408--419, doi:10.1016/j.atmosres.2011.09.004.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{distob}}, \code{\link{distmap}}, \code{\link{im}}, \code{\link{solutionset}}, \code{\link{deltametric}}, \code{\link{locmeasures2d}}, \code{\link{locmeasures2dPrep}}
}
\examples{
A <- B <- B2 <- matrix( 0, 10, 12)
A[2,3] <- 3
B[4,7] <- 400
B2[10,12] <- 17
hold <- locmeasures2dPrep("A", "B", thresholds=c(0.1,3.1,500))
metrV(hold)

hold2 <- locmeasures2dPrep("A", "B2", thresholds=c(0.1,3.1,500))
metrV( hold, hold2)

\dontrun{
data( pert000)
data(pert001)
testobj <- locmeasures2dPrep( "pert001", "pert000", thresholds=1e-8)
metrV( testobj) # compare to results in Fig. 3 (top right panel) of Zhu et al. (2011).

data( geom000)
data( geom001)
testobj <- locmeasures2dPrep( "geom001", "geom000", thresholds=0)
metrV( testobj)
# compare above to results in Fig. 2 (top right panel) of Zhu et al. (2011).
# Note that they differ wildly.  Perhaps because an actual elliptical area
# is taken in the paper instead of finding the values from the fields themselves?
}

}
\keyword{ math }
