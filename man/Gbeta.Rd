\name{Gbeta}
\alias{Gbeta}
\alias{GbetaIL}
\alias{G2IL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Spatial-Alignment Summary Measures
}
\description{
	Calculates the spatial alignment summary measures from Gilleland (2020)
}
\usage{
Gbeta(X, Xhat, threshold, beta, alpha = 0, rule = ">", ...)

GbetaIL(X, Xhat, threshold, beta, alpha = 0, rule = ">", w = 0.5, ...)

G2IL(X, Xhat, threshold, beta, alpha = 0, rule = ">", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X,Xhat}{ Observed and Forecast fields in the form of a matrix.  \code{X} and \code{Xhat} must have the same dimensions. May be class \dQuote{owin} objects for \code{Gbeta}.}
  \item{threshold}{
	  If \code{X} and \code{Xhat} are not binary \dQuote{owin} objects, then the threshold is used to define the binary field.  See \code{binarizer} for more information about this argument.
}
  \item{beta}{
	  single numeric defining the upper limit for y = y1y2 or y = y1y2(1+y3) determining the rate of decrease in the measure.  Default is half the domain size squared.
}
  \item{alpha}{
	  single numeric defining a lower limit for y = y1y2 or y = y1y2(1+y3) determining what constitutes a perfect match.  The default of zero requires the two fields to be identical, and is probably what is wanted in most cases.
}
  \item{rule}{
	  See \code{binarizer} for more information about this argument.
}
  \item{w}{
	  single numeric between zero and one describing how much weight to give to the first term in the definition of GbetaIL.  See details section below.
}
  \item{\dots}{
	  Not used.
}
}
\details{
	These summary measures were proposed in Gilleland (2020) and provide an index between zero (bad score) and one (perfect score) describing the closeness in spatial alignment between the two fields.  GbetaIL and G2IL also incorporate a distributional summary of the intensity errors.  Gbeta applied between two fields A and B is defined as

	Gbeta(A,B) = max( 1 - y/beta, 0),

	where y = y1 * y2, with y1 a measure of the overlap between A and B (if they overlap completely, y1 = 0); it is the number of points in AB^c and A^cB, where ^c denotes set complement. The term y2 = (MED(A,B) * nB + MED(B,A) * nA), where MED is the mean-error distance (see \code{locperf} for more about MED), with nA and nB representing the number of 1-valued grid points in the sets A and B, resp.  If alpha != 0, then the term y/beta is replaced with (y - alpha) / (beta - alpha).

	GbetaIL is defined to be

	GbetaIL(A,B) = w * Gbeta + (1-w) * theta(A,B),

	where theta is the maximum of zero and the linear correlation coefficient between the intensity values after having sorted them; this part is carried out via a call to \code{qqplot}.  If the number of points in the two sets differs (i.e., if nA != nB), then the larger set is linear interpolated to be the same size as the smaller set.  If both fields are empty, theta = 1.  If field A is empty, then theta = 1 - (nB / N), where N is the size of the domain.  Similarly, if field B is empty.  The rationale is that if one field is empty and the other has very few nonzero points, then the two fields are more similar than if the other field has many points.

	G2betaIL(A,B) = max{1 - (y1 * y2 *(1 + y3 ))/beta, 0},

	where y1 and y2 are as above and y3 is the mean-absolute difference between the sorted values from the sets A and B analogous as for GbetaIL.  If nA = nB = 0, then y3 = 0.  If only one of nA or nB is zero, then the absolute value of the maximum intensity of the other set is used.  The rationale being that this value represents the most egregious error so that if it is large, then the difference is penalized more.

	See Gilleland (2020) for more details about these measures.
}
\value{
	An object of class \dQuote{Gbeta} giving a single numeric giving the value of the summary (index) measures is returned, but with additional attributes that can be obtained using the attributes function, and are also displayed via the print function.  To remove the attributes, the function c can be used.  For Gbeta these include:
\item{components }{numeric vector giving nA, nB, nAB (the number of points in both sets), nA + nB - 2nAB (i.e., the overlap measure, y1), medAB, medBA, medAB * nnB, medBA * nA, and the two asymmetric versions of this measure (see Gilleland, 2020). }
\item{beta }{The value of beta used.}
\item{alpha}{The value of alpha used.}
\item{threshold}{ numeric vector giving the value of the threshold and the rule used.}
\item{data.name}{character vector giving the object names used for X and Xhat.}

The attributes for GbetaIL and G2IL are the same, but the components vector also includes the value of theta.  GbetaIL has an additional attribute called weights that gives w and 1 - w.
}
\references{
	Gilleland, E., 2020. Novel measures for summarizing high-resolution forecast performance. Submitted to Mon. Wea. Rev. on 31 March 2020. Re-submitted on 18 May 2020. Re-submitted to Advances in Statistical Climatology, Meteorology and Oceanography on 19 July 2020.  Temporarily available at: \url{https://ral.ucar.edu/staff/ericg/Gilleland2020.pdf}.
}
\author{
	Eric Gilleland
}


\seealso{
\code{\link{qqplot}}, \code{\link{binarizer}}, \code{\link{solutionset}}, \code{\link{locperf}}
}
\examples{
data( "obs0601" )
data( "wrf4ncar0531" )
res <- Gbeta( X = obs0601, Xhat = wrf4ncar0531, threshold = 2.1 )
c( res )
res
attributes( res )

GbetaIL( X = obs0601, Xhat = wrf4ncar0531, threshold = 2.1,
	beta = 601 * 501 )

G2IL( X = obs0601, Xhat = wrf4ncar0531, threshold = 2.1,
     beta = 601 * 501 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }% use one of  RShowDoc("KEYWORDS")
