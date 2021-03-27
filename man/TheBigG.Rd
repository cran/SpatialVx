\name{TheBigG}
\alias{TheBigG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Spatial Alignment Summary Measure Called G
}
\description{
	The spatial alignment summary measure, G, is a summary comparison for two gridded binary fields.
}
\usage{
TheBigG(X, Xhat, threshold, rule = ">", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X,Xhat}{
m by n matrices giving the \dQuote{observed} and forecast fields, respectively.
}
  \item{threshold,rule}{
The threshold and rule arguments to the \code{binarizer} function.
}
  \item{\dots}{
	  Not used.
}
}
\details{
	This function is an alternative version of Gbeta that does not require the user to select a parameter.  It is not informative about rare events relative to the domain size.  It is the cubed root of the product of two terms.  If A is the set of one-valued grid points in the binary version of \code{X} and B those for \code{Xhat}, then the first term is the size of the symmetric difference between A and B (i.e., an area with grid points squared as the units) and the second term is MED(A,B) * nB with MED(B,A) * nA, where MED is the mean-error distance and nA, nB are the numbers of grid points in each of A and B, respectively.  The second term has units of grid squares so that the product is units of grid squares cubed; hence, the reason for taking the cubed root for G.  The units for G are grid squares with zero being a perfect score and increasing scores imply worsening matches between the sets A and B.  See Gilleland (2021) for more details.
}
\value{
	An object of class \dQuote{TheBigG} is returned.  It is a single number giving the value of G but also has a list of attributes that can be accessed using the \code{attributes} function.  This list includes:
	\item{components}{A vector giving: nA, nB, nAB (number of points in the intersection), number of points in the symmetric difference, MED(A,B), MED(B,A), MED(A,B) * nB, MED(B,A) * nA, followed by the asymmetric versions of G for G(A,B) and G(B,A).}
	\item{threshold}{If a threshold is provided, then this component gives the threshold and rule arguments used.}
}
\references{
Gilleland, E. (2020) Novel measures for summarizing high-resolution forecast performance. \emph{Advances in Statistical Climatology, Meteorology and Oceanography}, \bold{7} (1), 13--34, doi: 10.5194/ascmo-7-13-2021.
}
\author{
	Eric Gilleland
}

\seealso{
\code{\link{Gbeta}}
}
\examples{
data( "obs0601" )
data( "wrf4ncar0531" )
res <- TheBigG( X = obs0601, Xhat = wrf4ncar0531, threshold = 2.1 )
res
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
