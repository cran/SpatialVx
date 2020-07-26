\name{disjointer}
\alias{disjointer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify Disjoint Sets of Connected Components
}
\description{
Identify disjoint sets of contiguous events in a binary field.  In many areas of research, this function finds connected components.
}
\usage{
disjointer(x, method = "C")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric matrix or other object that \code{as.im} from package \pkg{spatstat} works on.
}
  \item{method}{
Same argument as that in \code{connected} from package \code{spatstat}.
}
}
\details{
\code{disjointer} essentially follows the help file for \code{connected} to produce a list object where each component is an image describing one set of connected components (or blobs).  It is essentially a wrapper function to \code{connected}.  This function is mainly used internally by \code{FeatureFinder} and similar, but could be of use outside such functions.
}
\value{
An unnamed list object where each component is an image describing one set of connected components (or blobs).
}
\references{
Park, J.-M., Looney, C.G. and Chen, H.-C. (2000) Fast connected
     component labeling algorithm using a divide and conquer technique.
     Pages 373-376 in S.Y. Shin (ed) \emph{Computers and Their
     Applications:} Proceedings of the ISCA 15th International
     Conference on Computers and Their Applications, March 29--31, 2000,
     New Orleans, Louisiana USA. ISCA 2000, ISBN 1-880843-32-3.

Rosenfeld, A. and Pfalz, J.L. (1966) Sequential operations in
     digital processing.  \emph{Journal of the Association for Computing
     Machinery} \bold{13} 471--494.

}
\author{
Eric Gilleland
}

\seealso{
\code{\link{connected}}, \code{\link{FeatureFinder}}
}
\examples{
##
## For examples, see FeatureFinder
##
}
\keyword{ manip }
