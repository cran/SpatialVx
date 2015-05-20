\name{minboundmatch}
\alias{minboundmatch}
\title{
Minimum Boundary Separation Feature Matching
}
\description{
Match identified features within a spatial verification set via their minimum boundary separation.
}
\usage{
minboundmatch(x, type = c("single", "multiple"), mindist = Inf, verbose = FALSE, ...)
}
\arguments{
  \item{x}{
An object of class \dQuote{features}.
}
  \item{type}{
character string stating either \dQuote{single} or \dQuote{multiple}.  In the former case, each feature in one field will be matched to only one feature in the other, which will be taken to be the features who have the smallest minimum boundary separation.  In the case of \dQuote{multiple}, the \code{mindist} argument should be set to something small enough so that not every feature will be matched to every other feature.  Also, the \code{MergeForce} function may be useful in this case.
}
  \item{mindist}{
single numeric giving the minimum boundary separation distance (measured by grid squares) beyond which features should not be matched.
}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{\dots}{
Optional arguments to the \code{distmap} function from package \pkg{spatstat}.
}
}
\details{
the minimum boundary separation is calculated by first finding the distance map for every feature in the observed field, masking it by each feature in the forecast field, and then finding the minimum of the resulting masked distance map.  If \code{type} is \dQuote{single}, then the features are matched by the smallest minimum boundary separation per feature in each field.  If \code{type} is \dQuote{multiple}, then every feature is matched so long as their minimum boundary separation (measured in grid squares) is less than or equal to \code{mindist}.
}
\value{
A list object of class \dQuote{matched} is returned.  If the type argument is \dQuote{multiple}, then an implicite.merges component is included, which will work with the MergeForce function.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{deltamm}}, \code{\link{centmatch}}, \code{\link{MergeForce}}
}
\examples{
# see ?craer for examples.
}
\keyword{ manip }
\keyword{ math }
