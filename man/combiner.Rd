\name{combiner}
\alias{combiner}
\title{
Combine Features/Matched Objects
}
\description{
Combine two or more features or matched class objects into one object for aggregation purposes.
}
\usage{
combiner(...)
}
\arguments{
  \item{\dots}{
Two or more objects of class \dQuote{features} or \dQuote{matched} (can also be a list of these objects).
}
}
\details{
Useful for functions such as \code{compositer} and/or (coming soon) aggregating results for feature-based methods.
}
\value{
A list object of class \dQuote{combined} with the same components as the input arguments, but where some components (namely, X, Xhat, X.labeled, Y.labeled) are now arrays containing these values from each combined object.  The lists of lists contained in the X.feats and Y.feats include one long list of lists containing all of the individual features from each object.
}
\author{
Eric Gilleland
}

\seealso{
Functions that create objects of class \dQuote{features} and \dQuote{matched}: \code{\link{FeatureFinder}}, \code{\link{centmatch}}, and \code{\link{deltamm}}.
}
\examples{
# TO DO
}
\keyword{ manip }
