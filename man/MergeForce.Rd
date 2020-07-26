\name{MergeForce}
\alias{MergeForce}
\title{
Force Merges in Matched Feature Objects
}
\description{
Force merges in matched feature objects so that, among other things, subsequent analyses are quicker and cleaner.
}
\usage{
MergeForce(x, verbose = FALSE)
}
\arguments{
  \item{x}{
list object of class \dQuote{matched}.
}
  \item{verbose}{
logical, should progress information be printed to the screen.
}
}
\details{
Objects returned by functions such as \code{deltamm} and \code{centmatch} provide information necessary to merge and match features from \dQuote{features} objects.  In the case of \code{centmatch}, only implicit merges are given, and this function creates objects where the implicit merges are forced to be merged.  In the case of \code{deltamm}, a second pass through might yield better merges/matches in that without a second pass, only features in one field or the other can be merged and matched (not both simultaneously).  Using this function, and apssing the result back through \code{deltamm} can result in subsequent matches of merged features from both fields simultaneously.  Moreover, in some cases, it may be more computationally efficient to run this function once for subsequent analyses/plotting.
}
\value{
A list object of class \dQuote{matched} is returned containing several components and the same attributes as x.

    \item{match.message}{A character string stating how features were matched with (merged) apended.}

    \item{match.type}{character of length 2 naming the original matching function used and this function to note that the features have been forced to be merged/clustered together.  }

    \item{matches}{two-column matrix with forecast object numbers in the first column and corresponding matched observed features in the second column.  If no matches, this will have value integer(0) for each column giving a matrix with dimension 0 by 2.}

    \item{unmatched}{list with components X and Xhat giving the unmatched object numbers, if any, from the observed and forecast fields, resp.  If none, the value will be integer(0).}

    Note that all of the same list components of x are passed back, except for special information (which is usually no longer relevant) such as Q (deltamm), criteria, criteria.values, centroid.distances (centmatch)

Additionally, merges and/or implicit.merges (centmatch) are not included as they have been merged.

}
\author{
Eric Gilleland
}


\seealso{
For identifying features in a field: \code{\link{FeatureFinder}}

For merging and/or matching features: \code{\link{deltamm}}, \code{\link{centmatch}}, \code{\link{plot.matched}}
}
\examples{
x <- y <- matrix(0, 100, 100)
x[2:3,c(3:6, 8:10)] <- 1
y[c(4:7, 9:10),c(7:9, 11:12)] <- 1

x[30:50,45:65] <- 1
y[c(22:24, 99:100),c(50:52, 99:100)] <- 1

hold <- make.SpatialVx( x, y, field.type="contrived", units="none",
    data.name = "Example", obs.name = "x", model.name = "y" )

look <- FeatureFinder(hold, smoothpar=0.5)

look2 <- centmatch( look )

look2

look2 <- MergeForce( look2 )

look2

# plot( look2 )


}
\keyword{ manip }
