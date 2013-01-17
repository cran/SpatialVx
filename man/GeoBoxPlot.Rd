\name{GeoBoxPlot}
\alias{GeoBoxPlot}
\title{
Geographic Box Plot
}
\description{
Make a geographic box plot as detailed in Willmott et al. (2007).
}
\usage{
GeoBoxPlot(x, areas, ...)
}
\arguments{
  \item{x}{
numeric giving the values to be box-plotted.
}
  \item{areas}{
numeric of same length as x giving the associated areas for each value.
}
  \item{\dots}{
optional arguments to the \code{boxplot} function of R.  The argument plot is not allowed.
}
}
\details{
This function makes the geographic box plots described in Willmott et al. (2007) that calculates the five statistics in such a way as to account for the associated areas (e.g., over a grid where each grid box may have differing areas).

Missing values are not handled, and ideally should be handled before calling ths routine.

In future, this function may allow other options for \code{x} than currently, but for now, only numeric vectors are allowed.
}
\value{
List with the same components as returned by \code{boxplot}.
}
\references{
Willmott, C. J., Robeson, S. M. and Matsuura, K. (2007)  Geographic box plots.  \emph{Physical Geography}, \bold{28}, 331--344, doi:10.2747/0272-3646.28.4.331.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{boxplot}}
}
\examples{
##
## Reproduce the boxplots of Fig. 1 in Willmott et al. (2007).
##
x <- c(4,9,1,3,10,6,7)
a <- c(rep(1,4),2,1,3)
boxplot( x, at=1, xlim=c(0,3))
GeoBoxPlot(x, a, at=2, add=TRUE)
axis( 1, at=c(1,2), labels=c("Traditional", "Geographic"))
}
\keyword{ hplot }
