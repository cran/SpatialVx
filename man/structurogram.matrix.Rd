\name{structurogram.matrix}
\alias{structurogram.matrix}
\alias{plot.structurogram.matrix}
\title{
Structure Function for Gridded Fields
}
\description{
Calculates the structure function to the q-th order for gridded fields.
}
\usage{
structurogram.matrix(dat, q = 2, R = 5, dx = 1, dy = 1, zero.out = FALSE)

\method{plot}{structurogram.matrix}(x, ...)
}
\arguments{
  \item{dat}{
n by m matrix of numeric values defining a gridded spatial field (or image) such that distances can be determined from their positions in the matrix.
}
  \item{x}{list object output from \code{structurogram.matrix}}
  \item{q}{
numeric giving the order for the structure function (q = 2 yields the more common semi-variogram).
}
  \item{R}{
numeric giving the maximum radius for finding the structure differences assuming that the grid points are spaced one unit apart.  Default is to go to a radius of 5.
}
  \item{dx,dy}{
numeric giving the spacing of the grid points on the x- (y-) axis.  This is used to calculate the correct distance between grid points.
}

  \item{zero.out}{
logical, should zero-valued pixels be ignored?
}
  \item{\dots}{optional arguments to the \code{plot} function.}
}
\details{
This function is basically an exact copy of \code{variogram.matrix}, which itself is a copy of \code{vgram.matrix} from package \pkg{fields} (but allows and ignores missing values, in order to ignore zero-valued pixels and does not include Cressie's robust version of the variogram), whereby the differences are raised to a power of q instead of 2.  That is, it calculates the structure function given by Eq (4) in harris et al. (2001).  Namely,

	S_q(l_x,l_y) = <|R(x+l_x,y+l_y) - R(x,y)|^q>

where R is the field of interest, <> denotes the average over pixels in the image (note, in Harris et al. (2001), this is only over non-zero pixels, so is only equivalent to this equation if zero.out=TRUE), l_x and l_y are lags in the x and y directions, resp.  If q=2, then this is the semivariogram.

The \code{plot} method function makes two plots.  The first shows the structure by separation distance ignoring direction (circles) and all values (i.e., for each direction, dots).  The second shows the structure function values for separation distance and direction (see, e.g., \code{plot.vgram.matrix}).
}
\value{
 A list with the following components:
  \item{d }{numeric vector of distances for the differences (ignoring direction).}
  \item{vgram }{numeric vector giving the structure function values.  Note that the term 'vgram' is used here for compatibility with the plot.vgram.matrix function, which is employed by the plot method function used here.  This set of values ignores direction.}
  \item{d.full}{numeric vector of distances for all possible shifts up distance R.}
  \item{ind}{two column matrix giving the x- and y- increment used to compute shifts.}
  \item{vgram.full}{numeric vector giving the structure function for each direction in addition to separation distance.  Again, the word 'vgram' is used for compatibility with plot.vgram.matrix.}
Note that the plot method function does not return anything.
}
\references{
Harris, D., Foufoula-Georgiou, E., Droegemeier, K. K. and Levit, J. J. (2001)  Multiscale statistical properties of a high-resolution precipitation forecast.  \emph{J. Hydrometeorol.}, \bold{2}, 406--418.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{vgram.matrix}}, \code{\link{vgram}}, \code{\link{structurogram}}, \code{\link{plot.vgram.matrix}}
}
\examples{
data( "lennon" )
look <- structurogram.matrix(lennon, q=2)
plot(look)
# Compare the above with
\dontrun{
look2 <- vgram.matrix(lennon)
dev.new()
par(mfrow=c(1,2),bg="beige")
plot(look2$d, look2$vgram, xlab="separation distance", ylab="variogram")
points(look2$d.full, look2$vgram.full, pch=".")
plot.vgram.matrix(look2)

look <- structurogram.matrix(lennon, q=1)
plot(look)

look <- structurogram.matrix(lennon, q=1, zero.out=TRUE)
plot(look)
}
}
\keyword{ math }
