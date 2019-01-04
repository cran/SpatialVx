\name{structurogram}
\alias{structurogram}
\alias{plot.structurogram}
\title{
Structure Function for Non-Gridded Spatial Fields.
}
\description{
Computes pairwise differences (raised to the q-th power) as a function of distance.  Returns either raw values or statistics from binning.
}
\usage{
structurogram(loc, y, q = 2, id = NULL, d = NULL, lon.lat = FALSE, dmax = NULL,
    N = NULL, breaks = NULL)

\method{plot}{structurogram}(x, ...)
}
\arguments{
  \item{loc}{
numeric matrix where each row is the coordinate of a point in the field.
}
  \item{x}{list object returned by \code{structurogram} function.}
  \item{y}{
numeric vector giving the value of the field at each location.
}
  \item{q}{
numeric giving the value to which the paired differences should be raised.  Default (q=2) gives the usual semivariogram.
}
  \item{id}{
A 2 column matrix that specifies which variogram differnces to find.  If omitted all possible pairings are found.  This can used if the data has an additional covariate that determines proximity, for example a time window.
}
  \item{d}{
numeric matrix giving the distances among pairs (indexed by \code{id}).  If not included, these are determined directly from \code{loc}.
}
  \item{lon.lat}{
logical, are the coordinates longitude/latitude coordinates?  If so, distances are found using great-circle distance.
}
  \item{dmax}{
numeric giving the maximum distance for which to compute the structure function.
}
  \item{N}{
numeric giving the number of bins to use.
}
  \item{breaks}{
numeric vector giving bin boundaries for binning structure function values.  Need not be equally spaced, but must be ordered.
}
  \item{\dots}{optional arguments to plot function.}
}

\details{
This function is basically an exact copy of \code{vgram} from package \pkg{fields} whereby the differences are raised to a power of q instead of 2.  That is, it calculates the structure function given by Eq (4) in harris et al. (2001).  Namely,

        S_q(l_x,l_y) = <|R(x+l_x,y+l_y) - R(x,y)|^q>

where R is the field of interest, <> denotes the average over pixels in the image (note, in Harris et al. (2001), this is only over non-zero pixels, so is only equivalent to this equation if zero-valued points are first removed from y and loc), l_x and l_y are lags in the x and y directions, resp.  If q=2, then this is the semivariogram.

The \code{plot} method function plots the structure by separation distance (circles) along with a dark blue line giving the bin centers.
}
\value{
A list object of class \dQuote{structurogram} is returned with components:
  \item{d }{numeric vector giving the pair-wise distances.}
  \item{val }{numeric vector giving the structure function values for each distance.}
  \item{q}{numeric giving the value of q passed into the function.}
  \item{call}{Calling string}
  \item{stats}{Matrix of statistics for values in each bin.  Rows are the summaries returned by the stats function or describe (see package fields).  If either breaks or N arguments are not supplied then this component is not computed.}
   \item{centers}{numeric vector giving the bin centers.}
The plot method function does not return anything.
}
\references{
Harris, D., Foufoula-Georgiou, E., Droegemeier, K. K. and Levit, J. J. (2001)  Multiscale statistical properties of a high-resolution precipitation forecast.  \emph{J. Hydrometeorol.}, \bold{2}, 406--418.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{vgram}}, \code{\link{vgram.matrix}}, \code{\link{structurogram.matrix}}
}
\examples{
data( ozone2)
good<- !is.na(ozone2$y[16,])
x<- ozone2$lon.lat[good,] 
y<- ozone2$y[16,good]
     
look <- structurogram( x,y, N=15, lon.lat=TRUE)
plot(look)
# Compare above with results from example for function vgram from package fields.

look <- structurogram( x,y, N=15, lon.lat=TRUE, q=1)
plot(look)

}
\keyword{ math }
