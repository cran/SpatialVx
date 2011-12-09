\name{kernel2d}
\alias{kernel2d}
\title{Perform a onvolution smooth on a 2-d field.}
\description{Perform a onvolution smooth on a 2-d field.  Default is to take an average over all neighbors within (n-1)/2 grid points from each grid point.  Uses FFT with the convolution theorem for computational efficiency.}
\usage{
kernel2d(x, n, W = NULL, K = NULL, X = NULL, xdim = NULL, Nxy = NULL, setup = FALSE, verbose = FALSE)
}

\arguments{
  \item{x}{k X m matrix to be smoothed.}
  \item{n}{neighborhood length}
  \item{W}{(optional) 2k X 2m possibly complex matrix of scaled Fourier transformed kernel weights.  If NULL, these will be computed, but if passed, it will save one FFT in computation time.  This should not be given if setup is TRUE, or you will not get what you expect.}
  \item{K}{(optional) 2k X 2m matrix defining an alternative kernel to the default, which is a boxcar kernel that multiplies 1/(n^2) times the value at each grid point within (n-1)/2 of each grid point from the original field, and zero everywhere else.  This is an advanced option, and probably should be left NULL for most users.  If a different kernel is desired, you might try \code{image.smooth} from the fields package instead.}
  \item{X}{(optional) 2k X 2m matrix giving the Fourier transformed 'x'.  Can be used to save an FFT in computation, if this has already been calculated.}
  \item{xdim}{(optional) numeric vector of length 2 giving the dimensions of x.  Not really necessary, but as it will have already been calculated, seems silly to have to keep re-calculating it.  If NULL, it will be calculated here.}
  \item{Nxy}{(optional) total number of grid points of x.  Similar to xdim argument, not really necessary, and will be calculated if not passed.}
  \item{setup}{logical, should just the Fourier transformed kernel weights be returned for subsequent calls to this function?}
  \item{verbose}{logical, should progress information be printed to the screen?}
}
\details{This 2-d spatial kernel smoother applies the smoother of Roberts and Lean (2008) to a spatial field (see Hastie and Tibshirani, 1990 sec. 2.6; Ritter and Wilson, 2001, chapter 8; Barrett and Myers, 2004 for details about this type of convolution smoothing).  Specifically, if X is a matrix of grid points, then the returned field, denoted by Ebert (2008) as <X>s, is a smoothed field such that the value at each grid point '<X>s[i,j]' is given by: <X>s[i,j] = sum_k sum_l X[i + k - 1 - (n-1)/2, j + l - 1 - (n-1)/2]*K[i + k - 1 - (n-1)/2, j + l - 1 - (n-1)/2], where n is the neighborhood length, k,l = 1, ..., n, and K[i + k - 1 - (n-1)/2, j + l - 1 - (n-1)/2] (default) is constant, and equal to 1/n^2.  In order to be fast, loops are avoided.  Instead, the convolution theorem is applied with a Fast Fourier Transform (FFT).  If the weights 'W' are supplied, then you will save one FFT in computation time.

The convolution theorem says that the Fourier transform of a convolution between two functions f and g is equal to the product of the Fourier transformed functions.  That is, if F denotes the Fourier transform, and * the convolution operator, F( f*g ) = k F(f)F(g), where 'k' is a scaling factor.  The neighborhood smooth is given by a convolution between the field and a boxcar kernel (i.e., a square around a point with constant value 1/n^2).  Because of the FFT, this enables a fast way to compute this convolution.
 
In order to zero-pad the field, and perform a cyclic convolution, it is necessary to expand the field, 'x', and re-arrange the boxcar kernel (or else it will not be centered on the points).  If zero padding is not desired, then a field that has been extrapolated to an appropriate size in another way should be used, and a subset going back to the original size could be used.  Alternatively, a subset of an appropriate size could be taken from the resulting smoothed field in order to avoid edge effects.  The latter is probably a wise move.

}
\value{If setup is FALSE, then a k X m matrix giving the neighborhood smoothed field is returned.  Otherwise, a 2k X 2m possibly complex matrix giving the Fourier transformed kernel weights are returned, which can be used to save an FFT in computation time for subsequent calls to this function by supplying the W argument with this result.}
\references{
Barrett, H. H. and Myers, K. J., 2004: \emph{Foundations of Image Science}.  Wiley Series in Pure and Applied Optics, Editor: B. E. A. Saleh, Hoboken, New Jersey, 1540pp.

Ebert E. E., 2008. Fuzzy verification of high resolution gridded forecasts: A review and proposed framework. \emph{Meteorol. Appl.}, \bold{15}, 51--64. DOI: 10.1002/met.25.  Available at http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf

Hastie, T. J. and Tibshirani, R. J., 1990: \emph{Generalized Additive Models}. Chapman \& Hall/CRC Monographs on Statistics and Applied Probability 43, 335pp.

Ritter, G. X. and Wilson, J. N., 2001: \emph{Handbook of Computer Vision Algorithms in Image Algebra}.  2nd Edition, CRC Press, Boca Raton, Florida, U.S.A., 417pp.

Roberts, N. M. and H. W. Lean, 2008: Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events.  \emph{Mon. Wea. Rev.}, \bold{136}, 78--97.  DOI: 10.1175/2007MWR2123.1.
  
}
\author{Eric Gilleland}
\note{If n is 1, then the field is returned without applying any smoothing.}
\seealso{
\code{\link{fft}}, \code{\link{image.smooth}}, \code{\link{hoods2d}}, \code{\link{Fourier2d}}
}
\examples{
look <- matrix( 0, 10, 12)
look[4,7] <- 1
kernel2d( look, 3)
# The above returns the shape of the kernel applied, which
# is a square of length 3 centered on the grid point in look
# that has a value of 1.

# What happens if the 1 is on the edge?  the effect is zero padding:
look <- look*0
look[1,1] <- 1
kernel2d( look, 3)

# Suppose we want to do the above for several, say l, neighborhood lengths.
# We can save an FFT for l-1 of the convolutions.
look <- look*0
look[4,7] <- 1
lookFFT <- Fourier2d( look, xdim=c(10,12))
dim( lookFFT) # Note the dimension is twice that of look.
kernel2d( look, 3, X=lookFFT)
kernel2d( look, 5, X=lookFFT)

# Now, suppose we want to apply the same kernel smooth to different fields.
# We can save an FFT for each subsequent calculation as follows.
wg <- kernel2d( look, 3, setup=TRUE)
dim( wg) # Note the dimension is twice that of look.
kernel2d( look, 3, W=wg)
look <- look*0
look[8,5] <- 1
kernel2d( look, 3, W=wg)
look[5, 10] <- 1
kernel2d( look, 3, W=wg)

}
\keyword{ math }
\keyword{ manip }
