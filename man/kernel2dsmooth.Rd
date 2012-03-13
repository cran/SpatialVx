\name{kernel2dsmooth}
\alias{kernel2dsmooth}
\title{Perform a onvolution smooth on a 2-d field.}
\description{Perform a onvolution smooth on a 2-d field.  Default is to take an average over all neighbors within (n-1)/2 grid points from each grid point.  Uses FFT with the convolution theorem for computational efficiency.}
\usage{
kernel2dsmooth(x, kernel.type=NULL, K = NULL, W = NULL, X = NULL, xdim = NULL, Nxy = NULL, setup = FALSE, verbose = FALSE, ...)
}

\arguments{
  \item{x}{matrix to be smoothed.}
  \item{kernel.type}{(optional) character naming the kernel type accepted by \code{kernel2dmeitsjer}.  One and only one of this argument, K or W must be supplied.}
   \item{K}{(optional) matrix defining a kernel to be applied.  this function will expand, and flip the kernel about its center, so ideally it will have odd dimensions.}
  \item{W}{(optional) possibly complex matrix of scaled Fourier transformed kernel weights.  If NULL, these will be computed, but if passed, it will save one FFT in computation time.  This should not be given if setup is TRUE, or you will not get what you expect.  The dimensions are determined by this function, so it is probably best to supply this matrix only from a previous call with setup=TRUE.  The dimensions are chosen to optimize the FFT calculations (see Details section).}
  \item{X}{(optional) matrix giving the Fourier transformed 'x'.  Can be used to save an FFT in computation, if this has already been calculated.}
  \item{xdim}{(optional) numeric vector of length 2 giving the dimensions of x.  Not really necessary, but as it will have already been calculated, seems silly to have to keep re-calculating it.  If NULL, it will be calculated here.}
  \item{Nxy}{(optional) total number of grid points of x.  Similar to xdim argument, not really necessary, and will be calculated if not passed.}
  \item{setup}{logical, should just the Fourier transformed kernel weights, W, be returned for subsequent calls to this function?}
  \item{verbose}{logical, should progress information be printed to the screen?}
  \item{\dots}{optional arguments to \code{kernel2dmeitsjer} as required by the specific kernel given to kernel.type}
}
\details{This 2-d spatial kernel smoother applies a kernel smoother to a spatial field (see Hastie and Tibshirani, 1990 sec. 2.6; Ritter and Wilson, 2001, chapter 8; Barrett and Myers, 2004 for details about this type of convolution smoothing).  Specifically, if X is a matrix of grid points, then the returned field, denoted by Ebert (2008) as <X>s, is a smoothed field such that the value at each grid point '<X>s[i,j]' is given by: <X>s[i,j] = sum_k sum_l X[i + k - 1, j + l - 1]*K[i, j], where k,l = 1, ..., n, and K[i, j] is the kernel matrix.  In order to be fast, loops are avoided.  Instead, the convolution theorem is applied with a Fast Fourier Transform (FFT).  If the weights 'W' are supplied, then you will save one FFT in computation time.

The convolution theorem says that the Fourier transform of a convolution between two functions f and g is equal to the product of the Fourier transformed functions.  That is, if F denotes the Fourier transform, and * the convolution operator, F( f*g ) = k F(f)F(g), where 'k' is a scaling factor.  The neighborhood smooth is given by a convolution between the field and a boxcar kernel (i.e., a square around a point with constant value 1/n^2).  Because of the FFT, this enables a fast way to compute this convolution.
 
In order to zero-pad the field, and perform a cyclic convolution, it is necessary to expand the field, 'x', and re-arrange the kernel (or else it will not be centered on the points).  If zero padding is not desired, then a field that has been extrapolated to an appropriate size in another way should be used, and a subset going back to the original size could be used.  Alternatively, a subset of an appropriate size could be taken from the resulting smoothed field in order to avoid edge effects.  The latter is probably a wise move.  The image is expanded to the nearest power of two above the dimension of N=x + dimension of K - 1 in each direction, if N <= 1024, and N is rounded up to the nearest multiple of 512 otherwise.  This is to ensure that the FFT is fast.

In order to get the neighborhood type of smoothing of Roberts and Lean (2008) and Ebert (2008), use the boxcar kernel with the argument n giving the neighborhood length.  The resulting kernel is n X n with elements 1/n^2 at each point.  The result is that each grid point of the returned field is an average of the n^2 nearest neighbors.  Alternatively, one might prefer to use a disk kernel, which takes the radius, r, as an argument.  This gives a similar type of kernel, but ensures an average over a uniform distance from the center point.  The disk kernel is also that which is used in the smoothing step of Davis et al (2006a,2006b).  See the help file for \code{kernel2dmeitsjer} for other smoothing options.
}
\value{If setup is FALSE, then a k X m matrix giving the neighborhood smoothed field is returned.  Otherwise, a 2k X 2m possibly complex matrix giving the Fourier transformed kernel weights are returned, which can be used to save an FFT in computation time for subsequent calls to this function by supplying the W argument with this result.}
\references{
Barrett, H. H. and Myers, K. J., 2004: \emph{Foundations of Image Science}.  Wiley Series in Pure and Applied Optics, Editor: B. E. A. Saleh, Hoboken, New Jersey, 1540pp.

Davis CA, BG Brown, and RG Bullock, 2006a. Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.

Davis CA, BG Brown, and RG Bullock, 2006b. Object-based verification of precipitation forecasts, Part II: Application to convective rain systems. \emph{Mon. Wea. Rev.}, \bold{134}, 1785--1795.

Ebert E. E., 2008. Fuzzy verification of high resolution gridded forecasts: A review and proposed framework. \emph{Meteorol. Appl.}, \bold{15}, 51--64. DOI: 10.1002/met.25.  Available at http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf

Hastie, T. J. and Tibshirani, R. J., 1990: \emph{Generalized Additive Models}. Chapman \& Hall/CRC Monographs on Statistics and Applied Probability 43, 335pp.

Ritter, G. X. and Wilson, J. N., 2001: \emph{Handbook of Computer Vision Algorithms in Image Algebra}.  2nd Edition, CRC Press, Boca Raton, Florida, U.S.A., 417pp.

Roberts, N. M. and H. W. Lean, 2008: Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events.  \emph{Mon. Wea. Rev.}, \bold{136}, 78--97.  DOI: 10.1175/2007MWR2123.1.
  
}
\author{Eric Gilleland}
\note{If n is 1, then the field is returned without applying any smoothing.}
\seealso{
\code{\link{fft}}, \code{\link{image.smooth}}, \code{\link{hoods2d}}, \code{\link{Fourier2d}}, \code{kernel2dmeitsjer}
}
\examples{
look <- matrix( 0, 10, 12)
look[4,7] <- 1
kernel2dsmooth( look, kernel.type="boxcar", n=3)
# The above returns the shape of the kernel applied, which
# is a square of length 3 centered on the grid point in look
# that has a value of 1.

# What happens if the 1 is on the edge?  the effect is zero padding:
look <- look*0
look[1,1] <- 1
kernel2dsmooth( look, kernel.type="boxcar", n=3)

# Suppose we want to do the above for several, say l, neighborhood lengths.
# We can save an FFT for l-1 of the convolutions.
look <- look*0
look[4,7] <- 1
lookFFT <- Fourier2d( look, kdim=c(3,3))
dim( lookFFT) # Note the dimension is twice that of look.
kernel2dsmooth( look, kernel.type="boxcar", n=3, X=lookFFT)

# Now, suppose we want to apply the same kernel smooth to different fields.
# We can save an FFT for each subsequent calculation as follows.
wg <- kernel2dsmooth( look, kernel.type="boxcar", n=3, setup=TRUE)
dim( wg) # Note the dimension is twice that of look.
kernel2dsmooth( look, kernel.type="boxcar", n=3, W=wg)
look <- look*0
look[8,5] <- 1
kernel2dsmooth( look, kernel.type="boxcar", n=3, W=wg)
look[5, 10] <- 1
kernel2dsmooth( look, kernel.type="boxcar", n=3, W=wg)

}
\keyword{ math }
\keyword{ manip }
