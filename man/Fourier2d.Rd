\name{Fourier2d}
\alias{Fourier2d}
\title{
Fast Fourier Transform (FFT) of a 2-d Gridded Field
}
\description{
Function to compute the Fast Fourier Transform (FFT) of a gridded field.  The field is first expanded, and zero padded.
}
\usage{
Fourier2d(x, bigdim = NULL, kdim = NULL)
}
\arguments{
  \item{x}{n by m matrix giving the gridded field to transform.}
  \item{bigdim}{ (optional) numeric vector of length 2.  If NULL, it will be found from kdim.  One of this argument or kdim must be given.}
  \item{kdim}{ (optional) numeric vector of length 2 giving the dimension of the kernel to be applied.  If bigdim is NULL, then the optimal dimension for the FFT is calculated from the dimension of x and kdim.}
}
\details{
The gridded field is expanded to bigdim, which is the nearest power of 2 above N=dim( x) + kdim - 1 if N <= 1024, or to the nearest multiple of 512 otherwise, if bigdim is not supplied.  \code{x} is placed in the upper left corner of a matrix of zeros that has dimension bigdim.  Missing values are replaced by zero.  The FFT is conducted on the expanded/zero-padded gridded field.  This is intended to be used more internally for the \code{kernel2dsmooth} function in order to reduce the number of FFT's that need be calculated (e.g., if performed for the same threshold over several neighborhood lengths).  It is currently not used by SpatialVx function anymore, but may still be usefull.
}
\value{
A possibly complex matrix of size bigdim of the Fourier transformed gridded field is returned.
}

\author{
Eric Gilleland
}

\seealso{
\code{\link{fft}}, \code{\link{kernel2dsmooth}}, \code{\link{hoods2d}}, \code{\link{kernel2dmeitsjer}}
}
\examples{

look <- matrix( 0, 10, 12)
look[4,7] <- 1
lookFFT <- Fourier2d( look, kdim=c(3,3))
# Now, call 'kernel2dsmooth' with a neighborhood boxcar kernel that averages the
# nearest grid squares (i.e., neighborhood length of 3).  That is, a square
# of 1/(3^2) = 1/9 ~ 0.111111 with length 3 surrounding the 1 in 'look'.
kernel2dsmooth( look, kernel.type="boxcar", n=3, X=lookFFT)

# Note that the above could have been done with just:
kernel2dsmooth( look, kernel.type="boxcar", n=3)
# But, in the previous one, one less FFT was calculated.
}
\keyword{ manip }
\keyword{ math }
