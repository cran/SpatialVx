\name{Fourier2d}
\alias{Fourier2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Fast Fourier Transform (FFT) of a gridded field.
}
\description{
Function to compute the Fast Fourier Transform (FFT) of a gridded field.  the field is first expanded, and zero padded.
}
\usage{
Fourier2d(x, xdim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{n X m matrix giving the gridded field to transform.}
  \item{xdim}{ (optional) numeric vector of length 2.  If NULL, it will be found from dim( x).}
}
\details{
The gridded field is expanded to have twice its original dimension by placing \code{x} into the upper left corner of a matrix of zeros that has dimension 2n X 2m.  Missing values are replaced by zero.  The FFT is conducted on the expanded/zero-padded gridded field.  This is intended to be used more internally for the \code{kernel2d} function in order to reduce the number of FFT's that need be calculated (e.g., if performed for the same threshold over several neighborhood lengths).
}
\value{
A 2n X 2m possibly complex matrix of the Fourier transformed gridded field is returned.
}

\author{
Eric Gilleland
}

\seealso{
\code{\link{fft}}, \code{\link{kernel2d}}, \code{\link{hoods2d}}
}
\examples{

look <- matrix( 0, 10, 12)
look[4,7] <- 1
lookFFT <- Fourier2d( look)
# Now, call 'kernel2d' with a neighborhood boxcar kernel that averages the
# nearest grid squares (i.e., neighborhood length of 3).  That is, a square
# of 1/(3^2) = 1/9 ~ 0.111111 with length 3 surrounding the 1 in 'look'.
kernel2d( look, 3, X=lookFFT)

# Note that the above could have been done with just.
kernel2d( look, 3)
# But, in the previous one, one less FFT was calculated.
}
\keyword{ manip }
\keyword{ math }
