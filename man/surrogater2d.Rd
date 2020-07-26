\name{surrogater2d}
\alias{surrogater2d}
\alias{aaft2d}
\alias{fft2d}
\alias{mae}
\title{
Create Surrogate Fields
}
\description{
Create surrogate fields that have the same power spectrum and pdf as the original field.
}
\usage{
surrogater2d(Im, frac = 0.95, n = 10, lossfun = "mae", maxiter = 100, zero.down = TRUE,
    verbose = FALSE, ...)

aaft2d(Im, bigdim = NULL)

fft2d(x, bigdim = NULL, ...)

mae(x1, x2, ...)
}
\arguments{
  \item{Im}{
matrix from which surrogates are to be made.
}
  \item{x}{matrix to be Fourier transformed.}
   \item{x1,x2}{numeric or array of same dimensions giving the two fields over which to calculate the mean aboslute error.}
  \item{frac}{
single numeric giving the fraction of original amplitudes to maintain.
}
  \item{n}{
single numeric giving the number of surrogate fields to create (should be a whole number).
}
  \item{lossfun}{
character naming the loss function to use in computing the error between simulated surrogate fields in the iterative process.  Default is the mean absolute error given by the \code{mae} function detailed here.
}
  \item{maxiter}{
Maximum number of iterations allowed per surrogate.
}
  \item{zero.down}{
logical, does \code{Im} contain many zeros, and is otherwise positive?  If so, this sets negative numbers and unusually small numbers to zero.
}
   \item{bigdim}{numeric vector of length two giving the dimensions (larger than dimensions of \code{Im}) to compute the FFT's more efficiently (at least potentially).}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{\dots}{
additional arguments: in the case of \code{fft2d}, they are additional arguments to \code{fft} (i.e., to use inverse=TRUE), in the case of \code{surrogater2d}, they are additional arguments to the loss function given by \code{lossfun}, and in the case of \code{mae} (default), these are not used.
}
}
\details{
The \code{fft2d} function was written to simplify some of the code in \code{surrogater2d} and \code{aaft2d}.  It is simply a call to the \R function \code{fft}, but it first resets the dimensions to ones that should maximize the efficiency.  It will also return the dimensions if they are not passed in.

Surrogates are used in non-linear time series analysis to simulate similar time series for hypothesis testing purposes (e.g., Kantz and Schreiber, 1997).  Venugopal et al. (2005) use surrogates of two-dimensional fields as part of their Forecast Quality Index (FQI); which is the intention here.  Theiler et al. (1992) proposed a method known as the amplitude adjusted Fourier transform (AAFT) algorithm, and Schreiber and Schmitz (1996) proposed a modification to this approach in order to obtain surrogates with both the same power spectrum and pdf as the original series.

The AAFT method first renders the original data, denoted here as s_n, Gaussian via a rank ordering based on randomly generated Gaussian simulated data.  The resulting series, s_n'=g(s_n), is Gaussian and follows the same measured time evolution as s_n.  Next, phase randomized surrogates are made for s_n', call them s_n".  The rescaling g is then inverted by rank ordering s_n" according to the distribution of the original data, s_n.  This algorithm yields surrogates with the same pdf of amplitudes as s_n by construction, but typically not the same power spectra.  The algorithm proposed by Schreiber and Schmitz (1996) begins with the AAFT, and then iterates through a further algorithm as follows.

1. Hold a sorted list of s_n and the squared amplitudes of the Fourier transform of s_n, denote them by S2_k.

2. Take a random shuffle without replacement of the data, denote as s_n(0).

3. Take the Fourier transform of s_n(i).

4. Replace the S2_k(i) with S2_k.

5. Inverse the Fourier transform with the replaced amplitudes.

6. Rank order the series from 5 in order to assume exactly the values taken by s_n.

7. Check the accuracy of 6 using a loss function of some sort, and repeat steps 3 through 6 until a desired level of accuracy is achieved.
}
\value{
In the case of surrogater2d: A three dimesnional array of matrices with same dimension as Im, and third dimension giving the n surrogate fields.

In the case of aaft2d: A matrix of the same dimension as Im.

In the case of fft2d: If bigdim is NULL, a list object is returned with components fft and bigdim giving the FFT of x and the larger dimesnions used.  Otherwise, a matrix of dimension x is returned giving the FFT (or inverse FFT) of x.

In the case of mae: a single numeric giving the mean absolute error between x1 and x2.
}
\references{
Kantz, H. and Schreiber, T. (1997) \emph{Nonlinear time series analysis}.  Cambridge University Press, Cambridge, U.K., 304pp.

Schreiber, T. and Schmitz, A. (1996) Improved surrogate data for nonlinearity tests.  \emph{Physical Review Letters}, \bold{77}(4), 635--638.

Theiler, J., Eubank, S. Longtin, A. Galdrikian, B. and Farmer, J. D. (1992) \emph{Physica} (Amsterdam) \bold{58D}, 77.

Venugopal, V., Basu, S. and Foufoula-Georgiou, E. (2005) A new metric for comparing precipitation patterns with an application to ensemble forecasts.  \emph{J. Geophys. Res.}, \bold{110}, D08111, doi:10.1029/2004JD005395, 11pp.
}
\author{
Eric Gilleland, this code was adapted from matlab code written by Sukanta Basu (2007) available at: http://www.ral.ucar.edu/projects/icp/Software/FeaturesBased/FQI/Perturbed.m
}

\seealso{
\code{\link{fft}}, \code{\link{locmeasures2d}}, \code{\link{UIQI}}, \code{\link{ampstats}}
}
\examples{
data( "ExampleSpatialVxSet" )

x <- ExampleSpatialVxSet$vx

z <- surrogater2d( x, zero.down=FALSE, n=3)

\dontrun{
par( mfrow=c(2,2))
image.plot( look)
image.plot( look2[,,1])
image.plot( look2[,,2])
image.plot( look2[,,3])

  }
}
\keyword{ math }
\keyword{ manip }
