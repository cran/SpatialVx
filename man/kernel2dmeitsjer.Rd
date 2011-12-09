\name{kernel2dmeitsjer}
\alias{kernel2dmeitsjer}
\title{
Create a kernel matrix
}
\description{
Create a kernel matrix to be used in a 2-D convolution smooth (e.g., using kernel2dsmooth).
}
\usage{
kernel2dmeitsjer(type = "gauss", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
character name of the kernel to be created.
}
  \item{\dots}{
Other arguments to the specific kernel type.  See Details section below.
}
}
\details{
The specific types of kernels that can be made are as follows.  In each case, h=||x-x_c|| is the distance from the center of the kernel.  Every kernel that requires single numerics nx and ny be specified returns an nx X ny metrix.  Distances h are found by setting up a grid based on 1 to nx and 1 to ny denoting the points as (xgrid, ygrid), finding the center of the grid as (x.center, y.center)=(nx/2,ny/2), and then h = sqrt( (xgrid - x.center)^2 + (ygrid - y.center)^2).  For kernels that better reflect distance (e.g., using great-circle distance, anisotropic distances, etc.), the matrix h can be passed instead of nx and ny, but only for those kernels that take h as an argument.  In each case with sigma as an argument, sigma is the smoothing parameter.  There are many kernel functions allowed here, and not all of them make sense for every purpose.

"average" gives a kernel that will give an average of the nearest neighbors in each direction (can take an average grid points further in the x- direction than the y-direction, and vice versa).  Requires that 'nx' and 'ny' be specified, and the resulting kernel is defined by an 'nx X ny' matrix with each element equal to 1/(nx*ny).  If nx = ny, then the result is the same as the boxcar kernel below.

"boxcar" the boxcar kernel is an n X n matrix of 1/n^2.  This results in a neighborhood smoothing when used with 'kernel2dsmooth' giving the type of smoothed fields utilized, e.g., in Roberts and Lean (2008) and Ebert (2008).  Requires that n be specified.  Note that usually the boxcar is a square matrix of ones, which gives the sum of the nearest n^2 grid points.  This gives the average.

"cauchy" The Cauchy kernel is given by K(sigma)=1/(1+h^2/sigma).  Requires the arguments nx, ny and sigma.  See Souza (2010) for more details.

"disk" gives a circular averaging (pill-box) kernel (aka, disk kernel).  Similar to "average" or "boxcar", but this kernel accounts for a specific distance in all directions from the center (i.e., an average of grid squares within a circular radius of the central point).  This results in the convolution radius smoothing applied in Davis et al. (2006a,2006b).  Requires that 'r' (the desired radius) be supplied, and a square matrix of appropriate dimension is returned.

"epanechnikov" The Epanechnikov kernel is defined by max(0, 3/4*(1-h/(sigma^2))). See, e.g., Hastie and Tibshirani (1990).  Requires arguments 'nx', 'ny', and 'sigma'.

"exponential" The exponential kernel is given by K(sigma) = a*exp(-h/(2*sigma)).  Requires the arguments nx, ny and sigma, and optionally takes the argument 'a' (default is a=1).  An 'nx X ny' matrix is returned.  See Souza (2010) for more details.

"gauss" The Gaussian kernel defined by K(sigma) = 1/(2*pi*sigma^2)*exp(-h/(2*sigma)).  Requires the arguments 'nx', 'ny' and 'sigma' be specified.  The convolution with this kernel results in a Gaussian smoothed field as used in the practically perfect hindcast method of Brooks et al. (1998) (see also Ebert 2008) and studied by Sobash et al (2011) for spatial forecast verification purposes.  Returns an 'nx X ny' matrix.

"laplacian" Laplacian Operator kernel, which gives the sum of second partial derivatives for each direction.  It is often used for edge detection because it identifies areas of rapid intensity change.  Typically, it is first applied to a field that has been smoothed first by a Gaussian kernel smoother (or an approximation thereof; cf. type "LoG" below).  This method optionally the parameter 'alpha', which controls the shape of the Laplacian kernel, which must be between 0 and 1 (inclusive), or else it will be set to 0 (if < 0) or 1 (if > 1).  Returns a 3 X 3 kernel matrix.

"LoG" Laplacian of Gaussian kernel.  This combines the Laplacian Operator kernel with that of a Gaussian kernel.  The form is given by K(sigma) = -1/(pi*sigma^4)*exp(-h/(2*sigma^2))*(1-h/(2*sigma^2)).  Requires the arguments 'nx', 'ny' and 'sigma' be specified.  Returns an 'nx X ny' matrix.

"minvar" A minimum variance kernel, which is given by 3/8*(3 - 5*h/sigma^2) if h <= 1, and zero otherwise (see, e.g., Hastie and Tibshirani, 1990).  Requires the arguments 'nx', 'ny', and 'sigma' be specified.  Returns an 'nx X ny' matrix.

"multiquad" The multiquadratic kernel is similar to the rational quadratic kernel, and is given by K(a) = sqrt(h + a^2).  The inverse is given by 1/K(a).  Requires the arguments nx, ny and a be specified.  Optionally takes a logical named inverse determining whether to return the inverse multiquadratic kernel or not.

"prewitt" Prewitt filter kernel, which emphasizes horizontal (vertical) edges through approximation of a vertical (horizontal) gradient.  Optionally takes a logical argument named 'transpose', which if FALSE (default) emphasis is on horizontal, and if TRUE emphasis is on vertical.  Returns a 3 X 3 matrix whose first row is all ones, second row is all zeros, and third row is all negative ones for the transpose=FALSE case, and the transpose of this matrix in the transpose=TRUE case.

"power" The power kernel is defined by K(p) = -h^p.  The log power kernel is similarly defined as K(p) = -log(h^p+1).  Requires specification of the arguments nx, ny and p.  Alternatively takes the logical do.log to determine whether the log power kernel should be returned (TRUE) or not (FALSE).  Default if not passed is to do the power kernel.  Returns an nx X ny matrix.  See Souza (2010) for more details.

"radial" The radial kernel is returns a*|h|^(2*m-d)*log(|h|) if d is even and a*|h|^(2*m-d) otherwise.  Requires arguments a, m, d nx and ny.  Replaces any missing values with zero.

"ratquad" The rational quadratic kernel is used as an alternative to the Gaussian, and is given by K(a) = 1 - h/(h+a).  Requires the arguments nx, ny and a, and returns an 'nx X ny' matrix.  See Souza (2010) for more details.

"sobel" Same as prewitt, except that the elements 1,2 and 3,2 are replaced by two and neative two, resp.

"student" The generalized Student's t kernel is defined by K(p)=1/(1+h^p).  Requires the arguments nx, ny and p be specified.  Returns an nx X ny matrix.  See Souza (2010) for more details.

"unsharp" Unsharp contrast enhancement filter.  This is simply given by a 3 X 3 matrix of al zeros, except for a one in the center subtracted by a laplacian operator kernel matrix.  Requires the same arguments as for "laplacian".  Returns a 3 X 3 matrix.

"wave" The wave kernel is defined by K(phi) = phi/h * sin( h/phi).  Requires arguments nx, ny and phi be specified.  Returns an nx X ny matrix.
}
\value{
matrix of dimension determined by the specific type of kernel, and possibly user passed arguments giving the kernel to be used by \code{kernel2dsmooth}.
}
\references{
Brooks, H. E., K. M., and J. A. Hart, 1998: Objective limits on forecasting skill of rare events. \emph{19th Conf. Severe Local Storms}. Amer. Met. Soc., 552-555.

Davis CA, BG Brown, and RG Bullock, 2006a. Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.

Davis CA, BG Brown, and RG Bullock, 2006b. Object-based verification of precipitation forecasts, Part II: Application to convective rain systems. \emph{Mon. Wea. Rev.}, \bold{134}, 1785--1795.

Ebert E. E., 2008. Fuzzy verification of high resolution gridded forecasts: A review and proposed framework. \emph{Meteorol. Appl.}, \bold{15}, 51-64. DOI: 10.1002/met.25.  Available at \url{http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf}

Hastie, T. J. and Tibshirani, R. J., 1990: \emph{Generalized Additive Models}. Chapman \& Hall/CRC Monographs on Statistics and Applied Probability 43, 335pp.

Roberts, N. M. and H. W. Lean, 2008: Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events.  \emph{Mon. Wea. Rev.}, \bold{136}, 78-97.  DOI: 10.1175/2007MWR2123.1.

Sobash, R. A., J. S. Kain, D. R. Bright, A. R. Dean, M. C. Coniglio, and S. J. Weiss, 2011: Probabilistic forecast guidance for severe thunderstorms based on the identification of extreme phenomena in convection-allowing model forecasts. \emph{Wea. Forecasting}, \bold{26}, 714-728.

Souza, C\'{e}sar R., 2010.  \emph{Kernel Functions for Machine Learning Applications}. 17 Mar 2010. Web. \url{http://crsouza.blogspot.com/2010/03/kernel-functions-for-machine-learning.html}.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{fft}}, \code{\link{kernel2dsmooth}}, \code{\link{hoods2d}}
}
\examples{

x <- matrix( 0, 10, 12)
x[4,5] <- 1
kmat <- kernel2dmeitsjer( "average", nx=7, ny=5)
kernel2dsmooth( x, K=kmat)

##
## Can als call 'kernel2dsmooth' directly.
##
kernel2dsmooth( x, kernel.type="boxcar", n=5)
kernel2dsmooth( x, kernel.type="cauchy", sigma=20, nx=10, ny=12)
kernel2dsmooth( x, kernel.type="disk", r=3)
kernel2dsmooth( x, kernel.type="epanechnikov", nx=10, ny=12, sigma=4)
kernel2dsmooth( x, kernel.type="exponential", a=0.1, sigma=4, nx=10, ny=12)
kernel2dsmooth( x, kernel.type="gauss", nx=10, ny=12, sigma=4)
kernel2dsmooth( x, kernel.type="laplacian", alpha=0)
kernel2dsmooth( x, kernel.type="LoG", nx=10, ny=12, sigma=1) 
kernel2dsmooth( x, kernel.type="minvar", nx=10, ny=12, sigma=4)
kernel2dsmooth( x, kernel.type="multiquad", a=0.1, nx=10, ny=12)
kernel2dsmooth( x, kernel.type="power", p=0.5, nx=10, ny=12)
kernel2dsmooth( x, kernel.type="prewitt")
kernel2dsmooth( x, kernel.type="prewitt", transpose=TRUE)
kernel2dsmooth( x, kernel.type="radial", a=1, m=2, d=1, nx=10, ny=12)
kernel2dsmooth( x, kernel.type="ratquad", a=0.1, nx=10, ny=12)
kernel2dsmooth( x, kernel.type="sobel")
kernel2dsmooth( x, kernel.type="sobel", transpose=TRUE)
kernel2dsmooth( x, kernel.type="student", p=1.5, nx=10, ny=12)
kernel2dsmooth( x, kernel.type="unsharp", alpha=0)
kernel2dsmooth( x, kernel.type="wave", phi=45, nx=10, ny=12)

\dontrun{
data(lennon)
kmat <- kernel2dmeitsjer( "average", nx=7, ny=5)
lennon.smAvg <- kernel2dsmooth( lennon, K=kmat)
## Can also just make a call to kernel2dsmooth, which
## will call this function.
lennon.smBox <- kernel2dsmooth( lennon, kernel.type="boxcar", n=7)
lennon.smDsk <- kernel2dsmooth( lennon, kernel.type="disk", r=5)
par( mfrow=c(2,2), mar=rep(0.1,4))
image.plot( lennon, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smAvg, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smBox, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smDsk, col=tim.colors(256), axes=FALSE)

lennon.smEpa <- kernel2dsmooth( lennon, kernel.type="epanechnikov", nx=10, ny=10, sigma=20)
lennon.smGau <- kernel2dsmooth( lennon, kernel.type="gauss", nx=10, ny=10, sigma=20)
lennon.smMvr <- kernel2dsmooth( lennon, kernel.type="minvar", nx=10, ny=10, sigma=20)
par( mfrow=c(2,2), mar=rep(0.1,4))
image.plot( lennon, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smEpa, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smGau, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smMvr, col=tim.colors(256), axes=FALSE)

lennon.smLa0 <- kernel2dsmooth( lennon, kernel.type="laplacian", alpha=0)
lennon.smLap <- kernel2dsmooth( lennon, kernel.type="laplacian", alpha=0.999)
lennon.smLoG <- kernel2dsmooth( lennon, kernel.type="LoG", nx=10, ny=10, sigma=20)
par( mfrow=c(2,2), mar=rep(0.1,4))
image.plot( lennon, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smLa0, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smLap, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smLoG, col=tim.colors(256), axes=FALSE)

lennon.smPrH <- kernel2dsmooth( lennon, kernel.type="prewitt")
lennon.smPrV <- kernel2dsmooth( lennon, kernel.type="prewitt", transpose=TRUE)
lennon.smSoH <- kernel2dsmooth( lennon, kernel.type="sobel")
lennon.smSoV <- kernel2dsmooth( lennon, kernel.type="sobel", transpose=TRUE)
par( mfrow=c(2,2), mar=rep(0.1,4))
image.plot( lennon.smPrH, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smPrV, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smSoH, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smSoV, col=tim.colors(256), axes=FALSE)

lennon.smUsh <- kernel2dsmooth( lennon, kernel.type="unsharp", alpha=0.999)
par( mfrow=c(2,1), mar=rep(0.1,4))
image.plot( lennon, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smUsh, col=tim.colors(256), axes=FALSE)

lennon.smRad1 <- kernel2dsmooth( lennon, kernel.type="radial", a=2, m=2, d=1, nx=10, ny=10)
lennon.smRad2 <- kernel2dsmooth( lennon, kernel.type="radial", a=2, m=2, d=2, nx=10, ny=10)
par( mfrow=c(2,1), mar=rep(0.1,4))
image.plot( lennon.smRad1, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smRad2, col=tim.colors(256), axes=FALSE)

lennon.smRQd <- kernel2dsmooth( lennon, kernel.type="ratquad", a=0.5, nx=10, ny=10)
lennon.smExp <- kernel2dsmooth( lennon, kernel.type="exponential", a=0.5, sigma=20, nx=10, ny=10)
lennon.smMQd <- kernel2dsmooth( lennon, kernel.type="multiquad", a=0.5, nx=10, ny=10)
par( mfrow=c(2,2), mar=rep(0.1,4))
image.plot( lennon.smGau, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smRQd, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smExp, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smMQd, col=tim.colors(256), axes=FALSE)

lennon.smIMQ <- kernel2dsmooth( lennon, kernel.type="multiquad", a=0.5, nx=10, ny=10, inverse=TRUE)
par( mfrow=c(2,1), mar=rep(0.1,4))
image.plot( lennon.smMQd, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smIMQ, col=tim.colors(256), axes=FALSE)

lennon.smWav <- kernel2dsmooth( lennon, kernel.type="wave", phi=45, nx=10, ny=10)
par( mfrow=c(1,1), mar=rep(0.1,4))
image.plot( lennon.smWav, col=tim.colors(256), axes=FALSE)

lennon.smPow <- kernel2dsmooth( lennon, kernel.type="power", p=0.5, nx=10, ny=10)
lennon.smLpw <- kernel2dsmooth( lennon, kernel.type="power", p=0.5, nx=10, ny=10, do.log=TRUE)
par( mfrow=c(2,1), mar=rep(0.1,4))
image.plot( lennon.smPow, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smLpw, col=tim.colors(256), axes=FALSE)

lennon.smCau <- kernel2dsmooth( lennon, kernel.type="cauchy", sigma=20, nx=10, ny=10)
lennon.smStd <- kernel2dsmooth( lennon, kernel.type="student", p=1.5, nx=10, ny=10)
par( mfrow=c(2,1), mar=rep(0.1,4))
image.plot( lennon.smCau, col=tim.colors(256), axes=FALSE)
image.plot( lennon.smStd, col=tim.colors(256), axes=FALSE)
}
}
\keyword{ math }
\keyword{ manip }
