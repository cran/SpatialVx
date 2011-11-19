\name{SpatialVx-package}
\alias{SpatialVx-package}
\alias{SpatialVx}
\docType{package}
\title{
Spatial Forecast Verification
}
\description{
\pkg{SpatialVx} contains functions to perform several spatial forecast verification methods.  the initial release has only functionality to do most of the traditional and neighborhood smoothing methods.  Future releases will be considerably more complete.}

\details{

Primary functions include: \code{hoods2dPrep}, \code{hoods2d}, \code{pphindcast2d}, \code{kernel2d}, and \code{plot.hoods2d}.

All of the initial Spatial Forecast Verification Inter-Comparison Project (ICP, http://www.ral.ucar.edu/projects/icp) data sets used in the special collection of the Weather and Forecasting journal are included.  See the help file for 'obs0426', which will give information on all of the datasets included, as well as an example for plotting the data.

Ebert (2008) provides a nice review of these methods.  Roberts and Lean (2008) describes one of the methods, as well as the primary boxcar kernel smoothing method used throughout this package.  Gilleland et al. (2009, 2010) provides an overview of most of the various recently proposed methods, and Ahijevych et al. (2009) describes the data sets included in this package.  Some of these have been applied to the ICP test cases in Ebert (2009).
}

\author{
Eric Gilleland
}
\references{
Ahijevych, D., E. Gilleland, B.G. Brown, and E.E. Ebert, 2009. Application of spatial verification methods to idealized and NWP gridded precipitation forecasts. \emph{Wea. Forecasting}, \bold{24} (6), 1485--1497.

Ebert EE, 2008. Fuzzy verification of high resolution gridded forecasts: A review and proposed framework.  \emph{Meteorol. Appl.}, \bold{15}, 51--64. DOI: 10.1002/met.25 Available at http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf

Ebert, E. E., 2009: Neighborhood verification: A strategy for rewarding close forecasts.  \emph{Wea. Forecasting}, \bold{24}, 1498-1510, DOI: 10.1175/2009WAF2222251.1.

Gilleland, E., D. Ahijevych, B.G. Brown, B. Casati, and E.E. Ebert, 2009. Intercomparison of Spatial Forecast Verification Methods. \emph{Wea. Forecasting}, \bold{24}, 1416--1430, DOI: 10.1175/2009WAF2222269.1.

Gilleland, E., D.A. Ahijevych, B.G. Brown and E.E. Ebert, 2010: Verifying Forecasts Spatially. \emph{Bull. Amer. Meteor. Soc.}, October, 1365--1373.

Roberts, N. M. and H. W. Lean, 2008: Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events.  \emph{Mon. Wea. Rev.}, \bold{136}, 78--97. DOI: 10.1175/2007MWR2123.1.
}
\examples{
## The example below may take a few minutes to run,
## and therefore is commented out in order for R's
## package checks to run in a shorter amount of time.
## Compare the results to Ebert (2009) Fig. 5.
##
## Note that 'levels' refer to number of grid points,
## which for this example are ~4 km, so that the
## plot labels differ by a factor of 4 from Ebert (2009).
## Future versions of this package will allow for different
## labelling.
\dontrun{
data( pert004)
data(pert000)
hold <- hoods2dPrep( "pert004", "pert000", thresholds=c(1, 2, 5, 10, 20, 50), levels=c(1, 3, 5, 9, 17, 33, 65, 129, 257), units="mm/h")
look <- hoods2d( hold, which.methods=c("fss", "multi.event"), verbose=TRUE)
plot( look)
look2 <- pphindcast2d( hold, verbose=TRUE)
look2
}
}
\keyword{ package}
