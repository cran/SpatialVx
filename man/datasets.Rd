\name{obs0601}
\alias{obs0601}
\alias{wrf4ncar0531}
\alias{geom000}
\alias{geom001}
\alias{geom002}
\alias{geom003}
\alias{geom004}
\alias{geom005}
\alias{ICPg240Locs}
\docType{data}
\title{
Spatial Forecast Verification Methods Inter-Comparison Project (ICP) Test Cases and other example verification sets
}
\description{
Test cases used for the ICP.  In particular, those actually analyzed in the special collection of the journal, Weather and Forecasting.  Includes the nine \dQuote{real} cases, five simple geometric cases, and the seven perturbed \dQuote{real} cases.
}
\usage{
data( "obs0601" )
data( "wrf4ncar0531" )

data( "geom000" )
data( "geom001" )
data( "geom002" )
data( "geom003" )
data( "geom004" )
data( "geom005" )

data( "ICPg240Locs" )
}
\format{
  The format is:
 num [1:601, 1:501] 0 0 0 0 0 0 0 0 0 0 ...

The format is:
num [1:301101, 1:2] -110 -110 -110 -110 -110 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "lon" "lat"

}

\details{
One of the nine ICP \dQuote{real} cases is from one version of the Weather Research Forecast (WRF) model denoted \code{wrf4ncar} (see Kain et al. 2008; Ahijevych et al., 2009 for complete details), and the corresponding \dQuote{observed} field is stage II reanalyses denoted here by \dQuote{obs}.  The model is a 24-h forecast so that the valid time is for the next day (e.g., \code{obs0531} corresponds with \code{obs0601}).

These data are a subset from the 2005 Spring Program of the Storm Prediction Center/National Severe Storms Laboratory (SPC/NSSL, cf. Weiss et al., 2005; Kain et al., 2008).  Units for the real cases are in mm/h, and are on the NCEP g240 grid (~4-km resolution) with 601 X 501 grid points.  Both SPC and NSSL should be cited as sources for these cases, as well as Weiss et al. (2005) and possibly also Kain et al. (2008).  The data were made available to the ICP by M. E. Baldwin.

The five geometric cases are simple ellipses (each with two intensities) that are compared against the verification case (\code{geom000}) on the same NCEP g240 grid as the nine real cases.  See Ahijevych et al. (2009) for complete details.  Case \code{geom001} is exactly the same as \code{geom000}, but is displaced 50 grid points to the right (i.e., ~200 km too far east).  Case \code{geom002} is also identical to \code{geom000}, but displaced 200 grid points to the right.  case \code{geom003} is displaced 125 grid points to the right, and is also too big 9i.e., has a spatial extent, or coverage, bias).  Case \code{geom004} is also displaced 125 grid points to the right, but also has a different orientation (note, however, that it is not a true rotation of \code{geom000}).  Case \code{geom005} is displaced 125 grid points to the right, and has a huge spatial extent bias.  This last case is also the only one that actually overlaps with \code{geom000}, and therefore may be regarded by some as the best case.  It is certainly the case that comes out on top by the traditional verification statistics that are calculated on a grid point by grid point basis.  Ahijevych et al. (2009) should be cited if these geometric cases are used for publications, etc.


The longitude and latitude information for each grid (the NCEP g240 grid) is contained in the \code{ICPg240Locs} dataset.

Other data sets for the ICP can be obtained from the ICP web site (\url{https://ral.ucar.edu/projects/icp/}).  MesoVICT data sets are also available there.  All of the ICP test cases used to be available in this package, but had to be removed because of space concerns on CRAN.

}
\source{
\url{https://ral.ucar.edu/projects/icp/}
}
\references{
Ahijevych, D., Gilleland, E., Brown, B. G. and Ebert, E. E. (2009) Application of spatial verification methods to idealized and NWP gridded precipitation forecasts. \emph{Wea. Forecasting}, \bold{24} (6), 1485--1497.

Kain, J. S., Weiss, S. J., Bright, D. R., Baldwin, M. E. Levit, J. J. Carbin, G. W. Schwartz, C. S. Weisman, M. L. Droegemeier, K. K. Weber, and D. B. Thomas, K. W. (2008)  Some Practical Considerations Regarding Horizontal Resolution in the First Generation of Operational Convection-Allowing NWP. \emph{Wea. Forecasting}, \bold{23}, 931--952.

Weiss, S., Kain, J. Levit, J. Baldwin, M. E., Bright, D. Carbin, G. and Hart, J. (2005) NOAA Hazardous Weather Testbed. SPC/NSSL Spring Program 2005 Program Overview and Operations Plan. 61pp.
}
\examples{
\dontrun{
data( "obs0601" )
data( "wrf4ncar0531" )
data( "ICPg240Locs" )
## Plot verification sets with a map.
## Two different methods.

# First way does not preserve projections.
locr <- c( range( ICPg240Locs[,1]), range( ICPg240Locs[,2]))
zl <- range( c( c(obs0601), c( wrf4ncar0531) ) )
par( mfrow=c(2,1), mar=rep(0.1,4))
image( obs0601, axes=FALSE, col=c("grey", tim.colors(256)), zlim=zl)
par( usr=locr)
# if( map.available) map( add=TRUE, database="state") # from library( "maps" )
image( wrf4ncar0531, axes=FALSE, col=c("grey", tim.colors(256)), zlim=zl)
par( usr=locr)
# if( map.available) map( add=TRUE, database="state")
image.plot( obs0601, legend.only=TRUE, horizontal=TRUE, 
		col=c("grey", tim.colors(256)), zlim=zl)

# Second way preserves projections, but values are slighlty interpolated.
zl <- range( c( c(obs0601), c( wrf4ncar0531) ) )
par( mfrow=c(2,2), mar=rep(2.1,4))
image(as.image(c(t(obs0601)), x=ICPg240Locs, nx=601, ny=501, na.rm=TRUE), zlim=zl,
        col=c("grey", tim.colors(64)), axes=FALSE, main="Stage II Reanalysis 4/26/05 0000 UTC")
# map(add=TRUE, lwd=1.5)
# map(add=TRUE, database="state", lty=2)
image(as.image(c(t(wrf4ncar0531)), x=ICPg240Locs, nx=601, ny=501, na.rm=TRUE), zlim=zl,
        col=c("grey", tim.colors(64)), axes=FALSE, main="WRF NCAR valid 4/26/05 0000 UTC")
image.plot(obs0601, col=c("grey", tim.colors(64)), zlim=zl, legend.only=TRUE, horizontal=TRUE)
}
}


\keyword{datasets}
