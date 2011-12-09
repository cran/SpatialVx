\name{obs0426}
\alias{obs0426}
\alias{obs0513}
\alias{obs0514}
\alias{obs0518}
\alias{obs0519}
\alias{obs0525}
\alias{obs0601}
\alias{obs0603}
\alias{obs0604}
\alias{wrf2caps0425}
\alias{wrf2caps0512}
\alias{wrf2caps0513}
\alias{wrf2caps0517}
\alias{wrf2caps0518}
\alias{wrf2caps0524}
\alias{wrf2caps0531}
\alias{wrf2caps0602}
\alias{wrf2caps0603}
\alias{wrf4ncar0425}
\alias{wrf4ncar0512}
\alias{wrf4ncar0513}
\alias{wrf4ncar0517}
\alias{wrf4ncar0518}
\alias{wrf4ncar0524}
\alias{wrf4ncar0531}
\alias{wrf4ncar0602}
\alias{wrf4ncar0603}
\alias{wrf4ncep0425}
\alias{wrf4ncep0512}
\alias{wrf4ncep0513}
\alias{wrf4ncep0517}
\alias{wrf4ncep0518}
\alias{wrf4ncep0524}
\alias{wrf4ncep0531}
\alias{wrf4ncep0602}
\alias{wrf4ncep0603}
\alias{geom000}
\alias{geom001}
\alias{geom002}
\alias{geom003}
\alias{geom004}
\alias{geom005}
\alias{pert000}
\alias{pert001}
\alias{pert002}
\alias{pert003}
\alias{pert004}
\alias{pert005}
\alias{pert006}
\alias{pert007}
\alias{ICPg240Locs}
\docType{data}
\title{
Spatial Forecast Verification Methods Inter-Comparison Project (ICP) Test Cases and other example verification sets
}
\description{
Test cases used for the ICP.  In particular, those actually analyzed in the special collection of the journal, Weather and Forecasting.  Includes the nine "real" cases, five simple geometric cases, and the seven perturbed "real" cases.
}
\usage{data(obs0426)
	data( ICPg240Locs)
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
The nine "real" cases are forecast model output from three different versions of the Weather Research Forecast (WRF) model denoted "wrf2caps", "wrf4ncar" and "wrf4ncep" (see Kain et al. 2008; Ahijevych et al., 2009 for complete details), and the corresponding "observed" fields are stage II reanalyses denoted here by "obs".  The models are 24-h forecasts so that the valid time is for the next day (e.g., obs0426 corresponds with wrf2caps0425).  The final four digits for the "real" cases give the month and day of the forecast/observation.  These data were from the 2005 Spring Program of the Storm Prediction Center/National Severe Storms Laboratory (SPC/NSSL, cf. Weiss et al., 2005; Kain et al., 2008).  Units for the real cases are in mm/h, and are on the NCEP g240 grid (~4-km resolution) with 601 X 501 grid points.  Both SPC and NSSL should be cited as sources for these cases, as well as Weiss et al. (2005) and possibly also Kain et al. (2008).  The data were made available to the ICP by M. E. Baldwin.

The five geometric cases are simple ellipses (each with two intensities) that are compared against the "observed" case (geom000) on the same NCEP g240 grid as the nine real cases.  See Ahijevych et al. (2009) for complete details.  Case geom001 is exactly the same as geom000, but is displaced 50 grid points to the right (i.e., ~200 km too far east).  Case geom002 is also identical to geom000, but displaced 200 grid points to the right.  case geom003 is displaced 125 grid points to the right, and is also too big 9i.e., has a spatial extent, or coverage, bias).  Case geom004 is also displaced 125 grid points to the right, but also has a different orientation (note, however, that it is not a true rotation of geom000).  Case geom005 is displaced 125 grid points to the right, and has a huge spatial extent bias.  This last case is also the only one that actually overlaps with geom000, and therefore may be regarded by some as the best case.  It is certainly the case that comes out on top by the traditional verification statistics that are calculated on a grid point by grid point basis.  Ahijevych et al. (2009) should be cited if these geometric cases are used for publications, etc.

The seven perturbed cases (dented here by pert00X) are perturbations of pert000, which is adapted from wrf2caps0531 (again, see Ahijevych et al, 2009 for more details).  Case pert001 shifts pert000 three grid points to the right, and five grid points down (i.e., ~12 km to the east and ~20 km to the south).  Case pert002 is a shift six points to the right, and ten points down.  Case pert003 is a shift 12 points to the right and 20 points down.  Case 4 is a shift 24 points to the right, and 40 points down.  Case pert005 is a shift 48 points to the right and 80 points down.  Case pert006 is a shift 12 points to the right, and 20 points down, and the entire field is multiplied by 1.5.  Case pert007 has the same spatial displacements as pert003 and pert006, but also subtracts 1.27 from the entire field.  Ahijevych et al. (2009) should be cited if these perturbed cases are used for publications, etc.

The longitude and latitude information for each grid (the NCEP g240 grid) is contained in the ICPg240Locs dataset.

These constitute all of the official data test cases used for the first round of the ICP, and in the special collection of papers for \emph{Wea. Forecasting}.  Other data sets for the ICP can be obtained from the ICP web site (http://www.ral.ucar.edu/projects/icp).  Future data sets for the next round of the ICP will also be available there, and may potentially also be made available here.

}
\source{
http://www.ral.ucar.edu/projects/icp/datacases.html
}
\references{
Ahijevych, D., E. Gilleland, B.G. Brown, and E.E. Ebert, 2009. Application of spatial verification methods to idealized and NWP gridded precipitation forecasts. \emph{Wea. Forecasting}, \bold{24} (6), 1485--1497.

Kain, J. S., S. J. Weiss, D. R. Bright, M. E. Baldwin, J. J. Levit, G. W. Carbin, C. S. Schwartz, M. L. Weisman, K. K. Droegemeier, D. B. Weber, and K. W. Thomas, 2008.  Some Practical Considerations Regarding Horizontal Resolution in the First Generation of Operational Convection-Allowing NWP. \emph{Wea. Forecasting}, \bold{23}, 931--952.

Weiss, S., J. Kain, J. Levit, M. Baldwin, D. Bright, G. Carbin, and J. Hart, 2005: NOAA Hazardous Weather Testbed. SPC/NSSL Spring Program 2005 Program Overview and Operations Plan. 61pp. Available at: \url{http://www.spc.ncep.noaa.gov/exper/Spring_2005/2005_ops_plan.pdf}
}
\examples{
##
## Below is an example of one way to plot the ICP test cases
## with a map.  It does not give the correct projection of the
## grid.  In future, an example with the correct projection
## may be given.  The example is commented out in order for
## the R check system to run quickly, and the examples are not
## imperative here.
\dontrun{
data(obs0426)
data(wrf2caps0425)
data(wrf4ncar0425)
data(wrf4ncep0425)
data(ICPg240Locs)
require( fields)
map.available <- require( maps)
locr <- c( range( ICPg240Locs[,1]), range( ICPg240Locs[,2]))
zl <- range( c( c(obs0426), c( wrf2caps0425), c( wrf4ncar0425),
		c( wrf4ncep0425)))
par( mfrow=c(2,2), mar=rep(0.1,4))
image( obs0426, axes=FALSE, col=c("grey", tim.colors(256)), zlim=zl)
par( usr=locr)
if( map.available) map( add=TRUE, database="state")
image( wrf2caps0425, axes=FALSE, col=c("grey", tim.colors(256)), zlim=zl)
par( usr=locr)
if( map.available) map( add=TRUE, database="state")
image( wrf4ncar0425, axes=FALSE, col=c("grey", tim.colors(256)), zlim=zl)
par( usr=locr)
if( map.available) map( add=TRUE, database="state")
image( wrf4ncep0425, axes=FALSE, col=c("grey", tim.colors(256)), zlim=zl)
par( usr=locr)
if( map.available) map( add=TRUE, database="state")
image.plot( obs0426, legend.only=TRUE, horizontal=TRUE, 
		col=c("grey", tim.colors(256)), zlim=zl)}

}

\keyword{datasets}
