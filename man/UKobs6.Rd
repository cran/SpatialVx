\name{UKobs6}
\alias{UKobs6}
\alias{UKfcst6}
\alias{UKloc}
\docType{data}
\title{
Example Precipitation Rate Verification Set (NIMROD)
}
\description{
Example precipitation rate verification set from the very short-range mesoscale Numerical Weather Prediction (NWP) system used operationally at the UK Met Office.
}
\usage{data(UKobs6)
	data(UKloc)}
\format{
  The format is:
 chr "UKobs6"

  The format is:
  num [1:65536, 1:2] -11 -10.9 -10.9 -10.8 -10.7 ...
}
\details{
Precipitation rate (mm/h) verification set from the very short-range NWP system called NIMROD used operationally at the UK Met Office, and described in detail in Casati et al. (2004).  In particular, this is case 6 from Casati et al. (2004), showing a front timing error.  These data are made available for scientific purposes only.  Please cite the source
in any papers or presentations.  The proper reference is the U.K. Met Office.

Refer to Casati et al. (2004) for more information on these data.

The original lon/lat information is not available.  'UKloc' was created to match reasonably well with the figures in Casati et al. (2004), but should not be considered definite.
}
\source{
UK Met Office
}
\references{
Casati, B., Ross, G. and Stephenson, D. B. (2004) A new intensity-scale approach for the verification of spatial precipitation forecasts. \emph{Meteorol. Appl.}, \bold{11}, 141--154, doi:10.1017/S1350482704001239.
}
\examples{
data( "UKobs6" )
data( "UKfcst6" )
data( "UKloc" )

zl <- range(c(c(UKobs6), c(UKfcst6)))
par(mfrow=c(1,2))
image(UKobs6, col=c("grey", tim.colors(64)), zlim=zl, main="analysis", axes=FALSE)
par(usr=apply(UKloc, 2, range))
# map(add=TRUE) # from library( "maps" )
image.plot(UKfcst6, col=c("grey", tim.colors(64)), zlim=zl, main="forecast", axes=FALSE)
par(usr=apply(UKloc, 2, range))
# map(add=TRUE)
}
\keyword{datasets}
