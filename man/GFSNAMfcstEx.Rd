\name{GFSNAMfcstEx}
\alias{GFSNAMfcstEx}
\alias{GFSNAMobsEx}
\alias{GFSNAMlocEx}
\docType{data}
\title{
Example Verification Set
}
\description{ 
Example verification set of accumulated precipitation (mm) with 361 time points in addition to 2352 spatial locations on a grid.  Taken from a real, but unknown, weather model and observation analysis (one of GFS or NAM).  Accumulation is either 3-h or 24-h.
}
\usage{
data(GFSNAMfcstEx)
data(GFSNAMobsEx)
data(GFSNAMlocEx)
}
\format{
  The format is:
 num [1:2352, 1:361] 0 0 0 0 0 ...

  The format is:
 num [1:2352, 1:361] 0 0 0 0 0 ...

  The format is: A data frame with 2352 observations on the following 2 variables.
  \describe{
    \item{\code{Lat}}{a numeric vector of latitude coordinates for GFS/NAM example verification set.}
    \item{\code{Lon}}{a numeric vector of longitude coordinates for GFS/NAM example verification set.}
  }

}

\details{Example verification set with 2352 spatial locations over the United States, and 361 time points.  For both the forecast (\code{GFSNAMfcstEx}) and verification (\code{GFSNAMobsEx}), these are numeric matrices whose rows represent time, and columns represent space.  The associated lon/lat coordinates are provided by \code{GFSNAMlocEx} (2352 by 2 data frame with named components giving the lon and lat values).

Note that the available spatial locations are a subset of the original 70 X 100 60-km grid where each time point had no missing observations.  This example set is included with the package simply toi demonstrate some functionality that involves both space and time; though this is mostly a spatial-only package.
}

\examples{
data( "GFSNAMfcstEx" )
data( "GFSNAMobsEx" )
data( "GFSNAMlocEx" )

x <- colMeans(GFSNAMfcstEx, na.rm=TRUE)
y <- colMeans(GFSNAMobsEx, na.rm=TRUE)
look <- as.image(x - y, x=GFSNAMlocEx)
image.plot(look)
}
\keyword{datasets}
