\name{GFSNAMfcstEx}
\alias{GFSNAMfcstEx}
\alias{GFSNAMobsEx}
\alias{GFSNAMlocEx}
\docType{data}
\title{
GFS NAM example verification set
}
\description{ 
GFS NAM example verification set with 361 time points in addition to 2352 spatial locations on a grid.
}
\usage{data(GFSNAMfcstEx)
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

\details{GFS/NAM example verification set with 2352 spatial locations over the United States, and 361 time points.  For both the forecast (GFSNAMfcstEx) and verification (GFSNAMobsEx), these are numeric matrices whose rows represent time, and columns represent space.  The associated lon/lat coordinates are provided by GFSNAMlocEx (2352 X 2 data frame with named components giving the lon and lat values).}

\examples{
data(GFSNAMfcstEx)
data(GFSNAMobsEx)
data(GFSNAMlocEx)

x <- colMeans(GFSNAMfcstEx,na.rm=TRUE)
y <- colMeans(GFSNAMobsEx,na.rm=TRUE)
look <- as.image(x-y, x=GFSNAMlocEx)
image.plot(look)
}
\keyword{datasets}
