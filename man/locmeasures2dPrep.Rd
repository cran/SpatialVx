\name{locmeasures2dPrep}
\alias{locmeasures2dPrep}
\title{
Set up a list object for locmeasures2d
}
\description{
Sets up a list object to run with locmeasures2d.
}
\usage{
locmeasures2dPrep(Vx.name, Fcst.name, thresholds = NULL, k = NULL, alpha = 0.1,
    bdconst = NULL, p = 2, loc = NULL, qs = NULL, units = NULL)
}
\arguments{
  \item{Vx.name, Fcst.name}{
character giving the name of a verification and forecast field.  Must have same dimension as each other.
}
  \item{thresholds}{
numeric vector giving the threhsold values over which (inclusive) binary fields are created.  If NULL, quantiles for each field are calculated.  Specifically the 0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9 and 0.95 quantiles.
}
  \item{k}{
numeric vector for use with the partial Hausdorff distance.  For k that are whole numerics or integers >= 1, then the k-th highest value is returned by \code{locmeasures2d}.  If 0 <= k < 1, then the corresponding quantile is returned.
}
  \item{alpha}{numeric giving the alpha parameter for Pratt's Figure of Merit (FOM).  See the help file for \code{locperf} for more details.}
  \item{bdconst}{
numeric giving the cut-off value for Baddeley's delta metric.
}
  \item{p}{
numeric vector giving one or more values for the parameter p in Baddeley's delta metric.  Usually this is just 2.
}
  \item{loc}{
(optional) n by 2 matrix giving the lon/lat coordinates of the grid points (currently not used).
}
  \item{qs}{
character vector giving alternative names for the thresholds.  If \code{thresholds} is NULL and \code{qs} is NULL, then this will give the quantiles used.
}
  \item{units}{
character describing the units for the thresholds.
}
}
\value{
A list object withcomponents:
\item{data.name}{character vector giving the verification and forecast object names.}
\item{thresholds }{n by 2 matrix of thresholds.}
\item{qs}{ same as argument above.}
\item{loc}{ same as argument above.}
\item{units}{same as argument above.}
\item{xdim}{numeric vector of length 2 giving the dimensions of the forecast field (should be the same as the verification field).}
\item{bdconst}{same as argument above.}
\item{p}{same as argument above.}
\item{k}{same as argument above.}
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{locmeasures2d}}
}
\examples{
  data( geom000)
  data( geom001)
  hold <- locmeasures2dPrep( "geom000", "geom001", thresholds=c(0.1,50.1), k=4, units="in/100")
}
\keyword{ list }
