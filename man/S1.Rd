\name{S1}
\alias{S1}
\alias{S1.default} 
\alias{S1.SpatialVx}
\alias{ACC}
\alias{ACC.default}
\alias{ACC.SpatialVx}
\title{
S1 Score, Anomaly Correlation
}
\description{
Calculate the S1 score and anomaly correlation for a verification set.
}
\usage{
S1(x, ...)

\method{S1}{default}(x, ..., xhat, gradFUN = "KernelGradFUN")

\method{S1}{SpatialVx}(x, ..., xhat, gradFUN = "KernelGradFUN",
    time.point = 1, obs = 1, model = 1)

ACC(x, ...)

\method{ACC}{default}(x, ..., xhat, xclim = NULL, xhatclim = NULL)

\method{ACC}{SpatialVx}(x, ..., xclim = NULL, xhatclim = NULL,
    time.point = 1, obs = 1, model = 1)
}
\arguments{
  \item{x,xhat}{
  m by n matrices giving the verification and forecast fields, resp.

For \code{S1.SpatialVx} and \code{ACC.SpatialVx}, \code{x} is an object of class \dQuote{SpatialVx}.
}
  \item{xclim, xhatclim}{m by n matrices giving the climatologies for \code{X} and \code{Y}, resp.  If NULL, the result is simply a usual correlation.}
  \item{gradFUN}{
character identifying a function used to calculate the gradient fields for \code{X} and \code{Y}.  The default \code{KernelGradFUN} is to use a Laplacian of Gaussian kernel.
}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
    \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
  \item{\dots}{
optional arguments to the \code{gradFUN} function.  In the case of the default, the kernel can be changed (e.g., if only \dQuote{laplacian} is desired), and optional arguments to the \code{kernel2dmeitsjer} function (in this case, \code{nx}, \code{ny} and \code{sigma}).  Not used by \code{ACC}.
}
}
\details{
The S1 score is given by

	S1 = 100*sum(abs(DY_i - DX_i))/sum(max(abs(DY_i),abs(DX_i))),

where DY_i (DX_i)is the gradient at grid point i for the forecast (verification).  See Brown et al. (2012) and Thompson and Carter (1972) for more on this score.

The ACC is just the correlation between X - Xclim and Y - Yclim.
}
\value{
single numeric
}
\references{
Brown, B.G., Gilleland, E. and Ebert, E.E. (2012) Chapter 6: Forecasts of spatial fields. pp. 95--117, In \emph{Forecast Verification: A Practitioner's Guide in Atmospheric Science}, 2nd edition. Edts. Jolliffee, I. T. and Stephenson, D. B., Chichester, West Sussex, U.K.: Wiley, 274 pp.

Thompson, J. C. and Carter, G. M. (1972) On some characteristics of the S1 score.  \emph{J. Appl. Meteorol.}, \bold{11}, 1384--1385.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{kernel2dmeitsjer}}
}
\examples{
data( "UKobs6" )
data( "UKfcst6" )

S1( UKobs6, xhat = UKfcst6 )
ACC( UKobs6, xhat = UKfcst6 )

\dontrun{
data( "obs0601" )
data( "wrf4ncar0531" )
data( "ICPg240Locs" )

hold <- make.SpatialVx( obs0601, wrf4ncar0531, loc = ICPg240Locs,
    projection = TRUE, map = TRUE, loc.byrow = TRUE,
    field.type = "Precipitation", units = "mm/h",
    data.name = "ICP NSSL/SPC Spring 2005 Cases",
    obs.name = "obs0601", model.name = "wrf4ncar0531" )

plot( hold )

S1( hold )
ACC( hold )
}
}
\keyword{ math }
