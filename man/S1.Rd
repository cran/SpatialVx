\name{S1}
\alias{S1}
\alias{ACC}
\title{
S1 Score, Anomaly Correlation
}
\description{
Calculate the S1 score and anomaly correlation for a verification set.
}
\usage{
S1(X, Y, gradFUN = "KernelGradFUN", ...)
ACC(X, Y, Xclim = NULL, Yclim = NULL)
}
\arguments{
  \item{X,Y}{
  m by n matrices giving the verification and forecast fields, resp.
}
  \item{Xclim,Yclim}{m by n matrices giving the climatologies for \code{X} and \code{Y}, resp.  If NULL, the result is simply a usual correlation.}
  \item{gradFUN}{
character identifying a function used to calculate the gradient fields for \code{X} and \code{Y}.  The default \code{KernelGradFUN} is to use a Laplacian of Gaussian kernel.
}
  \item{\dots}{
optional arguments to the \code{gradFUN} function.  In the case of the default, the kernel can be changed (e.g., if only \dQuote{laplacian} is desired), and optional arguments to the \code{kernel2dmeitsjer} function (in this case, \code{nx}, \code{ny} and \code{sigma}).
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
data(UKobs6)
data(UKfcst6)
S1(UKobs6,UKfcst6)
}
\keyword{ math }
