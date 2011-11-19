\name{vxstats}
\alias{vxstats}
\title{
Calculate some common traditional forecast verification statistics.
}
\description{
Calculates some common traditional forecast verification statistics.
}
\usage{
vxstats(Fcst, Obs, which.stats = c("bias", "ts", "ets", "pod", "far", "f", "hk", "mse"), subset = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fcst}{
k X m matrix of forecast values.
}
  \item{Obs}{
k X m matrix of observed values.
}
  \item{which.stats}{
character vector giving the names of the desired statistics.  See Details below.
}
  \item{subset}{numeric vector indicating a subset of the verification set over which to calculate the verification statistics.}
}
\details{
The possible statistics that can be computed, as determined by 'which.stats' are:

"bias" the number of forecast events divided by the number of observed events (sometimes called frequency bias).

"ts" threat score, given by hits/(hits + misses + false alarms)

"ets" equitable threat score, given by (hits - hits.random)/(hits + misses + false alarms - hits.random), where 'hits.random' is the number of observed events times the number of forecast events divided by the total number of forecasts.

"pod" probability of detecting an observed event (aka, hit rate).  It is given by hits/(hits + misses).

"far" false alarm ratio, given by (false alarms)/(hits + false alarms).

"f" false alarm rate (aka probability of false detection) is given by (false alarms)/(correct rejections + false alarms).

"hk" Hansen-Kuipers Score is given by the difference between the hit rate ("pod") and the false alarm rate ("f").

 "mse" mean square error (not a contingency table statistic, but can be used with binary fields).  This is the only statistic that can be calculated here that does not require binary fields for 'Fcst' and 'Obs'.
}
\value{ A list with components determined by which.stats, which may include any or all of the following.
\item{bias }{numeric giving the frequency bias.}
\item{ts }{numeric giving the threat score.}
\item{ets}{numeric giving the equitable threat score, also known as the Gilbert Skill Score.}
\item{pod}{numeric giving the probability of decking an event, also known as the hit rate.}
\item{far}{numeric giving the false alarm ratio.}
\item{f}{numeric giving the false alarm rate.}
\item{hk}{numeric giving the Hanssen and Kuipers statistic.}
\item{mse}{numeric giving the mean square error.}
}
\author{Eric Gilleland}
\note{See the web page: http://www.cawcr.gov.au/projects/verification/ for more details about these statistics, and references.}
\section{Warning }{ It is up to the user to provide the appropriate type of fields for the given statistics to be computed.  For example, they must be binary for all types of which.stats except mse.} 

\seealso{
\code{\link{hoods2d}}
}
\examples{
# Calculate the traditional verification scores for the first geometric case
# of the ICP.
data(geom001)
data(geom000)
rmse <- sqrt(vxstats( geom001, geom000, which.stats="mse")$mse)
rmse
vxstats( geom001 > 0, geom000 > 0, which.stats=c("bias", "ts", "ets", "pod", "far", "f", "hk"))
}

\keyword{ math }
