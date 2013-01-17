\name{hump}
\alias{hump}
\docType{data}
\title{
Simulated Forecast and Verification Fields
}
\description{
Simulated forecast and verification fields for optical flow example
}
\usage{data(hump)}
\format{
  The format is:
List of 2
 $ initial: num [1:50, 1:50] 202 210 214 215 212 ...
 $ final  : num [1:50, 1:50] 244 252 257 258 257 ...
}
\details{
Although not identically the same as the data used in Fig. 1 of Marzban and Sandgathe (2010), these are forecast data simulated from the self-same distribution and perturbed in the same manner to get the observation.  The component initial is the forecast and final is the observation.

The forecast is on a 50 X 50 grid simulated from a bivariate Gaussian with standard deviation of 11 and centered on the coordinate (10, 10).  The observed field is the same as the forecast field, but shifted one grid length in each direction and has 60 added to it everywhere.
}
\references{
Marzban, C. and Sandgathe, S. (2010)  Optical flow for verification.  \emph{Wea. Forecasting}, \bold{25}, 1479--1494, doi:10.1175/2010WAF2222351.1.
}
\examples{
data(hump)
str(hump)
\dontrun{
initial <- hump$initial
final <- hump$final
look <- OF(final, initial, W=9, verbose=TRUE)
plot(look) # Compare with Fig. 1 in Marzban and Sandgathe (2010).
hist(look) # 2-d histogram.
plot(look, full=TRUE) # More plots.
}
}
\keyword{datasets}
