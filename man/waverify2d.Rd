\name{waverify2d}
\alias{waverify2d}
\alias{waverify2d.default}
\alias{waverify2d.SpatialVx}
\alias{mowaverify2d}
\alias{mowaverify2d.default}
\alias{mowaverify2d.SpatialVx}
\alias{plot.waverify2d}
\alias{print.waverify2d}
\title{
High-Resolution Gridded Forecast Verification Using Discrete Wavelet Decomposition
}
\description{
High-resolution gridded forecast verification using discrete wavelet decomposition.
}
\usage{
waverify2d(X, ..., Clim = NULL, wavelet.type = "haar", J = NULL)

\method{waverify2d}{default}(X, ..., Y, Clim = NULL, wavelet.type = "haar", J =
                 NULL, useLL = FALSE, compute.shannon = FALSE,
                 which.space = "field", verbose = FALSE)

\method{waverify2d}{SpatialVx}(X, ..., Clim = NULL, wavelet.type = "haar", J = NULL, 
    useLL = FALSE, compute.shannon = FALSE, which.space = "field", 
    time.point = 1, obs = 1, model = 1, verbose = FALSE)

mowaverify2d(X, ..., Clim = NULL, wavelet.type = "haar", J = 4)

\method{mowaverify2d}{default}(X, ..., Clim = NULL, Y, wavelet.type = "haar", J = 4,
    useLL = FALSE, compute.shannon = FALSE, which.space = "field", verbose = FALSE)

\method{mowaverify2d}{SpatialVx}(X, ..., Clim = NULL, wavelet.type = "haar", J = 4, 
    useLL = FALSE, compute.shannon = FALSE, which.space = "field", 
    time.point = 1, obs = 1, model = 1, verbose = FALSE)

\method{plot}{waverify2d}(x, ..., main1 = "X", main2 = "Y", main3 = "Climate", 
    which.plots = c("all", "dwt2d", "details", "energy", "mse", 
        "rmse", "acc"), separate = FALSE, col, horizontal = TRUE)

\method{print}{waverify2d}(x, ...)
}
\arguments{
  \item{X,Y,Clim}{
   m by n dyadic matrices (i.e., m = 2^M and n = 2^N, for M, N some integers) giving the verification and forecast fields (and optionally a climatology field), resp.  Alternatively, \code{X} may be a \dQuote{SpatialVx} object, in which case, \code{Y} is not given and in either case \code{Clim} must be provided if it is to be used.
}
   \item{x}{
list object of class \dQuote{waverify2d} as returned by \code{waverify2d}.
}
  \item{wavelet.type}{
character naming the type of wavelet to be used.  This is given as the \code{wf} argument to the \code{dwt.2d} function of package \pkg{waveslim}.
}
  \item{J}{
(optional) numeric integer giving the pre-determined number of levels to use.  If NULL, J is set to be log2(m) = M in \code{waverify2d} only.
}
  \item{useLL}{
logical, should the LL submatrix (i.e., the father wavelet or grand mean) be used to find the inverse DWT's for calculating the detail fields?
}
  \item{compute.shannon}{
logical, should the Shannon entropy be calculated for the wavelet decomposition?
}
  \item{which.space}{character (one of \dQuote{field} or \dQuote{wavelet}) naming from which space the detail fields should be used.  If \dQuote{field}, then it is in the original field (or image) space (i.e., the detail reconstruction), and if \dQuote{wavelet}, it will be done in the wavelet space (i.e., the detail wavelet coefficients).}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
  \item{main1,main2,main3}{optional characters naming each field to be used for the detail field plots and legend labelling on the energy plot.}
   \item{which.plots}{character vector describing which plots to make.  The default is to make all of them.  \dQuote{dwt2d} option uses \code{plot.dwt2d} from \pkg{waveslim}.  \dQuote{details} option makes image plots of the detail fields on which the various statistics are calculated.  The rest of the options give line plots showing the statistics.}
  \item{separate}{logical, should the plots be on their own devices (TRUE) or should some of them be put onto a single multi-panel device (FALSE, default)?}
    \item{col}{optional argument specifying the \code{col} argument in calls to functions like \code{image}.  Default is a concatenation of gray with \code{time.colors(64)}.}
    \item{horizontal}{logical, should the legend on image plots be horizontal (TRUE, placed at the bottom of the plot) or vertical (FALSE, placed at the right side of the plot)?}
  \item{verbose}{logical, should progress information be printed to the screen, including total run time?
}
\item{\dots}{
optional additonal plot or image.plot parameters.  If detail and energy, mse, rmse or acc plots are desired, must be applicable to both types of plots.  Not used by \code{print} method function.
}
}
\details{
This is a function to use discrete wavelet decomposition to analyze verification sets along the lines of Briggs and Levine (1997), as well as Casati et al. (2004) and Casati (2009).  In the originally proposed formulation of Briggs and Levine (1997), continuous verification statistics (namely, the anomaly correlation coefficient (ACC) and root mean square error (RMSE)) are calculated for detail fields obtained from wavelet decompositions of each of a forecast and verification field (and for ACC a climatology field as well).  Casati et al. (2004) introduced an intensity scale approach that applies 2-d DWT to binary (obtained from thresholding) difference fields (Forecast - Verification), and applying a skill score at each level based on the mean square error (MSE).  Casati (2009) extended this idea to look at the energy at each level as well.

This function makes use of the \code{dwt.2d} and \code{idwt.2d} functions from package \pkg{waveslim}, and \code{plot.waverify2d} uses the \code{plot.dwt.2d} function if \code{dwt2d} is selected through the \code{which.plots} argument.  See the help file for these functions, the references therein and the references herein for more on these approaches.

Generally, it is not necessary to use the father wavelet for the detail fields, but for some purposes, it may be desired.

\code{mowaverify2d} is very similar to \code{waverify2d}, but it allows fields to be non-dyadic (and may subsequently be slower).  It uses the \code{modwt.2d} and \code{imodwt.2d} functions from the package \pkg{waveslim}.  In particular, it performs a maximal overlap discrete wavelet transform on a matrix of arbitrary dimension.  See the help file and references therein for \code{modwt.2d} for more information, as well as Percival and Guttorp (1994) and Lindsay et al. (1996).

In Briggs and Levine (1997), they state that the calculations can be done in either the data (called field here) space or the wavelet space, and they do their examples in the field space.  If the wavelets are orthogonal, then the detail coefficeints (wavelet space), can be analyzed with the assumption that they are independent; whereas in the data space, they typically cannot be assumed to be independent.  Therefore, most statistical tests should be performed in the wavelet space to avoid issues arising from spatial dependence.
}
\value{
A list object of class \dQuote{waverify2d} with components:
  \item{J }{single numeric giving the number of levels.}
  \item{X.wave, Y.wave, Clim.wave }{objects of class \dQuote{dwt.2d} describing the wavelet decompositions for the verification and forecast fields (and climatology, if applicable), resp. (see the help file for dwt.2d from package waveslim for more about these objects).}
   \item{Shannon.entropy}{numeric matrix giving the Shannon entropy for each field.}
   \item{energy}{numeric matrix giving the energy at each level and field.}
   \item{mse,rmse}{numeric vectors of length J giving the MSE/RMSE for each level between the verification and forecast fields.}
   \item{acc}{If a climatology field is supplied, this is a numeric vector giving the ACC for each level.}
}
\references{
Briggs, W. M. and Levine, R. A. (1997) Wavelets and field forecast verification. \emph{Mon. Wea. Rev.}, \bold{125}, 1329--1341.

Casati, B., Ross, G. and Stephenson, D. B. (2004) A new intensity-scale approach for the verification of spatial precipitation forecasts. \emph{Meteorol. Appl.} \bold{11}, 141--154.

Casati, B. (2010) New Developments of the Intensity-Scale Technique within the Spatial Verification Methods Inter-Comparison Project. \emph{Wea. Forecasting} \bold{25}, (1), 113--143, doi:10.1175/2009WAF2222257.1.

Lindsay, R. W., Percival, D. B. and Rothrock, D. A. (1996)  The discrete wavelet transform and the scale analysis of the surface properties of sea ice.  \emph{IEEE Transactions on Geoscience and Remote Sensing}, \bold{34} (3), 771--787.

Percival, D. B. and Guttorp, P. (1994)  Long-memory processes, the Allan variance and wavelets.  In \emph{Wavelets in Geophysics}, Foufoula-Georgiou, E. and Kumar, P., Eds., New York: Academic, pp. 325--343.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{dwt.2d}}, \code{\link{idwt.2d}}, \code{\link{hoods2d}}
}
\examples{
grid<- list( x= seq( 0,5,,64), y= seq(0,5,,64))
obj<-Exp.image.cov( grid=grid, theta=.5, setup=TRUE)
look<- sim.rf( obj) 
look[ look < 0] <- 0
look <- zapsmall( look)

look2 <- sim.rf( obj) 
look2[ look2 < 0] <- 0 
look2 <- zapsmall( look2)

res <- waverify2d(look, Y=look2)
plot(res)
summary(res)

\dontrun{
data( "UKobs6" )
data( "UKfcst6" )

look <- waverify2d(UKobs6, Y=UKfcst6)

plot(look, which.plots="energy")
look2 <- mowaverify2d(UKobs6, UKfcst6, J=8)
plot(look2, which.plots="energy")

plot(look, main1="NIMROD Analysis", main2="NIMROD Forecast")

plot(look2, main1="NIMROD Analysis", main2="NIMROD Forecast")

# Alternative using "SpatialVx" object.
data( "UKloc" )

hold <- make.SpatialVx( UKobs6, UKfcst6, loc = UKloc,
    map = TRUE, field.type = "Rainfall", units = "mm/h",
    data.name = "Nimrod", obs.name = "Obs 6",
    model.name = "Fcst 6" )

look <- waverify2d(hold)

plot(look, which.plots="details")


}
}
\keyword{ math }
