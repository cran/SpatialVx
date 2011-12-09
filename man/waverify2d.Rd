\name{waverify2d}
\alias{waverify2d}
\alias{mowaverify2d}
\alias{plot.waverify2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High-Resolution Gridded Forecast Verification Using Discrete Wavelet Decomposition
}
\description{
High-resolution gridded forecast verification using discrete wavelet decomposition.
}
\usage{
waverify2d(X, Y, Clim = NULL, wavelet.type = "haar", J = NULL, useLL = FALSE, compute.shannon = FALSE, which.space = "field", verbose = FALSE)
mowaverify2d(X, Y, Clim = NULL, wavelet.type = "haar", J = 4, useLL = FALSE, compute.shannon = FALSE, which.space = "field", verbose = FALSE)
\method{plot}{waverify2d}(x, main1 = "X", main2 = "Y", main3 = "Climate", which.plots = c("all", "dwt2d", "details", "energy", "mse", "rmse", "acc"), separate = FALSE, ...)
}
\arguments{
  \item{X,Y,Clim}{
   m X n dyadic matrices (i.e., m = 2^M and n = 2^N, for M, N some integers) giving the verification and forecast fields (and optionally a climatology field), resp.
}
   \item{x}{
list object of class "waverify2d" as returned by \code{waverify2d}.
}
  \item{wavelet.type}{
character naming the type of wavelet to be used.  This is given as teh \code{wf} argument to the \code{dwt.2d} function of package \pkg{waveslim}.
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
  \item{which.space}{character (one of "field" or "wavelet") naming from which space the detail fields should be used.  If "field", then it is in the original field (or image) space (i.e., the detail reconstruction), and if "wavelet", it will be done in the wavelet space (i.e., the detail wavelet coefficients).}
  \item{main1,main2,main3}{optional characters naming each field to be used for the detail field plots and legend labelling on the energy plot.}
   \item{which.plots}{character vector describing which}
  \item{separate}{logical, should the plots be on their own devices (TRUE) or should some of them be put onto a single multi-panel device (FALSE, default)?}
  \item{verbose}{logical, should progress information be printed to the screen, including total run time?
}
\item{\dots}{
optional additonal plot or image.plot parameters.  If detail and energy, mse, rmse or acc plots are desired, must be applicable to both types of plots.
}
}
\details{
This is a function to use discrete wavelet decomposition to analyze verification sets along the lines of Briggs and Levine (1997), as well as Casati et al. (2004) and Casati (2009).  In the originally proposed formulation of Briggs and Levine (1997), continuous verification statistics (namely, the anomaly correlation coefficient (ACC) and root mean square error (RMSE)) are calculated for detail fields obtained from wavelet decompositions of each of a forecast and verification field (and for ACC a climatology field as well).  Casati et al. (2004) introduced an intensity scale approach that applies 2-d DWT to binary (obtained from thresholding) difference fields (Forecast - Verification), and applying a skill score at each level based on the mean square error (MSE).  Casati (2009) extended this idea to look at the energy at each level as well.

This function makes use of the dwt.2d and idwt.2d functions from package waveslim, and plot.waverify2d uses the plot.dwt.2d function if dwt2d is selected through the which.plots argument.  See the help file for these functions, the references therein and the references herein for more on these approaches.

Generally, it is not necessary to use the father wavelet for the detail fields, but for some purposes, it may be desired.

mowaverify2d is very similar to waverify2d, but it allows fields to be non-dyadic (and may subsequently be slower).  It uses the modwt.2d and imodwt.2d functions from the package waveslim.  In particular, it performs a maximal overlap discrete wavelet transform on a matrix of arbitrary dimension.  See the help file and references therein for modwt.2d for more information, as well as Percival and Guttorp (1994) and Lindsay et al. (1996).

In Briggs and Levine (1997), they state that the calculations can be done in either the data (called field here) space or the wavelet space, and they do their examples in the field space.  If the wavelets are orthogonal, then the detail coefficeints (wavelet space), can be analyzed with the assumption that they are independent; whereas in the data space, they typically cannot be assumed to be independent.  Therefore, most statistical tests should be performed in the wavelet space to avoid issues arising from spatial dependence.
}
\value{
A list object of class "waverify2d" with components:
  \item{J }{single numeric giving the number of levels.}
  \item{X.wave, Y.wave, Clim.wave }{objects of class "dwt.2d" describing the wavelet decompositions for the verification and forecast fields (and climatology, if applicable), resp. (see the help file for dwt.2d from package waveslim for more about these objects).}
   \item{Shannon.entropy}{numeric matrix giving the Shannon entropy for each field.}
   \item{energy}{numeric matrix giving the energy at each level and field.}
   \item{mse,rmse}{numeric vectors of length J giving the MSE/RMSE for each level between the verification and forecast fields.}
   \item{acc}{If a climatology field is supplied, this is a numeric vector giving the ACC for each level.}
}
\references{
Briggs, W. M. and R. A. Levine, 1997. Wavelets and field forecast verification. \emph{Mon. Wea. Rev.}, \bold{125}, 1329--1341.

Casati B, G Ross, and DB Stephenson, 2004. A new intensity-scale approach for the verification of spatial precipitation forecasts. \emph{Meteorol. Appl.} \bold{11}, 141--154.

Casati, B., 2010: New Developments of the Intensity-Scale Technique within the Spatial Verification Methods Inter-Comparison Project. \emph{Wea. Forecasting} \bold{25}, (1), 113--143, DOI: 10.1175/2009WAF2222257.1.

Lindsay, R. W., D. B. Percival, and D. A. Rothrock, 1996.  The discrete wavelet transform and the scale analysis of the surface properties of sea ice.  \emph{IEEE Transactions on Geoscience and Remote Sensing}, \bold{34} (3), 771--787.

Percival, D. B. and P. Guttorp, 1994.  Long-memory processes, the Allan variance and wavelets.  In \emph{Wavelets in Geophysics}, E. Foufoula-Georgiou and P. Kumar, Eds., New York: Academic, pp. 325--343.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{dwt.2d}}, \code{\link{idwt.2d}}, \code{\link{hoods2d}}, \code{\link{hoods2dPrep}}
}
\examples{
data(UKobs6)
data(UKfcst6)
look <- waverify2d(UKobs6, UKfcst6)
plot(look, which.plots="energy")
look2 <- mowaverify2d(UKobs6, UKfcst6, J=8)
plot(look2, which.plots="energy")

\dontrun{
pdf("dyadicDWTex.pdf")
plot(look, main1="NIMROD Analysis", main2="NIMROD Forecast")
dev.off()

pdf("nondyadicMODWTex.pdf")
plot(look2, main1="NIMROD Analysis", main2="NIMROD Forecast")
dev.off()

data(pert000)
data(pert004)
look <- mowaverify(pert000, pert004, J=8, verbose=TRUE) # Slow, but does not require fields to be dyadic.
plot(look, which.plots="energy") # Also can just do plot(look), but should print to a pdf file (e.g., using pdf()).

# Try one with some kind of climatology field.  Here using surrogater2d function.
data(UKloc)
hold <- surrogater2d(UKobs6, n=1, maxiter=50, verbose=TRUE)
hold <- matrix(hold, 256, 256)
image(hold, col=c("grey",tim.colors(64)), axes=FALSE)
image.plot(UKloc, col=c("grey",tim.colors(64)), legend.only=TRUE, horizontal=TRUE)

look <- waverify2d(UKobs6, UKfcst6, hold)
pdf("waveletEx.pdf")
plot(look)
dev.off()
}
}
\keyword{ math }
