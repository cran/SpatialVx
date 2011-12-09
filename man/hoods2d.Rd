\name{hoods2d}
\alias{hoods2d}
\alias{plot.hoods2d}
\title{
Calculate various neighborhood verification statistics for a gridded verification set.
}
\description{
Calculates most of the various neighborhood verification statistics for a gridded verification set as reviewed in Ebert (2008).
}
\usage{
hoods2d(obj, which.methods = c("mincvr", "multi.event", "fuzzy", "joint", "fss", "pragmatic"),
	verbose = FALSE)
\method{plot}{hoods2d}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{list object of class "hoods2dPrep" as returned by the \code{hoods2dPrep} function.}
  \item{which.methods}{character vector giving the names of the methods.  Default is for the entire list to be executed.  See Details section for specific option information.}
  \item{verbose}{logical, should progress information be printed to the screen?  Will also give the amount of time (in hours, minutes, or seconds) that the function took to run.}
  \item{x}{ list object output from \code{hoods2d}.}
  \item{\dots}{not used.}
}
\details{ \code{hoods2d} uses an object from the function 'hoods2dPrep' that includes most of the options utilized by this function, including the thresholds and neighborhood lengths (levels) to be used.  The neighborhood methods (cf. Ebert 2008, 2009; Gilleland et al., 2009, 2010) apply a (kernel) smoothing filter (cf. Hastie and Tibshirani, 1990) to either the raw forecast (and possibly also the observed) field(s) or to the binary counterpart(s) determined by thresholding.

The specific smoothing filter applied for these methods could be of any type, but those described in Ebert (2008) are generally taken to be "neighborhood" filters.  In some circles, this is referred to as a convolution filter with a boxcar kernel.  Because the smoothing filter can be represented this way, it is possible to use the convolution theorem with the Fast Fourier Transform (FFT) to perform the neighborhood smoothing operation very quickly. The particular approach used here "zero pads" the field, and replaces all missing values with zero as well, which is also the approach proposed in Roberts and Lean (2008).  If any missing values are introduced after the convolution, they are removed.

If zero-padding is undesirable, then two options are available:  1. Give a subset to the 'hoods2dPrep' function (e.g., some tile within the domain) so that the final statistics are calculated only on this subset, or  2. Extrapolate the fields before applying this function (and 'hoods2dPrep').  In the case of 2, you might want to also give it the subset (e.g., to give it only the original un-extrapolated fields).
  
To simplify the notation for the descriptions of the specific methods employed here, the notation of Ebert (2008) is adopted.  That is, if a method uses neighborhood smoothed observations (NO), then the neighborhood smoothed observed field is denoted <X>s, and the associated binary field, by <Ix>s.  Otherwise, if the observation field is not smoothed (denoted by SO in Ebert, 2008), then simply X or Ix are used.  Similarly, for the forecast field, <Y>s or <Iy>s are used for neighborhood smoothed forecast fields (NF).  If it is the binary fields that are smoothed (e.g., the original fields are thresholded before smoothing), then the resulting fields are denoted <Px>s and <Py>s, resp.  Below, NO-NF indicates that a neighborhood smoothed observed field (<Yx>s, <Ix>s, or <Px>s) is compared with a neighborhood smoothed forecast field, and SO-NF indicates that the observed field is not smoothed.

Options for 'which.methods' include:

"mincvr": (NO-NF) The minimum coverage method compares <Ix>s and <Iy>s by thresholding the neighborhood smoothed fields <Px>s and <Py>s (i.e., smoothed versions of Ix and Iy) to obtain <Ix>s and <Iy>s.  Indicator fields <Ix>s and <Iy>s are created by thresholding <Px>s and <Py>s by frequency threshold \code{Pe} given by the \code{obj} argument.  Scores calculated between <Ix>s and <Iy>s include: probability of detecting an event (pod, also known as the hit rate), false alarm ratio (far) and ets (cf. Ebert, 2008, 2009).

"multi.event": (SO-NF) The Multi-event Contingency Table method compares the binary observed field Ix against the smoothed forecast indicator field, <Iy>s, which is determined similarly as for "mincvr" (i.e., using Pe as a threshold on <Py>s).  The hit rate and false alarm rate (F) are calculated (cf. Atger, 2001).

"fuzzy": (NO-NF) The fuzzy logic approach compares <Px>s to <Py>s by creating a new contingency table where hits = sum_i min(<Px>s_i,<Py>s_i), misses = sum_i min(<Px>s_i,1-<Py>s_i), false alarms = sum_i min(1-<Px>s_i,<Py>s_i), and correct negatives = sum_i min(1-<Px>s_i,1-<Py>s_i) (cf. Ebert 2008).

"joint": (NO-NF) Similar to "fuzzy" above, but hits  = sum_i prod(<Px>s_i,<Py>s_i), misses = sum_i prod(<Px>s_i,1-<Py>s_i), false alarms = sum_i prod(1-<Px>s_i,<Py>s_i), and correct negatives = sum_i prod(1-<Px>s_i,1-<Py>s_i) (cf. Ebert, 2008).

"fss": (NO-NF) Compares <Px>s and <Py>s directly using a Fractions Brier and Fractions Skill Score (FBS and FSS, resp.), where FBS is the mean square difference between <Px>s and <Py>s, and the FSS is one minus the FBS divided by a reference MSE given by the sum of the sum of squares of <Px>s and <Py>s individually, divided by the total (cf. Roberts and Lean, 2008).

"pragmatic": (SO-NF) Compares Ix with <Py>s, calculating the Brier and Brier Skill Score (BS and BSS, resp.), where the reference forecast used for the BSS is taken to be the mean square error between the base rate and Ix (cf. Theis et al., 2005).

}
\value{
 A list object of class "hoods2d" with components determined by the \code{which.methods} argument.  Each component is itself a list object containing relevant components to the given method.  For example, hit rate is abbreviated pod here, and if this is an output for a method, then there will be a component named pod (all lower case).  The Gilbert Skill Score is abbreviated 'ets' (equitable threat score; again all lower case here).  The list components will be some or all of the following.
    \item{mincvr }{list with components: pod, far and ets}
    \item{multi.event }{list with components: pod, f and hk}
    \item{fuzzy}{list with components: pod, far and ets}
    \item{joint }{list with components: pod, far and ets}
    \item{fss}{list with components: fss, fss.uniform, fss.random}
    \item{pragmatic}{list with components: bs and bss}
}
\references{
Atger, F., 2001: Verification of intense precipitation forecasts from single models and ensemble prediction systems.  \emph{Nonlin. Proc. Geophys.}, \bold{8}, 401--417.

Ebert, E. E., 2008: Fuzzy verification of high resolution gridded forecasts: A review and proposed framework.  \emph{Meteorol. Appl.}, \bold{15}, 51--64. DOI: 10.1002/met.25 Available at \url{http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf}

Ebert, E. E., 2009: Neighborhood verification: A strategy for rewarding close forecasts.  \emph{Wea. Forecasting}, \bold{24}, 1498--1510, DOI: 10.1175/2009WAF2222251.1.

Gilleland, E., D. Ahijevych, B.G. Brown, B. Casati, and E.E. Ebert, 2009. Intercomparison of Spatial Forecast Verification Methods. \emph{Wea. Forecasting}, \bold{24}, 1416--1430, DOI: 10.1175/2009WAF2222269.1.

Gilleland, E., D.A. Ahijevych, B.G. Brown and E.E. Ebert, 2010: Verifying Forecasts Spatially. \emph{Bull. Amer. Meteor. Soc.}, October, 1365--1373.

Hastie, T. J. and Tibshirani, R. J., 1990: \emph{Generalized Additive Models}. Chapman \& Hall/CRC Monographs on Statistics and Applied Probability 43, 335pp.

Roberts, N. M. and H. W. Lean, 2008: Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events.  \emph{Mon. Wea. Rev.}, \bold{136}, 78--97. DOI: 10.1175/2007MWR2123.1.

Theis, S. E., A. Hense, U. Damrath, 2005: Probabilistic precipitation forecasts from a deterministic model: A pragmatic approach.  \emph{Meteorol. Appl.}, \bold{12}, 257--268.

Yates, E., S. Anquetin, V. Ducrocq, J.-D. Creutin, D. Ricard, and K. Chancibault, 2006: Point and areal validation of forecast precipitation fields.  \emph{Meteorol. Appl.}, \bold{13}, 1--20.

Zepeda-Arce, J., E. Foufoula-Georgiou, K. K. Droegemeier, 2000: Space-time rainfall organization and its role in validating quantitative precipitation forecasts.  \emph{J. Geophys. Res.}, \bold{105}(D8), 10,129--10,146.
}
\author{
Eric Gilleland
}
\note{Thresholded fields are taken to be >= the threshold.}

\seealso{\code{\link{fft}}, \code{\link{kernel2dsmooth}}, \code{\link{plot.hoods2d}}, \code{\link{vxstats}} }

\examples{
x <- y <- matrix( 0, 50, 50)
x[ sample(1:50,10), sample(1:50,10)] <- rexp( 100, 0.25)
y[ sample(1:50,20), sample(1:50,20)] <- rexp( 400)
hold <- hoods2dPrep("y", "x", thresholds=c(0.1, 0.5), levels=c(1, 3, 20))
look <- hoods2d( hold, which.methods=c("multi.event", "fss"))
look
\dontrun{
plot(look)

data(geom001)
data(geom000)
data(ICPg240Locs)
hold <- hoods2dPrep( "geom001", "geom000", thresholds=c(0.01,50.01), levels=c(1, 3, 9, 17, 33, 65, 129, 257),
	loc=ICPg240Locs, units="in/100")
look <- hoods2d(hold, verbose=TRUE)
plot( look) # Might want to use 'pdf' to print these out so that all of them can be observed.
data(pert004)
data(pert000)
hold <- hoods2dPrep( "pert004", "pert000", thresholds=c(1,10,50), levels=c(1, 3, 17, 33, 65, 129, 257),
	loc=ICPg240Locs, units="mm/h")
look <- hoods2d( hold, verbose=TRUE)
plot( look)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
