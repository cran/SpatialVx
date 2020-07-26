\name{hoods2d}
\alias{hoods2d}
\alias{plot.hoods2d}
\alias{print.hoods2d}
\title{
Neighborhood Verification Statistics for a Gridded Verification Set.
}
\description{
Calculates most of the various neighborhood verification statistics for a gridded verification set as reviewed in Ebert (2008).
}
\usage{
hoods2d(object, which.methods = c("mincvr", "multi.event",
                 "fuzzy", "joint", "fss", "pragmatic"), time.point = 1,
                 obs = 1, model = 1, Pe = NULL, levels = NULL, max.n =
                 NULL, smooth.fun = "hoods2dsmooth", smooth.params =
                 NULL, rule = ">=", verbose = FALSE)

\method{plot}{hoods2d}(x, ..., add.text = FALSE)

\method{print}{hoods2d}(x, ...)
}
\arguments{
  \item{object}{list object of class \dQuote{SpatialVx}.}
  \item{which.methods}{character vector giving the names of the methods.  Default is for the entire list to be executed.  See Details section for specific option information.}
  \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
  \item{Pe}{(optional) numeric vector of length q >= 1 to be applied to the fields sPy and possibly sPx (see details section).  If NULL, then it is taken to be the most relaxed requirement (i.e., that an event occurs at least once in a neighborhood) of Pe=1/(nlen^2), where nlen is the length of the neighborhood.}
   \item{levels}{numeric vector giving the successive values of the smoothing parameter.  For example, for the default method, these are the neighborhood lengths over which the levels^2 nearest neighbors are averaged for each point.  Values should make sense for the specific smoothing function.  For example, for the default method, these should be odd integers.}
    \item{max.n}{(optional) single numeric giving the maximum neighborhood length to use.  Only used if levels are NULL.}
    \item{smooth.fun}{character giving the name of a smoothing function to be applied.  Default is an average over the n^2 nearest neighbors, where n is taken to be each value of the \code{levels} argument.}
    \item{smooth.params}{list object containing any optional arguments to \code{smooth.fun}.  Use NULL if none.}
  \item{rule}{character string giving the threshold rule to be applied.  See help file for \code{thresholder} function for more information.}
  \item{verbose}{logical, should progress information be printed to the screen?  Will also give the amount of time (in hours, minutes, or seconds) that the function took to run.}
  \item{x}{list object output from \code{hoods2d}.}
    \item{add.text}{logical, should the text values of FSS be added to its quilt plot?}
  \item{\dots}{not used.}
}
\details{
\code{hoods2d} uses an object of class \dQuote{SpatialVx} that includes some information utilized by this function, including the thresholds to be used.  The neighborhood methods (cf. Ebert 2008, 2009; Gilleland et al., 2009, 2010) apply a (kernel) smoothing filter (cf. Hastie and Tibshirani, 1990) to either the raw forecast (and possibly also the observed) field(s) or to the binary counterpart(s) determined by thresholding.

The specific smoothing filter applied for these methods could be of any type, but those described in Ebert (2008) are generally taken to be \dQuote{neighborhood} filters.  In some circles, this is referred to as a convolution filter with a boxcar kernel.  Because the smoothing filter can be represented this way, it is possible to use the convolution theorem with the Fast Fourier Transform (FFT) to perform the neighborhood smoothing operation very quickly. The particular approach used here \dQuote{zero pads} the field, and replaces all missing values with zero as well, which is also the approach proposed in Roberts and Lean (2008).  If any missing values are introduced after the convolution, they are removed.

To simplify the notation for the descriptions of the specific methods employed here, the notation of Ebert (2008) is adopted.  That is, if a method uses neighborhood smoothed observations (NO), then the neighborhood smoothed observed field is denoted <X>s, and the associated binary field, by <Ix>s.  Otherwise, if the observation field is not smoothed (denoted by SO in Ebert, 2008), then simply X or Ix are used.  Similarly, for the forecast field, <Y>s or <Iy>s are used for neighborhood smoothed forecast fields (NF).  If it is the binary fields that are smoothed (e.g., the original fields are thresholded before smoothing), then the resulting fields are denoted <Px>s and <Py>s, resp.  Below, NO-NF indicates that a neighborhood smoothed observed field (<Yx>s, <Ix>s, or <Px>s) is compared with a neighborhood smoothed forecast field, and SO-NF indicates that the observed field is not smoothed.

Options for which.methods include:

\dQuote{mincvr}: (NO-NF) The minimum coverage method compares <Ix>s and <Iy>s by thresholding the neighborhood smoothed fields <Px>s and <Py>s (i.e., smoothed versions of Ix and Iy) to obtain <Ix>s and <Iy>s.  Indicator fields <Ix>s and <Iy>s are created by thresholding <Px>s and <Py>s by frequency threshold \code{Pe} given by the \code{obj} argument.  Scores calculated between <Ix>s and <Iy>s include: probability of detecting an event (pod, also known as the hit rate), false alarm ratio (far) and ets (cf. Ebert, 2008, 2009).

\dQuote{multi.event}: (SO-NF) The Multi-event Contingency Table method compares the binary observed field Ix against the smoothed forecast indicator field, <Iy>s, which is determined similarly as for \dQuote{mincvr} (i.e., using Pe as a threshold on <Py>s).  The hit rate and false alarm rate (F) are calculated (cf. Atger, 2001).

\dQuote{fuzzy}: (NO-NF) The fuzzy logic approach compares <Px>s to <Py>s by creating a new contingency table where hits = sum_i min(<Px>s_i,<Py>s_i), misses = sum_i min(<Px>s_i,1-<Py>s_i), false alarms = sum_i min(1-<Px>s_i,<Py>s_i), and correct negatives = sum_i min(1-<Px>s_i,1-<Py>s_i) (cf. Ebert 2008).

\dQuote{joint}: (NO-NF) Similar to \dQuote{fuzzy} above, but hits  = sum_i prod(<Px>s_i,<Py>s_i), misses = sum_i prod(<Px>s_i,1-<Py>s_i), false alarms = sum_i prod(1-<Px>s_i,<Py>s_i), and correct negatives = sum_i prod(1-<Px>s_i,1-<Py>s_i) (cf. Ebert, 2008).

\dQuote{fss}: (NO-NF) Compares <Px>s and <Py>s directly using a Fractions Brier and Fractions Skill Score (FBS and FSS, resp.), where FBS is the mean square difference between <Px>s and <Py>s, and the FSS is one minus the FBS divided by a reference MSE given by the sum of the sum of squares of <Px>s and <Py>s individually, divided by the total (cf. Roberts and Lean, 2008).

\dQuote{pragmatic}: (SO-NF) Compares Ix with <Py>s, calculating the Brier and Brier Skill Score (BS and BSS, resp.), where the reference forecast used for the BSS is taken to be the mean square error between the base rate and Ix (cf. Theis et al., 2005).

}
\value{
 A list object of class \dQuote{hoods2d} with components determined by the \code{which.methods} argument.  Each component is itself a list object containing relevant components to the given method.  For example, hit rate is abbreviated pod here, and if this is an output for a method, then there will be a component named pod (all lower case).  The Gilbert Skill Score is abbreviated 'ets' (equitable threat score; again all lower case here).  The list components will be some or all of the following.
    \item{mincvr }{list with components: pod, far and ets}
    \item{multi.event }{list with components: pod, f and hk}
    \item{fuzzy}{list with components: pod, far and ets}
    \item{joint }{list with components: pod, far and ets}
    \item{fss}{list with components: fss, fss.uniform, fss.random}
    \item{pragmatic}{list with components: bs and bss}

New attributes are added giving the values for some of the optional arguments: levels, max.n, smooth.fun, smooth.params and Pe.
}
\references{
Atger, F. (2001) Verification of intense precipitation forecasts from single models and ensemble prediction systems.  \emph{Nonlin. Proc. Geophys.}, \bold{8}, 401--417.

Ebert, E. E. (2008) Fuzzy verification of high resolution gridded forecasts: A review and proposed framework.  \emph{Meteorol. Appl.}, \bold{15}, 51--64. doi:10.1002/met.25 

Ebert, E. E. (2009) Neighborhood verification: A strategy for rewarding close forecasts.  \emph{Wea. Forecasting}, \bold{24}, 1498--1510, doi:10.1175/2009WAF2222251.1.

Gilleland, E., Ahijevych, D., Brown, B. G., Casati, B. and Ebert, E. E. (2009) Intercomparison of Spatial Forecast Verification Methods. \emph{Wea. Forecasting}, \bold{24}, 1416--1430, doi:10.1175/2009WAF2222269.1.

Gilleland, E., Ahijevych, D. A., Brown, B. G. and Ebert, E. E. (2010) Verifying Forecasts Spatially. \emph{Bull. Amer. Meteor. Soc.}, October, 1365--1373.

Hastie, T. J. and Tibshirani, R. J. (1990) \emph{Generalized Additive Models}. Chapman \& Hall/CRC Monographs on Statistics and Applied Probability 43, 335pp.

Roberts, N. M. and Lean, H. W. (2008) Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events.  \emph{Mon. Wea. Rev.}, \bold{136}, 78--97. doi:10.1175/2007MWR2123.1.

Theis, S. E., Hense, A. Damrath, U. (2005) Probabilistic precipitation forecasts from a deterministic model: A pragmatic approach.  \emph{Meteorol. Appl.}, \bold{12}, 257--268.

Yates, E., Anquetin, S. Ducrocq, V., Creutin, J.-D., Ricard, D. and Chancibault, K. (2006) Point and areal validation of forecast precipitation fields.  \emph{Meteorol. Appl.}, \bold{13}, 1--20.

Zepeda-Arce, J., Foufoula-Georgiou, E., Droegemeier, K. K. (2000) Space-time rainfall organization and its role in validating quantitative precipitation forecasts.  \emph{J. Geophys. Res.}, \bold{105}(D8), 10,129--10,146.
}
\author{
Eric Gilleland
}
\note{Thresholded fields are taken to be >= the threshold.}

\seealso{\code{\link{fft}}, \code{\link{kernel2dsmooth}}, \code{\link{plot.hoods2d}}, \code{\link{vxstats}}, \code{\link{thresholder}} }

\examples{
x <- y <- matrix( 0, 50, 50)
x[ sample(1:50,10), sample(1:50,10)] <- rexp( 100, 0.25)
y[ sample(1:50,20), sample(1:50,20)] <- rexp( 400)
hold <- make.SpatialVx( x, y, thresholds = c(0.1, 0.5), field.type = "Random Exp. Var." )
look <- hoods2d( hold, which.methods=c("multi.event", "fss"), levels=c(1, 3, 19))
look
plot(look)

\dontrun{
data( "UKobs6" )
data( "UKfcst6" )
data( "UKloc" )

hold <- make.SpatialVx( UKobs6, UKfcst6, thresholds = c(0.01, 20.01),
    projection = TRUE, map = TRUE, loc = UKloc, loc.byrow = TRUE,
    field.type = "Precipitation", units = "mm/h",
    data.name = "Nimrod", obs.name = "Observations 6",
    model.name = "Forecast 6" )

hold
plot(hold)
hist(hold, col="darkblue")

look <- hoods2d(hold, which.methods=c("multi.event", "fss"), 
    levels=c(1, 3, 5, 9, 17), verbose=TRUE)
plot(look)

data( "geom001" )
data( "geom000" )
data( "ICPg240Locs" )

hold <- make.SpatialVx( geom000, geom001, thresholds = c(0.01, 50.01),
    projection = TRUE, map = TRUE, loc = ICPg240Locs, loc.byrow = TRUE,
    field.type = "Geometric Objects Pretending to be Precipitation",
    units = "mm/h", data.name = "ICP Geometric Cases", obs.name = "geom000",
    model.name = "geom001" )

look <- hoods2d(hold, levels=c(1, 3, 9, 17, 33, 65, 129, 257), verbose=TRUE)

plot( look) # Might want to use 'pdf' first.

    }
}
\keyword{ math }
