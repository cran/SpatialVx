\name{FeatureSuite}
\alias{FeatureSuite}
\alias{FeatureFunPrep}
\alias{plot.FeatureSuite}
\alias{print.FeatureSuite}
\alias{summary.FeatureSuite}
\title{
Wrapper function to perform feature-based spatial verification
}
\description{
Wrapper function to perform feature-based spatial verification.
}
\usage{
FeatureSuite(object, fun.object, time.point = 1, model = 1, verbose = FALSE)

FeatureFunPrep(identfun = NULL, identfun.args = NULL, mergefun = NULL,
    mergefun.args = NULL, matchfun = NULL, matchfun.args = NULL,
    analysisfun = NULL, analysisfun.args = NULL)

\method{plot}{FeatureSuite}(x, ...)

\method{print}{FeatureSuite}(x, ...)

\method{summary}{FeatureSuite}(object, ...)
}
\arguments{
  \item{object}{
\code{FeatureSuite}: list object of class \dQuote{SpatialVx}.

\code{summary}: list object of class \dQuote{FeatureSuite}.
}
  \item{fun.object}{
list object as output from the \code{FeatureFunPrep} function.
}
   \item{x}{list object of class \dQuote{FeatureSuite}.}
   \item{identfun,mergefun,matchfun,analysisfun}{character specifying a function for performing feature: (i) identification, (ii) merging, (iii) matching and (iv) analyzing the resulting (possibly merged and/or matched) features.  The arguments \code{mergefun} and \code{matchfun} may be left NULL, in which case no merging and/or matching of features will be carried out. }
   \item{identfun.args,mergefun.args,matchfun.args,analysisfun.args}{list objects containing relevant optional arguments to \code{identfun}, \code{mergefun}, \code{matchfun}, and \code{analysisfun}.  these may be left NULL.  See Details section for specific arguments that each function must take.  Such arguments should not be supplied here.}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{model}{numeric indicating which forecast model to select for the analysis.}
  \item{verbose}{
logical, should progress information be printed to the screen (including total lapsed time)?
}
  \item{\dots}{optional arguments to the \code{summary} method function.  This function calls \code{summary} method functions for objects of class feature and matched (both of which can be merged and/or matched), so these optional arguments may be any of those that are optional for those method functions.  In particular, \code{silent} (logical, should the information be printed to the screen or not?), and, for example, see the help file for \code{summary.FeatureMatchAnalyzer} for optional arguments if \code{FeatureMatchAnalyzer} is used for analyzing matched features.}
}
\details{
Several methods for spatial forecast verification have focused on identifying features within the fields of the verification set, and then analyzing them usually after merging and matching features (Ebert and McBride, 2000; Davis et al., 2006a, 2006b; Micheas et al., 2007; Wernli et al., 2008; Lack et al., 2010).  These functions are very generic functions that facilitate feature-based verification.  Some functions for performing specific methods are supplied such as, \code{convthresh}.

\code{Fcst.name} and \code{Vx.name} should point to gridded fields that can be coerced to class \dQuote{im} of package \pkg{spatstat}.  The functions must take \code{object} (same as \code{prep.object} above) as an input argument, as well as the following arguments (even if they are simply ignored within the function).  Additonal arguments are allowed, and can be accessed via the \code{identfun.args}, \code{mergefun.args}, etc. arguments above.

\code{identfun}: object (same as \code{prep.object} above) is the only required input argument.

\code{mergefun}: If supplied, must take x (output from \code{identfun}).

\code{matchfun}: If supplied, must take x (output from \code{identfun}), y (the output from \code{mergefun}).

\code{analysisfun}: x (output from \code{identfun}), y (output from \code{mergefun}), matches (output from \code{matchfun}).

The values returned by each of the above functions may vary substantially provided they all work together.

\code{summary.FeatureMatchAnalyzer} can be used to determine total interest values for matched features.
}
\value{
FeatureSuite returns a list object with numbered components depending on the numbers of thresholds (if any) supplied with output from the analysisfun function as well as additional components:
\item{prep.obj}{character naming the prep object used.}
\item{fun.object}{character naming the function obect supplied to the fun.object argument.}

FeatureFunPrep returns a list object with components:

   \item{identify}{list object with components fun and args.}
   \item{merge}{list object with components fun and args.}
   \item{match}{list object with components fun and args.}
   \item{analysis}{list object with components fun and args.}

In each case, fun is a character naming a function, and args is a list object specifying any additional arguments beyond those required as stated in the Details section above.

summary.FeatureSuite invisibly returns a list object with possible components (depending on which functions are used):
\item{feature.summary}{returned value from the summary method function for objects of class \dQuote{features}.}
\item{analysis.summary}{returned value from the analysis function's summary method.  See summary.FeatureMatchAnalyzer for an example.}

}
\references{
Davis, C. A., Brown, B. G. and Bullock, R. G. (2006a) Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.

Davis, C. A., Brown, B. G. and Bullock, R. G. (2006b) Object-based verification of precipitation forecasts, Part II: Application to convective rain systems. \emph{Mon. Wea. Rev.}, \bold{134}, 1785--1795.

Ebert, E. E. and McBride, J. L. (2000) Verification of precipitation in weather systems: Determination of systematic errors.  \emph{J. Hydrology}, \bold{239}, 179--202.

Lack, S., Limpert, G. L. and Fox, N. I. (2010) An object-oriented multiscale verification scheme. \emph{Wea. Forecasting}, \bold{25}, 79--92, doi:10.1175/2009WAF2222245.1.

Micheas, A. C., Fox, N. I., Lack, S. A. and Wikle, C. K. (2007) Cell identification and verification of QPF ensembles using shape analysis techniques. \emph{J. Hydrology}, \bold{343}, 105--116.

Wernli, H., Paulat, M., Hagen, M. and Frei, C. (2008) SAL--A novel quality measure for the verification of quantitative precipitation forecasts.  \emph{Mon. Wea. Rev.}, \bold{136}, 4470--4487.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{connected}}, \code{\link{tiles}}, \code{\link{tess}}, \code{\link{as.im}}, \code{\link{convthresh}}, \code{\link{make.SpatialVx}}
}
\examples{
x <- y <- matrix(0, 10, 12)
x[2:3,c(3:6, 8:10)] <- 1
y[c(4:7, 9:10),c(7:9, 11:12)] <- 1

hold <- make.SpatialVx(x, y, field.type="contrived",
    units="units", data.name=c("Example", "x", "y"))

# Do SAL
holdfun <- FeatureFunPrep(identfun="threshfac", analysisfun="saller") 
look <- FeatureSuite(hold, holdfun)
summary(look)

\dontrun{
data(geom000)
data(geom002)
data(ICPg240Locs)

hold <- make.SpatialVx(geom000, geom002, projection=TRUE, map=TRUE,
    loc=ICPg240Locs, field.type="Precipitation", units="mm/h",
    data.name=c("Geometric ICP Cases", "geom000", "geom002"))

holdfun <- FeatureFunPrep(identfun="threshfac",
    identfun.args=list(thresh=0.1), analysisfun="saller")

look <- FeatureSuite(hold, holdfun)
summary(look)

## Try to make a good example!!
data(UKfcst6)
data(UKobs6)
data(UKloc)

hold <- make.SpatialVx(UKobs6, UKfcst6, map=TRUE, loc=UKloc,
    field.type="Rainfall", units="mm/h",
    data.name=c("Nimrod", "obs 6", "fcst 6"))

holdfun <- FeatureFunPrep(identfun="convthresh",
    identfun.args=list(smoothpar=1.5, thresh=1),
    matchfun="centmatch", analysisfun="FeatureMatchAnalyzer")

look <- FeatureSuite(hold, holdfun, verbose=TRUE)
summary(look)

# use the pdf function to print the following plots to a file
# in order to see all of them.

plot(look)
}
}
\keyword{ manip }
