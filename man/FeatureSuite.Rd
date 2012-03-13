\name{FeatureSuite}
\alias{FeatureSuite}
\alias{FeatureSuitePrep}
\alias{FeatureFunPrep}
\title{
Wrapper function to perform feature-based spatial verification
}
\description{
Wrapper function to perform feature-based spatial verification.
}
\usage{
FeatureSuite(prep.object, fun.object, verbose = FALSE)
FeatureSuitePrep(Fcst.name, Vx.name, loc = NULL, units = NULL)
FeatureFunPrep(identfun = NULL, identfun.args = NULL, mergefun = NULL, mergefun.args = NULL, matchfun = NULL, matchfun.args = NULL, analysisfun = NULL, analysisfun.args = NULL)
}
\arguments{
  \item{prep.object}{
list object as output from the \code{FeatureSuitePrep} function.
}
  \item{fun.object}{
list object as output from the \code{FeatureFunPrep} function.
}
   \item{Fcst.name, Vx.name}{character giving the name of the fields in the verification set.}
   \item{loc}{(optional) numeric n*m X 2 matrix giving lon/lat coordinates for the verification set.  Currently, this is not used for anything, but may have a use in the future, e.g., for adding a map to plots of neighborhood smoothed fields.}
  \item{units}{(optional) character giving the name of the units of the verification set.}
   \item{identfun,mergefun,matchfun,analysisfun}{character specifying a function for performing feature: (i) identification, (ii) merging, (iii) matching and (iv) analyzing the resulting (possibly merged and/or matched) features.  The arguments \code{mergefun} and \code{matchfun} may be left NULL, in which case no merging and/or matching of features will be carried out. }
   \item{identfun.args,mergefun.args,matchfun.args,analysisfun.args}{list objects containing relevant optional arguments to \code{identfun}, \code{mergefun}, \code{matchfun}, and \code{analysisfun}.  these may be left NULL.  See Details section for specific arguments that each function must take.  Such arguments should not be supplied here.}
  \item{verbose}{
logical, should progress information be printed to the screen (including total lapsed time)?
}
}
\details{
Several methods for spatial forecast verification have focused on identifying features within the fields of the verification set, and then analyzing them usually after merging and matching features (Ebert and McBride, 2000; Davis et al., 2006a, 2006b; Micheas et al., 2007; Wernli et al., 2008; Lack et al., 2010).  These functions are very generic functions that facilitate feature-based verification.  Some functions for performing specific methods are supplied such as, convthresh.

Fcst.name and Vx.name should point to gridded fields that can be coerced to class "im" of package spatstat.  The functions must take object (same as prep.object above) as an input argument, as well as the following arguments (even if they are simply ignored within the function).  Additonal arguments are allowed, and can be accessed via the identfun.args, mergefun.args, etc. arguments above.

identfun: object (same as prep.object above) is the only required input argument.

mergefun: If supplied, must take x (output from identfun).

matchfun: If supplied, must take x (output from identfun), y (the output from mergefun).

analysisfun: x (output from identfun), y (output from mergefun), matches (output from matchfun).

The values returned by each of the above functions may vary substantially provided they all work together.

}
\value{
FeatureSuite returns a list object with numbered components depending on the numbers of thresholds (if any) supplied with output from the analysisfun function as well as additional components:
\item{prep.obj}{character naming the prep object used.}
\item{fun.object}{character naming the function obect supplied to the fun.object argument.}

FeatureSuitePrep returns a list object with components:
   \item{Fcst.name}{character giving the name of the forecast matrix of the verification set.}
   \item{Vx.name}{character giving the name of the observed matrix of the verification set.}
   \item{xdim}{numeric vector of length 2 giving the dimensions of the verification set.}
   \item{loc}{Same as in Arguments above.}
   \item{units}{Same as in Arguments above.}

FeatureFunPrep returns a list object with components:

   \item{identify}{list object with components fun and args.}
   \item{merge}{list object with components fun and args.}
   \item{match}{list object with components fun and args.}
   \item{analysis}{list object with components fun and args.}

In each case, fun is a character naming a function, and args is a list object specifying any additional arguments beyond those required as stated in the Details section above.

}
\references{
Davis CA, BG Brown, and RG Bullock, 2006a. Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.

Davis CA, BG Brown, and RG Bullock, 2006b. Object-based verification of precipitation forecasts, Part II: Application to convective rain systems. \emph{Mon. Wea. Rev.}, \bold{134}, 1785--1795.

Ebert, E. E. and J. L. McBride, 2000: Verification of precipitation in weather systems: Determination of systematic errors.  \emph{J. Hydrology}, \bold{239}, 179--202.

Lack, S., G.L. Limpert and N.I. Fox, 2010: An object-oriented multiscale verification scheme. \emph{Wea. Forecasting}, \bold{25}, 79--92, DOI: 10.1175/2009WAF2222245.1

Micheas AC, NI Fox, SA Lack, and CK Wikle. 2007. Cell identification and verification of QPF ensembles using shape analysis techniques. \emph{J. Hydrology}, \bold{343}, 105--116.

Wernli, H., M. Paulat, M. Hagen, and C. Frei, 2008: SAL--A novel quality measure for the verification of quantitative precipitation forecasts.  \emph{Mon. Wea. Rev.}, \bold{136}, 4470--4487.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{connected}}, \code{\link{tiles}}, \code{\link{tess}}, \code{\link{as.im}}, \code{\link{convthresh}}
}
\examples{
x <- y <- matrix(0, 10, 12)
x[2:3,c(3:6, 8:10)] <- 1
y[c(4:7, 9:10),c(7:9, 11:12)] <- 1

hold <- FeatureSuitePrep("y", "x")

# Do SAL
holdfun <- FeatureFunPrep(identfun="salIDfun", analysisfun="saller") 
look <- FeatureSuite( hold, holdfun)

\dontrun{
data(geom000)
data(geom002)
hold <- FeatureSuitePrep("geom002", "geom000")
holdfun <- FeatureFunPrep(identfun="salIDfun", identfun.args=list(thresh=0.1), analysisfun="saller")

look <- FeatureSuite(hold, holdfun)

## Try to make a good example!!
hold <- FeatureSuitePrep("pert004", "pert000")
holdfun <- FeatureFunPrep(identfun="convthresh", identfun.args=list(smoothpar=5.5, thresh=0.01), analysisfun="FeatureMatchAnalyzer")
look <- FeatureSuite(hold, holdfun)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
