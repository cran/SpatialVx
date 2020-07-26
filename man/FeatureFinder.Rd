\name{FeatureFinder}
\alias{FeatureFinder}
\alias{plot.features}
\alias{print.features}
\alias{summary.features}
\alias{plot.summary.features}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Threshold-based Feature Finder
}
\description{
Identify spatial features within a verification set using a threshold-based method.
}
\usage{
FeatureFinder(object, smoothfun = "disk2dsmooth", do.smooth = TRUE,
    smoothpar = 1, smoothfunargs = NULL, thresh = 1e-08, idfun = "disjointer",
    min.size = 1, max.size = Inf, fac = 1, zero.down = FALSE, time.point = 1,
    obs = 1, model = 1, ...)

\method{plot}{features}(x, ..., type = c("both", "obs", "model"))

\method{print}{features}(x, ...)

\method{summary}{features}(object, ...)

\method{plot}{summary.features}(x, ...)

}
\arguments{
  \item{object}{
An object of class \dQuote{SpatialVx}.
}
  \item{x}{list object of class \dQuote{features} as returned by \code{FeatureFinder}.}
  \item{smoothfun}{
character naming a 2-d smoothing function from package \pkg{smoothie}.  Not used if \code{do.smooth} is FALSE.
}
  \item{do.smooth}{
logical, should the field first be smoothed before trying to identify features (resulting field will not be smoothed, this is just for identifying features).  Default is to do convolution smoothing using a disc kernel as is recommended by Davis et al (2006a).
}
  \item{smoothpar}{
numeric of length one or two giving the smoothing parameter for \code{smoothfun}.  If length is two, the first value is applied to the forecast field and the second to the verification field.  The default smooth function (\code{smoothfun} argument) is the disk kernel smoother, so this argument gives the radius of the disk.  See the help file for \code{hoods2dsmooth} from package \pkg{smoothie} for other smoother choices; this argument corresponds to the \code{lambda} argument of these functions.}

  \item{smoothfunargs}{
list object with named additional arguments to \code{smoothfun}.
}
  \item{thresh}{
numeric vector of length one or two giving the threshold over which (inclusive) features should be identified.  If different thresholds are used for the forecast and verification fields, then the first element is the threshold for the forecast, and the second for the verification field.
}
  \item{idfun}{
character naming the function used to identify (and label) individual features in the thresholded, and possibly smoothed, fields.  Must take an argument 'x', the thresholded, and possibly smoothed, field.
}
  \item{min.size}{
numeric of length one or two giving the minimum number of contiguous grid points exceeding the threshold in order to be included as a feature (can be used to exclude any small features).  Default does not exclude any features.  If length is two, first value applies to the forecast and second to the verification field.
}
  \item{max.size}{
numeric of length one or two giving the maximum number of contiguous grid points exceeding the threshold in order to be included as a feature (can be used to exclude large features, if the need be).  Default does not exclude any features.  If length is two, then the first value applies ot the forecast field, and the second to the verification.
}
  \item{fac}{
numeric of length one or two giving a factor by which to multiply the R quantile in determining the threshold from the fields.  For example, ~ 1/15 is suggested in Wernli et al (2008, 2009).  If length is two, then the first value applies to the threshold of the forecast and the second to that of the verification field.
}
  \item{zero.down}{
logical, should negative values and relatively very small values be set to zero after smoothing the fields?  For thresholds larger than such values, this argument is moot.  'zapsmall' is used to set the very small positive values to zero.
}
  \item{time.point}{
numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.
}
  \item{obs, model}{
numeric indicating which observation/forecast model to select for the analysis.
}
\item{type}{
character string stating which features to plot (observed, forecast or both).  If both, a panel of two plots will be made side-by-side.
}
  \item{\dots}{
\code{FeatureFinder}: additional arguments to \code{idfun}.

\code{plot}: optional arguments to \code{image.plot} for the color legend bar only.

Not used by the \code{print} or \code{summary} functions.

The 'summary' method function can take the argument: 'silent'-logical, should information be printed to the screen (FALSE) or not (TRUE).
}
}
\details{
\code{FeatureFinder} applies for finding features based on three proposed methods from different papers; and also allows for combinations of the methods.  The methods include: the convolution-threshold approach of Davis et al. (2006a,b), which uses a disc kernel convolution smoother to first smooth the fields, then applies a threshold to remove low-intensity areas.  Feautres are identified by groups of contiguous \dQuote{events} (or connected components in the computer vision/image analysis literature) using \code{idfun}.  Nachamkin (2009) and Lack et al (2010) further require that features have at least \code{min.size} connected components in order to be considered a feature (in order to remove very small areas of threshold excesses).  Wernli et al. (2009) modify the threshold by a factor (see the \code{fac} argument).

In addition to the above options, it is also possible to remove features that are too large, as for some purposes, it is the small-scale features that are of interest, and sometimes the larger features can cause problems when merging and matching features across fields.
}
\value{
\code{FeatureFinder} returns a list object of class \dQuote{features} with comopnents:

  \item{data.name}{character vector naming the verification and forecast (R object) fields, resp.}
  \item{X.feats,Y.feats}{The identified features for the verification and forecast fields as returned by the idfun function.}
  \item{X.labeled,Y.labeled}{matrices of same dimension as the forecast and verification fields giving the images of the convolved and thresholded verification and forecast fields, but with each individually identified object labeled 1 to the number of objects in each field.}
  \item{identifier.function,identifier.label}{character strings naming the function and giving the long name (for use with plot method function).}

An additional attribute, named \dQuote{call}, is given.  This attribute shows the original function call, and is used mainly by the print function..

The plot method functions do not return anything.

The summary method function for objects of class \dQuote{features} returns a list with components:

\item{X,Y}{matrices whose rows are objects and columns are properties: centroidX and centroidY (the x- and y- coordinates for the feature centroids), area (the area of each feature in squared grid points), the orientation angle for the fitted major axis, the aspect ratio, Intensity0.25 and Intensity0.9 (the lower quartile and 0.9 quantile of intensity values for each feature).}

}
\references{
Davis, C. A., Brown, B. G. and Bullock, R. G. (2006a) Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. _Mon. Wea. Rev._, *134*, 1772-1784.

Davis, C. A., Brown, B. G. and Bullock, R. G. (2006b) Object-based verification of precipitation forecasts, Part II: Application to convective rain systems. _Mon. Wea. Rev._, *134*, 1785-1795.

Lack, S. A., Limpert, G. L. and Fox, N. I. (2010) An object-oriented multiscale verification scheme.  _Wea.  Forecasting_, *25*, 79-92, doi:10.1175/2009WAF2222245.1.

Nachamkin, J. E. (2009) Application of the composite method to the spatial forecast verification methods intercomparison dataset.  _Wea. Forecasting_, *24*, 1390-1400, doi:10.1175/2009WAF2222225.1.

Wernli, H., Paulat, M. Hagen, M. and Frei, C. (2008) SAL-A novel quality measure for the verification of quantitative precipitation forecasts.  _Mon. Wea. Rev._, *136*, 4470-4487.

Wernli, H., Hofmann, C. and Zimmer, M. (2009) Spatial forecast verification methods intercomparison project: Application of the SAL technique.  _Wea. Forecasting_, *24*, 1472-1484, doi:10.1175/2009WAF2222271.1.
}
\author{
Eric Gilleland
}
\note{
This function replaces the now deprecated functions: \code{convthresh}, \code{threshsizer} and \code{threshfac}.
}

\seealso{

Functions used in identifying the features (mostly from package \pkg{spatstat}:

\code{\link{connected}}, \code{\link{as.im}}, \code{\link{tess}}, \code{\link{tiles}}, \code{\link{owin}}, \code{\link{make.SpatialVx}}, \code{\link{disjointer}}

Functions that work on the resulting \dQuote{features} objects for merging and/or matching features within/across fields:

\code{\link{centmatch}}, \code{\link{deltamm}}, \code{\link{minboundmatch}}

To force merges (implicit or otherwise; recommended): \code{\link{MergeForce}}
}
\examples{

data( "ExampleSpatialVxSet" )

x <- ExampleSpatialVxSet$vx
xhat <- ExampleSpatialVxSet$fcst

hold <- make.SpatialVx( x, xhat, field.type = "simulated",
		       units = "none", data.name = "Example",
		       obs.name = "x", model.name = "xhat" )

look <- FeatureFinder( hold, smoothpar = 0.5, thresh = 1 )

par( mfrow=c(1,2))
image.plot(look$X.labeled)
image.plot(look$Y.labeled)

\dontrun{
x <- y <- matrix(0, 100, 100)
x[2:3,c(3:6, 8:10)] <- 1
y[c(4:7, 9:10),c(7:9, 11:12)] <- 1
     
x[30:50,45:65] <- 1
y[c(22:24, 99:100),c(50:52, 99:100)] <- 1
     
hold <- make.SpatialVx( x, y, field.type = "contrived", units = "none",
         data.name = "Example", obs.name = "x", model.name = "y" )
     
look <- FeatureFinder(hold, smoothpar=0.5) 

par( mfrow=c(1,2))
image.plot(look$X.labeled)
image.plot(look$Y.labeled)
     
look2 <- centmatch(look)
     
FeatureTable(look2)
     
look3 <- deltamm( look, N = 201, verbose = TRUE ) 
FeatureTable( look3 )


# data( "pert000" )
# data( "pert004" )
# data( "ICPg240Locs" )
     
# hold <- make.SpatialVx( pert000, pert004,
#    loc = ICPg240Locs, projection = TRUE, map = TRUE, loc.byrow = TRUE,
#    field.type = "Precipitation", units = "mm/h",
#    data.name = "ICP Perturbed Cases", obs.name = "pert000",
#    model.name = "pert004" )
     
# look <- FeatureFinder(hold, smoothpar=10.5, thresh = 5)
# plot(look)

# look2 <- deltamm( look, N = 701, verbose = TRUE )

# look2 <- MergeForce( look2 )

# plot(look2)

# summary( look2 )

# Now remove smallest features ( those with fewer than 700 grid squares).

# look <- FeatureFinder( hold, smoothpar = 10.5, thresh = 5, min.size = 700 )

# look # Now only two features.

# plot( look )

# Now remove the largest features (those with more than 1000 grid squares). 

# look <- FeatureFinder( hold, smoothpar = 10.5, thresh = 5, max.size = 1000 )

# look

# plot( look )

# Remove any features smaller than 700 and larger than 2000 grid squares).

# look <- FeatureFinder( hold, smoothpar = 10.5, thresh = 5,
 #    min.size = 700, max.size = 2000 )

# look

# plot( look )

# Find features according to Wernli et al. (2008).
# look <- FeatureFinder( hold, thresh = 5, do.smooth = FALSE, fac = 1 / 15 )

# look

# plot( look )

# Now do a mix of the two types of methods.
# look <- FeatureFinder( hold, smoothpar = 10.5, thresh = 5, fac = 1 / 15 )

# look

# plot( look )

}
}
\keyword{ math }
\keyword{ manip }
\keyword{ graphs }
