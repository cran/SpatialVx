\name{make.SpatialVx}
\alias{make.SpatialVx}
\alias{hist.SpatialVx}
\alias{plot.SpatialVx}
\alias{print.SpatialVx}
\alias{summary.SpatialVx}
\title{
Spatial Verification Sets -- SpatialVx Object
}
\description{
A list object containing the verification sets of spatial verification and forecast fields with pertinent information.
}
\usage{
make.SpatialVx(X, Xhat, thresholds = NULL, loc = NULL, projection =
                 FALSE, subset = NULL, time.vals = NULL, reg.grid =
                 TRUE, map = FALSE, loc.byrow = FALSE, field.type = "",
                 units = "", data.name = "", obs.name = "X", model.name
                 = "Xhat", q = c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75,
                 0.9, 0.95), qs = NULL)

\method{hist}{SpatialVx}(x, ..., time.point = 1, obs = 1, model = 1,
                 threshold.num = NULL)

\method{plot}{SpatialVx}( x, ..., time.point = 1, obs = 1, model = 1,
    col = c( "gray", tim.colors( 64 ) ), zlim, mfrow = c(1, 2) )

\method{print}{SpatialVx}(x, ...)

\method{summary}{SpatialVx}(object, ...)
}
\arguments{
  \item{X}{An n X m  matrix or n X m X T array giving the verification field of interest.  If an array, T is the number of time points.
}
  \item{Xhat}{
An n X m matrix or n X m X T array giving the forecast field of interest, or a list of such matrices/arrays with each component of the list an n X m matrix or n X m X T array defining a separate forecast model.
}
  \item{thresholds}{
single numeric, numeric vector, or Nu X Nf matrix, where Nu are the number of thresholds and Nf the number of forecast models plus one (for the verification) giving the threshold values of interest for the verification set or components of the set.  If NULL (default), then thresholds will be calculated as the quantiles (defined through argument \code{q}) of each field.  If a single numeric or a numeric vector, then an n X 2 matrix will be created (with column names \dQuote{X} and \dQuote{Xhat} where each column is identical.  Otherwise, different thresholds may be applied to each of the verification and forecast fields.  For example, if quantiles are used for thresholds, then each field will have their own unique thresholds.
}
  \item{loc}{
If lon/lat coordinates are available, then this is an n * m X 2 matrix giving the lon/lat coordinates of each grid point or location.  Should follow the convention used by the \pkg{maps} package.
}
  \item{projection}{
logical, are the grids projections onto the globe?  If so, when plotting, it will be attempted to account for this by using the \code{poly.image} function from package \code{fields}.  In this case, each column of \code{loc} will be converted to a matrix using \code{byrow} equal to the value of \code{loc.byrow}.
}
  \item{subset}{
vector identifying which specific grid points should be included (if not all of them).  This argument may be ignored by most functions and is included for possible future functionality.
}
  \item{time.vals}{
If more than one time point is available in the set (i.e., the set is of n X m X T arrays, with T > 1), then this argument can be used to define the time points.  If missing, the default will yield the vector \code{1:T}.  But, it is possible to include actual time information.  This is also a forward looking feature that may or may not have any subsequent functionality.
}
  \item{reg.grid}{
logical, is the verification set on a regular grid?  This is another feature intended for possible future functionality.  Most functions in this package assume the set is on a regular grid.
}
  \item{map}{
logical, should the plot function attempt to place a map onto the plot?  Only possible if the \code{loc} argument is given.
}
  \item{field.type, units}{
character used for plot labelling and printing information to the screen.  Describes what variable and in what units the field represents.
}
  \item{data.name, obs.name, model.name}{
character vector describing the verification set overall, the observation(s) and the model(s), resp.
}
  \item{q}{
numeric vector giving the values of quantiles to be used for thresholds.  Only used if \code{thresholds} is NULL.
}
  \item{qs}{
character vector describing the quantiles used.  Again, only used if \code{thresholds} is NULL.  This is for subsequent plot/print labelling purposes.
}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.  May also be a function name, in which case the function is applied at each grid point individually across time.  }
  \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
   \item{col, zlim}{optional arguments to \code{image}, and/or \pkg{fields} functions \code{poly.image} and \code{image.plot}}
  \item{mfrow}{ optional argument to change the mfrow argument for the graphic device.  Default is one row with two plots (obs and model).  If null, then the mfrow argument will not be changed. }
    \item{x, object}{list object of class \dQuote{SpatialVx}.}
    \item{loc.byrow}{logical determining whether to set up the location matrices using \code{byrow} = TRUE or FALSE (for use with \code{poly.image}.}
  \item{threshold.num}{If not null, then the threshold index to apply a threshold to the fields before creating the histogram.}
\item{\dots}{\code{hist} method: optional arguments to \code{hist}.

\code{plot} method: if \code{time.point} is a function, then these allow for optional arguments to this function to be passed.

\code{print} and \code{summary} methods: Not used.
}
}
\details{
This function merely describes a spatial verification set that includes the actual data as well as numerous attributes that are used by several of the subsequent functions that might be employed.  In many cases, the attribute information may be passed on to output from other functions for plot labelling and printing purposes (e.g., in order to identify the verification set, time point(s), etc.).

All (or perhaps most) subsequent functions in this package utilize objects of this class and the information contained in the attributes.  This function simply gathers information and data sets into a particular form.

The plot method function attempts to create an image plot of each field in the set (at each time point).  If projection is TRUE, then it will attempt to preserve the projection (via \code{poly.image} of package \pkg{fields}).  It will also add white contour lines showing the thresholds.  If map is TRUE and \code{loc} was supplied, then a map will also be added, if possible.
}
\value{
A list object with two (unnamed) components:
  \item{1 }{matrix or array (same as input argument) giving the observation}
  \item{2 }{Either a matrix or array (same as input argument) or a list of such objects if more than one forecast model.}

Several attributes are also included among the following:

   \item{xdim}{numeric of length 2 or 3 giving the dimensions of the verification set (i.e., m, n and T, if relevant).}
   \item{time}{vector giving the time values}
   \item{thresholds}{matrix giving the thresholds for each field.  If there is more than one forecast, and they use the same threshold, this matrix may have only two columns.}
   \item{udim}{the dimensions of the thresholds matrix.}
   \item{loc}{nm X 2 matrix giving the locations.  If loc was not given, this will be c(rep(1:n, m), rep(1:m, each=n)).}
   \item{subset}{If given, this is a numeric vector describing a subset of loc to be used.}
   \item{data.name, obs.name, model.name}{character vector giving the names of the data sets (same as input arguments).}
   \item{nforecast}{single numeric giving the number of different forecast models contained in the object.}
   \item{field.type, units}{character strings, same as input arguments.}
   \item{projection}{logical, is the grid a projection?}
   \item{reg.grid}{logical, is the grid a regular grid?}
   \item{map}{logical, should a map be added to image plots of the data?}
   \item{qs}{character vector giving the names of the threshold quantiles.}
   \item{msg}{A message involving the data name, field type and units for adding info to plots, etc.}
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{hoods2d}}, \code{\link{poly.image}}
}
\examples{
data( "UKobs6" )
data( "UKfcst6" )
data( "UKloc" )

hold <- make.SpatialVx( UKobs6, UKfcst6, thresholds = c(0.01, 20.01),
    loc = UKloc, field.type = "Precipitation", units = "mm/h",
    data.name = "Nimrod", obs.name = "Observations 6", model.name = "Forecast 6",
    map = TRUE)

hold

plot( hold )

hist( hold )

hist( hold, threshold.num = 2 )

}
\keyword{ manip }
\keyword{ hplot }
