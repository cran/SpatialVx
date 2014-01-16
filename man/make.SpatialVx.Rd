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
make.SpatialVx(X, Xhat, thresholds = NULL,
    loc = NULL, projection = FALSE, subset = NULL,
    time.vals = NULL, reg.grid = TRUE, map = FALSE,
    loc.byrow = FALSE, field.type = "", units = "",
    data.name = c("X", "Xhat"),
    q = c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95),
    qs = NULL)

\method{hist}{SpatialVx}(x, ..., time.point = 1, model = 1)

\method{plot}{SpatialVx}(x, ..., set.pw = FALSE, time.point = 1, model = 1,
                 col, zlim, horizontal = TRUE)

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
logical, is the verification set on a regular grid?  This is another feature intended for possible future funcitonality.  Most functions in this package assume the set is on a regular grid.
}
  \item{map}{
logical, should the plot function attempt to place a map onto the plot?  Only possible if the \code{loc} argument is given.
}
  \item{field.type, units}{
character used for plot labelling and printing information to the screen.  Describes what variable and in what units the field represents.
}
  \item{data.name}{
character vector describing the verification set.  May have length equal to the number of forecasts and verification fields or it may additionally have as the first entry an overall name for the set.  For example, for the geometric cases from the ICP, one might have a vector that begins with \dQuote{ICP geometric cases} followed by \dQuote{geom000}, \dQuote{geom001}, \dQuote{geom002}, etc.  Subsequent plot labels, etc., may (or may not) make use of this information.
}
  \item{q}{
numeric vector giving the values of quantiles to be used for thresholds.  Only used if \code{thresholds} is NULL.
}
  \item{qs}{
character vector describing the quantiles used.  Again, only used if \code{thresholds} is NULL.  This is for subsequent plot/print labelling purposes.
}
   \item{set.pw}{logical, should a panel of plots be determined and set by the function.}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
  \item{model}{numeric indicating which forecast model to select for the analysis.}
   \item{col, zlim, horizontal}{optional arguments to \code{image}, and/or \pkg{fields} functions \code{poly.image} and \code{image.plot}}
    \item{x, object}{list object of class \dQuote{SpatialVx}.}
    \item{loc.byrow}{logical determining whether to set up the location matrices using \code{byrow} = TRUE or FALSE (for use with \code{poly.image}.}
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
   \item{data.name}{character vector giving (same as input argument).}
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
data(UKobs6)
data(UKfcst6)
data(UKloc)

hold <- make.SpatialVx(UKobs6, UKfcst6, thresholds=c(0.01, 20.01),
    loc=UKloc, field.type="Precipitation", units="mm/h",
    data.name=c("Nimrod", "Observations 6", "Forecast 6"), map=TRUE)

hold

plot(hold)

\dontrun{
# Stage 2 Observation files from
# NSSL/NCEP Spring 2005 Forecast Experiment.
data(obs0426)
data(obs0513)
data(obs0514)
data(obs0518)
data(obs0519)
data(obs0525)
data(obs0601)
data(obs0603)
data(obs0604)

st2 <- array(c(c(obs0426), c(obs0513), c(obs0514), c(obs0518),
    c(obs0519), c(obs0525), c(obs0601), c(obs0603), c(obs0604)),
    dim=c(601, 501, 9))

rm(obs0426, obs0513, obs0514, obs0518, obs0519, obs0525, obs0601,
    obs0603, obs0604)

# wrf2caps
data(wrf2caps0425)
data(wrf2caps0512)
data(wrf2caps0513)
data(wrf2caps0517)
data(wrf2caps0518)
data(wrf2caps0524)
data(wrf2caps0531)
data(wrf2caps0602)
data(wrf2caps0603)

# wrf4ncar
data(wrf4ncar0425)
data(wrf4ncar0512)
data(wrf4ncar0513)
data(wrf4ncar0517)
data(wrf4ncar0518)
data(wrf4ncar0524)
data(wrf4ncar0531)
data(wrf4ncar0602)
data(wrf4ncar0603)

# wrf4ncep
data(wrf4ncep0425)
data(wrf4ncep0512)
data(wrf4ncep0513)
data(wrf4ncep0517)
data(wrf4ncep0518)
data(wrf4ncep0524)
data(wrf4ncep0531)
data(wrf4ncep0602)
data(wrf4ncep0603)

wrf2caps <- array(c(c(wrf2caps0425), c(wrf2caps0512),
    c(wrf2caps0513), c(wrf2caps0517), c(wrf2caps0518),
    c(wrf2caps0524), c(wrf2caps0531), c(wrf2caps0602),
    c(wrf2caps0603)),
    dim=c(601, 501, 9))

wrf4ncar <- array(c(c(wrf4ncar0425), c(wrf4ncar0512),
    c(wrf4ncar0513), c(wrf4ncar0517), c(wrf4ncar0518),
    c(wrf4ncar0524), c(wrf4ncar0531), c(wrf4ncar0602),
    c(wrf4ncar0603)),
    dim=c(601, 501, 9))

wrf4ncep <- array(c(c(wrf4ncep0425), c(wrf4ncep0512),
    c(wrf4ncep0513), c(wrf4ncep0517), c(wrf4ncep0518),
    c(wrf4ncep0524), c(wrf4ncep0531), c(wrf4ncep0602),
    c(wrf4ncep0603)),
    dim=c(601, 501, 9))

rm(wrf2caps0425, wrf2caps0512, wrf2caps0513, wrf2caps0517,
    wrf2caps0518, wrf2caps0524, wrf2caps0531, wrf2caps0602,
    wrf2caps0603, wrf4ncar0425, wrf4ncar0512, wrf4ncar0513,
    wrf4ncar0517, wrf4ncar0518, wrf4ncar0524, wrf4ncar0531,
    wrf4ncar0602, wrf4ncar0603, wrf4ncep0425, wrf4ncep0512,
    wrf4ncep0513, wrf4ncep0517, wrf4ncep0518, wrf4ncep0524,
    wrf4ncep0531, wrf4ncep0602, wrf4ncep0603)

fcst <- list(wrf2caps, wrf4ncar, wrf4ncep)

rm(wrf2caps, wrf4ncar, wrf4ncep)

# Now, create the object.
data(ICPg240Locs)

ICPreal <- make.SpatialVx(st2, fcst, thresholds=c(0.1, 20.1),
    loc=ICPg240Locs, projection=TRUE, loc.byrow = TRUE,
    time.vals=c(2005042600, 2005051300, 2005051400, 2005051800,
        2005051900, 2005052500, 2005060100, 2005060300, 2005060400),
    map=TRUE, field.type="Precipitation", units="mm/h",
    data.name=c("Stage II Analysis", "WRF 2 CAPS", "WRF 4 NCAR",
        "WRF 4 NCEP"))

pdf("ICPreal.pdf")
plot(ICPreal, set.pw=FALSE)
dev.off()

}
}
\keyword{ manip }
\keyword{ hplot }
