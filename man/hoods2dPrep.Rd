\name{hoods2dPrep}
\alias{hoods2dPrep}
\alias{plot.hoods2dPrep}
\alias{hist.hoods2dPrep}
\title{Set up a list object that can be used by the hoods2d function.}
\description{Creates a list object with all of the information needed to run the hoods2d function, as well as other neighborhood approach functions.}
\usage{
hoods2dPrep(Fcst.name, Vx.name, thresholds = NULL, Pe = NULL, levels = NULL, max.n = NULL,
		subset = NULL, loc = NULL, qs = NULL, field.type = "", units = NULL,
		smooth.fun = "hoods2dsmooth", smooth.params = NULL, mapit = FALSE)
\method{plot}{hoods2dPrep}(x, ...)
\method{hist}{hoods2dPrep}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fcst.name}{character giving the name of a forecast matrix object.}
  \item{Vx.name}{character giving the name of an observed matrix object.  The matrix should have same dimension as that for Fcst.name.}
  \item{x}{list object as returned by \code{hoods2dPrep}.}
  \item{thresholds}{(optional) numeric vector of length q >= 1, or or q X 2 matrix giving the thresholds to be applied.  If a matrix, the first column are thresholds to be applied to the forecast field, and the second column for the observed field.  If NULL, then thresholds are taken to be the quantiles: 0, .10, .25, .33, .50, .66, .75, .90, .95 of each field, resp.}
  \item{Pe}{(optional) numeric vector of length q >= 1 to be applied to the fields sPy and possibly sPx (see help file for \code{hoods2d}).  If NULL, then it is taken to be the most relaxed requirement (i.e., that an event occurs at least once in a neighborhood) of 'Pe'=1/(nlen^2), where nlen is the length of the neighborhood.}
  \item{levels}{numeric vector giving the successive values of the smoothing parameter.  For example, for the default method, these are the neighborhood lengths over which the levels^2 nearest neighbors are averaged for each point.  Values should make sense for the specific smoothing function.  For example, for the default method, these should be odd integers.}
  \item{max.n}{(optional) single numeric giving the maximum neighborhood length to use.  Only used if levels are not passed.}
  \item{subset}{(optional) numeric vector indicating a subset of points over which to calculate the final statistics.  If given, then fields will be smoothed over the entire forecast/observed field, but resulting statistics will only be computed on the subset.  The subset numbers should be within the range 1 to n*m, where the forecast (observed) field has dimension n X m.  The subset is applied to the concatenation of the fields, which concatenates the columns together.}
  \item{loc}{(optional) numeric n*m X 2 matrix giving lon/lat coordinates for the verification set.  Currently, this is not used for anything, but may have a use in the future, e.g., for adding a map to plots of neighborhood smoothed fields.}
  \item{qs}{(optional) character vector giving the values of the thresholds.  This is primarily used, e.g., if the thresholds differ for the forecast and observed fields (e.g., if quantiles are used instead of hard values).  This information is used by \code{hoods2dPlot}.}
  \item{field.type}{character naming the type of field given by the verification set (e.g., "precipitation", "wind speed", etc.).}
  \item{units}{(optional) character giving the name of the units of the verification set.}
  \item{smooth.fun}{character giving the name of a smoothing function to be applied.  Default is an average over the n^2 nearest neighbors, where n is taken to be each value of the \code{levels} argument.}
  \item{smooth.params}{additional arguments to \code{kernel2dsmooth} (\code{kenrel2dmeitsjer}) for the specific smoothing function.  The default method does not take any extra parameters.}
  \item{mapit}{logical stating whether or not the \code{loc} argument corresponds to lon/lat coordinates (TRUE) or not (FALSE).  Used by plot method function to add a map to the plots if TRUE.  Note that maps placed on the plot will generally not have the correct projections.}
  \item{\dots}{optional arguments to the \code{hist} function.  Not used by the plot method function.}
}
\details{This function is simply used to simplify later analyses and plotting.  It gathers all of the information needed to run \code{hoods2d}, \code{plot.hoods2d} and \code{hoods2dPlot}.

The plot method function makes image plots of the verification set with gray contours showing the chosen thresholds for each field.  It will try to add a map if the map component is TRUE and the loc coordinates are lon/lat and in the correct format (see the map function from package maps).

The hist method function creates several histograms in two columns.  The first column correspond to the verification field and the second to the forecast field.  The first row shows the histogram for all values in the fields, and each successive row shows histograms for the field only where values are >= each threshold from the thresholds component of the hoods2dPrep object.  That is, (N+1) X 2 panels of histograms are displayed where N is the number of thresholds.
}
\value{A list object with components:
   \item{Fcst.name}{character giving the name of the forecast matrix of the verification set.}
   \item{Vx.name}{character giving the name of the observed matrix of the verification set.}
   \item{xdim}{numeric vector of length 2 giving the dimensions of the verification set.}
   \item{Nxy}{numeric giving the total number of grid points in the verification set (i.e., \code{prod(xdim)}).}
   \item{levels}{numeric vector giving the desired neighborhood lengths.}
   \item{max.n}{numeric giving the largest neighborhood length to be applied.}
   \item{thresholds}{q X 2 matrix whose columns give the thresholds for the forecast and observed fields, reps.}
   \item{qs}{Same as in Arguments above.}
   \item{Pe}{numeric vector of length q >= 1 (i.e., same as number of rows of threshold) giving the thresholds for the smoothed binary fields.}
   \item{subset}{Same as in Arguments above.}
   \item{loc}{Same as in Arguments above.}
   \item{field.type}{Same as in Arguments above.}
   \item{units}{Same as in Arguments above.}
   \item{map}{Same as mapit in Arguments above.}
  The plot and hist method functions do not return anything.  Plots are created.
}
\author{Eric Gilleland}

\seealso{
\code{\link{hoods2d}}
}
\examples{
data(UKobs6)
data(UKfcst6)
data(UKloc)
hold <- hoods2dPrep( "UKfcst6", "UKobs6", thresholds=c(0.01,20.01), levels=c(1, 3, 5, 9, 17),
                loc=UKloc, field.type="QPF", units="mm/h", mapit=TRUE)
plot(hold)
hist(hold, col="darkblue")

\dontrun{
data(geom001)
data(geom000)
data(ICPg240Locs)
hold <- hoods2dPrep( "geom001", "geom000", thresholds=c(0.01,50.01), levels=c(1, 3, 5, 9, 17, 33, 65, 129, 257),
		loc=ICPg240Locs, field.type="QPF", units="in/100", mapit=TRUE)
plot(hold)
hist(hold, col="darkblue")
}
}
\keyword{ list }
