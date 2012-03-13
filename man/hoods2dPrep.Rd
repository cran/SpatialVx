\name{hoods2dPrep}
\alias{hoods2dPrep}
\title{Set up a list object that can be used by the hoods2d function.}
\description{Creates a list object with all of the information needed to run the hoods2d function, as well as other neighborhood approach functions.}
\usage{
hoods2dPrep(Fcst.name, Vx.name, thresholds = NULL, Pe = NULL, levels = NULL, max.n = NULL,  subset = NULL, loc = NULL, qs = NULL, units = NULL, smooth.fun = "hoods2dsmooth", smooth.params = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fcst.name}{character giving the name of a forecast matrix object.}
  \item{Vx.name}{character giving the name of an observed matrix object.  The matrix should have same dimension as that for Fcst.name.}
  \item{thresholds}{(optional) numeric vector of length q >= 1, or or q X 2 matrix giving the thresholds to be applied.  If a matrix, the first column are thresholds to be applied to the forecast field, and the second column for the observed field.  If NULL, then thresholds are taken to be the quantiles: 0, .10, .25, .33, .50, .66, .75, .90, .95 of each field, resp.}
  \item{Pe}{(optional) numeric vector of length q >= 1 to be applied to the fields sPy and possibly sPx (see help file for \code{hoods2d}).  If NULL, then it is taken to be the most relaxed requirement (i.e., that an event occurs at least once in a neighborhood) of 'Pe'=1/(nlen^2), where nlen is the length of the neighborhood.}
  \item{levels}{numeric vector giving the successive values of the smoothing parameter.  For example, for the default method, these are the neighborhood lengths over which the levels^2 nearest neighbors are averaged for each point.  Values should make sense for the specific smoothing function.  For example, for the default method, these should be odd integers.}
  \item{max.n}{(optional) single numeric giving the maximum neighborhood length to use.  Only used if levels are not passed.}
  \item{subset}{(optional) numeric vector indicating a subset of points over which to calculate the final statistics.  If given, then fields will be smoothed over the entire forecast/observed field, but resulting statistics will only be computed on the subset.  The subset numbers should be within the range 1 to n*m, where the forecast (observed) field has dimension n X m.  The subset is applied to the concatenation of the fields, which concatenates the columns together.}
  \item{loc}{(optional) numeric n*m X 2 matrix giving lon/lat coordinates for the verification set.  Currently, this is not used for anything, but may have a use in the future, e.g., for adding a map to plots of neighborhood smoothed fields.}
  \item{qs}{(optional) character vector giving the values of the thresholds.  This is primarily used, e.g., if the thresholds differ for the forecast and observed fields (e.g., if quantiles are used instead of hard values).  This information is used by \code{hoods2dPlot}.}
  \item{units}{(optional) character giving the name of the units of the verification set.}
  \item{smooth.fun}{character giving the name of a smoothing function to be applied.  Default is an average over the n^2 nearest neighbors, where n is taken to be each value of the \code{levels} argument.}
  \item{smooth.params}{additional arguments to \code{kernel2dsmooth} (\code{kenrel2dmeitsjer}) for the specific smoothing function.  The default method does not take any extra parameters.}
}
\details{This function is simply used to simplify later analyses and plotting.  It gathers all of the information needed to run \code{hoods2d}, \code{plot.hoods2d} and \code{hoods2dPlot}.}
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
   \item{units}{Same as in Arguments above.}
}
\author{Eric Gilleland}

\seealso{
\code{\link{hoods2d}}
}
\examples{
data(geom001)
data(geom000)
data(ICPg240Locs)
hold <- hoods2dPrep( "geom001", "geom000", thresholds=c(0.01,50.01), levels=c(1, 3, 5, 9, 17, 33, 65, 129, 257), loc=ICPg240Locs, units="in/100")
}
\keyword{ list }
