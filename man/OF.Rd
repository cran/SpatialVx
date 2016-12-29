\name{OF}
\alias{OF}
\alias{OF.default}
\alias{OF.SpatialVx}
\alias{plot.OF}
\alias{print.OF}
\alias{hist.OF}
\alias{summary.OF}
\title{
Optical Flow Verification
}
\description{
Perform verification using optical flow as described in Marzban and Sandgathe (2010).
}
\usage{
OF(x, ...)

\method{OF}{default}(x, ..., xhat, W = 5, grads.diff = 1, center = TRUE, 
    cutoffpar = 4, verbose = FALSE)

\method{OF}{SpatialVx}(x, ..., time.point = 1, obs = 1, model = 1, W = 5, grads.diff = 1,
    center = TRUE, cutoffpar = 4, verbose = FALSE)

\method{plot}{OF}(x, ...)

\method{print}{OF}(x, ...)

\method{hist}{OF}(x, ...)

\method{summary}{OF}(object, ...)
}
\arguments{
  \item{x, xhat}{
Default: m by n matrices describing the verification and forecast fields, resp.  The forecast field is considered the initial field that is morphed into the final (verification) field.

\code{OF.SpatialVx}: list object of class \dQuote{SpatialVx}.

\code{plot}, \code{hist} and \code{print} methods: list object as returned by \code{OF}.
}
  \item{object}{
list object as returned by \code{OF}.
}
  \item{W}{
numeric/integer giving the window size (should be no smaller than 5).
}
  \item{grads.diff}{
1 or 2 describing whether to use first or second differences in finding the first derivative.
}
  \item{center}{
logical, should the fields be centered before performing the optical flow?
}
  \item{cutoffpar}{
numeric, set to NaN everything exceeding median +/- \code{cutoffpar}*sd.
}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
    \item{time.point}{numeric or character indicating which time point from the \dQuote{SpatialVx} verification set to select for analysis.}
    \item{obs, model}{numeric indicating which observation/forecast model to select for the analysis.}
  \item{\dots}{
For \code{OF}: optional arguments to the \code{optim} function (cannot be \code{par}, \code{fn}, \code{gr} or \code{method}).  See details section for \code{plot} and \code{hist} method functions.  Not used by the \code{summary} method function.
}
}
\details{
Estimates the optical flow of the forecast field into the verification field.  Letting I_o(x,y) and I_f(x,y) represent the intensities of each field at coordinate (x,y), the collection of pairs (dx, dy) is the optical flow field, where:

I_o(x,y) ~ I_f(x,y) + [partial(I_f) wrt x]*dx + [partial(I_f) wrt y]*dy.

The procedure follows that proposed by Lucas and Kanade (1981) whereby for some window, W, it is assumed that all dx (dy) are assumed constant, and least squares estimation is used to estimate dx and dy (see Marzban and Sandgathe, 2010 for more on this implementation).  This function iteratively calls optflow for each window in the field.

The above formulation is linear in the parameters.  Marzban and Sandgathe (2010) also introduce an additive error component, which leads to a nonlinear version of the above.  Namely,

I_o(x,y) ~ I_f(x,y) + [partial(I_f) wrt x]*dx + [partial(I_f) wrt y]*dy + A(x,y).

See Marzban and Sandgathe for more details.

The plot method function can produce a figure like that of Fig. 1, 5, and 6 in Marzban and Sandgathe (2010) or with option \code{full=TRUE}, even more plots.  Optional arguments that may be passed in via the ellipses include: \code{full} (logical, produce a figure analogous to Fig. 1, 5 and 6 from Marzban and Sandgathe (2010) (FALSE/default) or make more plots (TRUE)), \code{scale} (default is 1 or no scaling, any numeric value by which the fields are divided/scaled before plotting), \code{of.scale} (default is 1, factor by which display vectors can be magnified), \code{of.step} (plot OF vectors every of.step, default is 4), \code{prop} (default is 2, value for \code{prop} argument in the call to \code{rose.diag} from package \pkg{CircStats}), \code{nbins} (default is 40, number of bins to use in the call to \code{rose.diag}).

The \code{hist} method function produces a two-dimensional histogram like that of Fig. 3 and 7 in Marzban and Sandgathe (2010).  It can also take various arguments passed via the ellipses.  They include: \code{xmin}, \code{xmax}, \code{ymin}, \code{ymax} (lower and upper bounds for the histogram breaks in the x- (angle) and y- (magnitude/displacement error) directions, resp.  Defaults to (0,360) and (0,4)), nbreaks (default is 100, the number of breaks to use).

The \code{summary} method mostly uses the \code{stats} function from package \pkg{fields} to summarize results of the errors, but also uses \code{circ.summary} from package \pkg{CircStats} for the angular errors.
}
\value{
OF returns a list object of class \dQuote{OF} with components:
  \item{data}{list with components x and xhat containing the data.}
  \item{data.name}{character vector giving the names of the verification and forecast fields.}
  \item{call}{object of class \dQuote{call} giving the original function call.}
  \item{rows,cols}{numeric vector giving the rows and columns used for finding the centers of windows.  Needed by the plot and hist method functions.}
  \item{err.add.lin}{m by n matrix giving the linear additive errors (intensities).}
  \item{err.mag.lin}{m by n matrix giving the linear magnitude (displacement) errors.}
  \item{err.ang.lin}{m by n matrix giving the linear angular errors.}
  \item{err.add.nlin,err.mag.nlin,err.ang.nlin}{same as above but for nonlinear errors.}
  \item{err.vc.lin,err.vr.lin,err.vc.nlin,err.vr.nlin}{m by n matrices giving the x- and y- direction movements for the linear and nonlinear cases, resp.}
The hist method function invisibly returns a list object of class \dQuote{OF} that contains the same object that was passed in along with new components:
 \item{breaks}{a list with components x and y giving the breaks in each direction}
 \item{hist.vals}{itself a list with components xb, yb (the number of breaks -1 used for each direction), and nb (the histogram values for each break)}
The plot and summary mehtod functions do not return anything.
}
\references{
Lucas, B D. and Kanade, T. (1981)  An iterative image registration technique with an application to stereo vision.  \emph{Proc. Imaging Understanding Workshop}, DARPA, 121--130.

Marzban, C. and Sandgathe, S. (2010)  Optical flow for verification.  \emph{Wea. Forecasting}, \bold{25}, 1479--1494, doi:10.1175/2010WAF2222351.1.
}
\author{
Caren Marzban, marzban \dQuote{at} u.washington.edu, with modifications by Eric Gilleland
}

\seealso{
\code{\link{optflow}}, \code{\link{optim}}, \code{\link{circ.summary}}
}
\examples{

\dontrun{
data(hump)
initial <- hump$initial
final <- hump$final
look <- OF(final, xhat=initial, W=9, verbose=TRUE)
plot(look) # Compare with Fig. 1 in Marzban and Sandgathe (2010).
par(mfrow=c(1,1))
hist(look) # 2-d histogram.
plot(look, full=TRUE) # More plots.
summary(look)

# Another way to skin the cat.
hold <- make.SpatialVx( final, initial, field.type = "Bi-variate Gaussian",
    obs.name = "final", model.name = "initial" )

look2 <- OF(hold, W=9, verbose=TRUE)
plot(look2)
par(mfrow=c(1,1))
hist(look2)
plot(look2, full=TRUE)
summary(look2)
}
}
\keyword{ math }
\keyword{ regression }
