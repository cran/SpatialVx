\name{interester}
\alias{interester}
\alias{print.interester}
\alias{summary.interester}
\title{
Feature Interest
}
\description{
Calculate interest maps for specific feature comparisons and compute the total interest, as well as median of maximum interest.
}
\usage{
interester(x, properties = c("cent.dist", "angle.diff", "area.ratio", "int.area",
    "bdelta", "haus", "ph", "med", "msd", "fom", "minsep"),
    weights = c(0.24, 0.12, 0.17, 0.12, 0, 0, 0, 0, 0, 0, 0.35),
    b1 = c(35, 30, 0, 0, 0.5, 35, 20, 40, 120, 1, 40),
    b2 = c(100, 90, 0.8, 0.25, 85, 400, 200, 200, 400, 0.25, 200),
    verbose = FALSE, ...)

\method{print}{interester}(x, ...)

\method{summary}{interester}(object, ...,
    min.interest = 0.8, long = TRUE, silent = FALSE)
}
\arguments{
  \item{x}{
\code{interester}: object of class \dQuote{features} or \dQuote{matched}.

\code{print}: object of class \dQuote{interester}.
}
  \item{object}{object of class \dQuote{interester}.}
  \item{properties}{
character vector naming which properties from \code{FeatureComps} should be considered in the total interest.
}
  \item{weights}{
numeric of length equal to the length of \code{properties}.  Weights equal to zero will result in the removal of those properties from \code{properties} so that the default computes only those values utilized in Davis et al (2009).
}
  \item{b1,b2}{
All interest maps (except that for \dQuote{fom}) are piecewise linear, and of the form: f(x) = 1,0 (depending on the property) if x < \code{b1}, x <= \code{b1}, x > \code{b2} or x >= \code{b2} (see details for more information).
}
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{min.interest}{numeric between zero and one giving the desired minimum value of interest.  Only used for display purposes.  If \code{long} is TRUE, a dashed line is printed where the total interest falls below this value.  If \code{long} is FALSE, only the total interest values above this value are displayed.  All values are returned invisibly regardless of the values of \code{min.interest} and \code{long}.}
  \item{long}{logical, should all interest values be displayed (TRUE) or only those above \code{min.interest} (FALSE)?}
  \item{silent}{logical, should summary information be displayed to the screen (FALSE)?}
  \item{\dots}{
\code{interester}: optional arguments to \code{FeatueComps}.

Not used by \code{print} or \code{summary}.
}
}
\details{
This function calculates the feature interest according to the MODE algorithm described in Davis et al (2009).  Properties that can be computed are those available in \code{FeatureComps}, except for \dQuote{bearing}.  Interest maps are computed according to piece-wise linear functions (except for \dQuote{fom}) depending on the property.  For all properties besides \dQuote{area.ratio}, \dQuote{int.area} and \dQuote{fom}, the interest maps are of the form:

f(x) = 1, if x <= b1

f(x) = a0 + a1 * x, if x > b1 and x <= b2, where a1 = -1/(b2 - b1) and a0 = 1 - a1 * b1

f(x) = 0, if x > b2

For properties \dQuote{area.ratio} and \dQuote{int.area}, the interest maps are of the form:

f(x) = 0, if x < b1

f(x) = a0 + a1 * x, if x >= b1 and x < b2, where a1 = 1/(b2 - b1) and a0 = 1 - a1 * b2

f(x) = 1, if x >= b2

Finally, for \dQuote{fom}, a function that tries to give as much weight to values near one is applied.  It is given by:

f(x) = b1 * exp(-0.5 * ((x - 1) / b2)^4)

The default values for b1 and b2 will not necessarily give the same results as in Davis et al (2009), but also, the distance map for their intersection area ratio differs from that here.  The interest function for FOM is further restricted to fall within the interval [0, 1], so care should be taken if b1 and/or b2 are changed for this function.

The \code{interester} function calculates the individual interest values for each property and each pair of features, and returns both these individual interest values as well as a matrix of total interest.  The \code{print} function will print the entire matrix of individual interest values if there are fewer than twenty pairs of features, and will print their summary otherwise.  The \code{summary} function will order the total interest from highest to lowest and print this information (along with which feature pairs correspond to the total interest value).  It will also calculate the median of maximum interest (MMI) as suggested by Davis et al (2009).  If there is only one feature in either field, then this value will just be the maximum total interest.

The centroid distance property is less meaningful if the sizes of the two features differ greatly, and therefore, the interest value for this property is further multiplied by the area ratio of the two features.  Similarly, angle difference is less meaningful if one or both of the features are circular in shape.  Therefore, this property is multiplied by the following factor, following Davis et al (2009) Eq (A1), where r1 and r2 are the aspect ratios (defined as the length of the minor axis divided my that of the major axis) of the two features, respectively.

    sqrt( [ (r1 - 1)^2 / (r1^2 + 1) ]^0.3 * [ (r2 - 1)^2 / (r2^2 + 1) ]^0.3 )

The \code{print} function displays either the individualinterest values for each property and feature pairings, or more often, a summary of these values (if the display would otherwise be too large).  It also shows a matrix whose rows are observed features and columns forecast features, with the total interest values therein associated.

\code{summary} shows the sorted total interest from highest to lowest for each pair.  A dashed line separates the values above \code{min.interest} from those below, and if \code{long} is TRUE, then values below that line are not displayed.  It also reports the median of maximum interest (MMI) defined by Davis et al (2009) as an overall feature-based summary of forecast performance.  It is derived by collecting the row maxima and column maxima from the total interest matrix, shown by the \code{print} command, into a vector, and then finding the median of this vector.

}

\value{

A list of class \dQuote{interester} is returned with components:
  \item{interest }{matrix whose named rows correspond to the each property that was calculated and whose columns are feature pairings.  The values are the interest calculated for the specific property and pair of features.}
  \item{total.interest }{matrix of total interest for each pair of features where rows are observed features and columns are forecast features.}

If no features are available in either field, NULL is returned.

Nothing is returned by print.

summary invisibly returns a list object of class \dQuote{summary.interester} with components:

  \item{sorted.interest}{similar to the interest component from the value returned by interester, but sorted from highest to lowest interest, along with the feature number information for each field.}
  \item{mmi}{the median of maximum interest value.}

}

\references{
Davis, C. A., Brown, B. G., Bullock, R. G. and Halley Gotway, J. (2009) The Method for Object-based Diagnostic Evaluation (MODE) applied to numerical forecasts from the 2005 NSSL/SPC Spring Program.  \emph{Wea. Forecsting}, \bold{24}, 1252--1267, DOI: 10.1175/2009WAF2222241.1.
}
\author{
Eric Gilleland
}
\note{
The terminology used for features within the entire \pkg{SpatialVx} package attempts to avoid conflict with terminology used by R.  So, for example, the term property is used in lieu of \dQuote{attributes} so as not to be confused with R object attributes.  The term \dQuote{feature} is used in place of \dQuote{object} to avoid confusion with an R object, etc.
}

\seealso{
Identifying features: \code{\link{FeatureFinder}}

Functions for calculating the properties: \code{\link{FeatureComps}}, \code{\link{FeatureProps}}
}
\examples{

x <- y <- matrix(0, 100, 100)
x[ 2:3, c(3:6, 8:10) ] <- 1
y[ c(4:7, 9:10), c(7:9, 11:12) ] <- 1
     
x[ 30:50, 45:65 ] <- 1
y[ c(22:24, 99:100), c(50:52, 99:100) ] <- 1
     
hold <- make.SpatialVx(x, y, field.type = "contrived", units = "none",
    data.name = "Example", obs.name = "x", model.name = "y" )
     
look <- FeatureFinder(hold, smoothpar = 0.5) 
     
look2 <- interester(look)
look2
summary(look2)

}
\keyword{ manip }
