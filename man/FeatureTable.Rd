\name{FeatureTable}
\alias{FeatureTable}
\alias{ci.FeatureTable}
\alias{print.FeatureTable}
\alias{summary.FeatureTable}
\title{
Feature-Based Contingency Table
}
\description{
Create a feature-based contingency table from a matched object and calculate some summary scores with their standard errors.
}
\usage{
FeatureTable(x, fudge = 1e-08, hits.random = NULL, correct.negatives = NULL, fA = 0.05)

\method{ci}{FeatureTable}(x, alpha = 0.05, ...)

\method{print}{FeatureTable}(x, ...)

\method{summary}{FeatureTable}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{FeatureTable}: An object of class \dQuote{matched} (e.g., from \code{deltamm} or \code{centmatch}.

\code{print}: An object of class \dQuote{FeatureTable}.
}
   \item{object}{An object of class \dQuote{FeatureTable}.}
  \item{fudge}{
    value added to denominators of scores to ensure no division by zero.  Set to zero if this practice is not desired.
}
  \item{hits.random}{
    If a different value for random hits in the caluclation of GSS than provided is desired, it can be given here.  Default uses Eq (3) from Davis et al (2009).
}
  \item{correct.negatives}{

    If a different value for correct negatives than provided is desired, it can be given here.  Default uses Eq (4) from Davis et al (2009).

}
  \item{fA}{
numeric between zero and 1 giving the fraction of area occupied for the purpose of matching as in Davis et al (2009).
}
  \item{alpha}{numeric between zero and one giving the (1 - alpha) * 100 percent confidence level.}
  \item{\dots}{Additional arguments to \code{ci}.}
}
\details{

This function takes an object of class \dQuote{matched} and calculates a contingency table based on matched and unmatched objects.  If no value for correct negatives is given, then it will also determine them based on Eq (3) from Davis et al (2009).  The following contingency table scores and their standard errors (based on their usual traditional version) are returned.  It should be noted that the standard errors may not be entirely meaningful because they do not capture the uncertainty associated with identiying, merging and matching features within the fields.  Neverhteless, they are calculated here for investigative purposes.  Note that hits are determined by number of matched objects, which for some matching algorithms can mean that features are matched more than once (e.g., if using \code{centmatch}).  In essence, this fact may artificially increase thenumber of hits.  On the other hand, situations exist where such handling may be more appropriate than not having duplicate matches.

hits are determined by the total number of matched features.

false alarms are the total number of unmatched forecast features.

misses are the total number of unmatched observed features.

correct negatives are less obviously defined.  If the user does not supply a value, then these are calculated according to Eq (4) in Davis et al (2009).

GSS: Gilbert skill score (aka Equitable Threat Score) based on Eq (2) of Davis et al (2009).

POD: probability of detecting an event (aka the hit rate).

false alarm rate: (aka probability of false detection) is the ratio of false alarms to the number of false alarms and correct negtives.

FAR: the false alarm ratio is the ratio of false alarms to the total forecast events (in this case, the total number of forecast features in the field).

HSS: Heidke skill score

The \code{print} method function simply calls \code{summary}, which prints the feature-based contingency table in addition to calling \code{ci}.  The confidence intervals are based on the normal approximation method using the estimated standard errors, which themselves are suspicious.  In any case, the intervals can give a feel for some of the uncertainty associated with the scores, but should not be considered as solid.

}
\value{
A list with inherited attributes from x and components:
  \item{estimates }{named numeric vector giving the estimated scores.}
  \item{se }{named numeric vector giving the estimated standard errors of the scores.}
  \item{feature.contingency.table}{named numeric vector giving the feature-based contingency table.}
}
\references{
Davis, C. A., Brown, B. G., Bullock, R. G. and Halley Gotway, J. (2009) The Method for Object-based Diagnostic Evaluation (MODE) applied to numerical forecasts from the 2005 NSSL/SPC Spring Program.  \emph{Wea. Forecsting}, \bold{24}, 1252--1267, DOI: 10.1175/2009WAF2222241.1.
}
\author{
Eric Gilleland
}
\note{
Standard error estimates are based on the univariate equivalent formulations, which do not account for uncertainties introduced in the feature identification, merging/clustering and matching.  They should not be considered as legitimate, and resulting confidence intervals should be mistrusted. 
}

\seealso{
To identify features in the fields: \code{\link{FeatureFinder}}

To match (and merge) features: \code{\link{centmatch}}, \code{\link{deltamm}}
}
\examples{
##
## See help file for 'deltamm' for examples.
##
}
\keyword{ htest }
\keyword{ manip }
