\name{iwarper}
\alias{iwarper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Warping By Hand
}
\description{
Instigate an image warp by selecting control points in the zero- and one-energy images by hand.
}
\usage{
iwarper(x0, x1, nc = 4, labcol = "magenta", 
    col = c("gray", tim.colors(64)), zlim, cex = 2, alwd = 1.25, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0, x1}{
Numeric matrices giving the zero- and one-energy images.  The \code{x1} image is ultimately warped into the \code{x0} image.
}
  \item{nc}{
integer giving the number of control points to select.
}
  \item{labcol}{
character describing the color to use when labeling the control points as they are selected.
}
  \item{col}{
The color scheme to use in plotting the images.
}
  \item{zlim}{
Range of values for the color scheme in plotting the images.
}
  \item{cex}{
The usual \code{cex} parameter.  See the help file for \code{par} for more information.
}
  \item{alwd}{
line width for the arrows added to the deformation plot.
}
  \item{\dots}{
Optional arguments to \code{poly.image}.
}
}
\details{
A pair-of-thin-plate-splines image warp mapping is estimated by hand.  See Dryden and Mardia (1998) Chapter 10 for more information.
}
\value{
A list object of class \dQuote{iwarped}

    \item{Im0, Im1, Im1.def}{The zero- and one-energy images and the deformed one-energy image, resp.}
    \item{p0, p1}{The nc by 2 column matrices of hand-selected zero- and one-energy control points, resp.} 
    \item{warped.locations, s}{Two-column matrices giving the entire set of warped locations and original locations, resp.}
    \item{theta}{The matrices defining the image warp, L, iL and B, where the last is the bending energy, and the first two are nc + 3 by nc + 3 matrices describing the control points and inverse control-point matrices.}

}
\references{
Dryden, I. L. and K. V. Mardia (1998) \emph{Statistical Shape Analysis}.  Wiley, New York, NY, 347pp.
}
\author{
Eric Gilleland
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{warper}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ spatial }
% \keyword{ models }

