\name{warper}
\alias{warper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Warp
}
\description{
Estimate an image warp
}
\usage{
warper(Im0, Im1, p0, init, s, imethod = "bicubic", lossfun = "Q", 
    lossfun.args = list(beta = 0, Cmat = NULL), grlossfun = "defaultQ", 
    lower, upper, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Im0, Im1}{
Numeric matrices giving the zero- and one-energy images.  The \code{Im1} image is ultimately warped into the \code{Im0} image.
}

  \item{p0}{
nc by 2 matrix giving the zero-energy control points.
}
  \item{init}{
nc by 2 matrix giving an initial estimate of the one-energy control points.
}
  \item{s}{
Two-column matrix giving the full set of locations.  Works best if these are integer-valued coordinate indices.
}
  \item{imethod}{
character giving he interpolation method to use.  May be one of "round", "bilinear" or "bicubic".
}
  \item{lossfun}{
Function giving the loss function over which to optimize the warp.  Default is \code{Q}, see \code{args{Q}} to see the required arguments for this function.
}
  \item{lossfun.args}{
A list giving optional arguments to \code{lossfun}.
}
  \item{grlossfun}{
(optional) function giving the gradient of the loss function given by \code{lossfun}.
}
  \item{lower, upper}{
(optional) arguments to the \code{nlminb} function which is used to optimize the loss function.
}
 
  \item{verbose}{
logical, should progress information be printed to the screen?
}
  \item{\dots}{
Optional arguments to \code{nlminb}.
}
}
\details{
A pair-of-thin-plate-splines image warp is estimated by optimizing a loss function using nlminb.  It can be very difficult to get a good estimate.  It is suggested, therefore, to obtain good initial estimates for the one-energy control points.  The function \code{iwarper} can be useful in this context.
}
\value{ A list object of class \dQuote{warped} is returned with components:
    \item{Im0, Im1, Im1.def}{Matrices giving the zero- and one-energy images and the deformed one-energy image, resp.}
    \item{p0, p1}{zero- and one-energy control points, resp.}
    \item{sigma}{Estimated standard error of the mean difference between the zero-energy and deformed one-energy images.}            "warped.locations" "init"            
    \item{s, imethod, lossfun, lossfun.args}{Same as input arguments.}
    \item{theta}{The matrices defining the image warp, L, iL and B, where the last is the bending energy, and the first two are nc + 3 by nc + 3 matrices describing the control points and inverse control-point matrices.}
    \item{arguments}{Any arguments passed via \ldots}
    \item{fit}{The output from nlminb.}
    \item{proc.time}{The process time.}
}
\references{
Dryden, I. L. and K. V. Mardia (1998) \emph{Statistical Shape Analysis}.  Wiley, New York, NY, 347pp.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{iwarper}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ spatial }
% \keyword{ model }
