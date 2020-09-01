\name{SpatialVx-package}
\alias{SpatialVx-package}
\alias{SpatialVx}
\docType{package}
\title{
Spatial Forecast Verification
}
\description{
\pkg{SpatialVx} contains functions to perform many spatial forecast verification methods.
}

\details{

Primary functions include:

0. \code{make.SpatialVx}: An object that contains the verification sets and pertinent information.

1. Filter Methods:

1a. Neighborhood Methods:

Neighborhood methods generally apply a convolution kernel smoother to one or both of the fields in the verificaiton set, and then apply the traditional scores.  Most of the methods reviewed in Ebert (2008, 2009) are included in this package.  The main functions are:

\code{hoods2d}, \code{pphindcast2d}, \code{kernel2dsmooth}, and \code{plot.hoods2d}.

1b. Scale Separation Methods:

Scale separation refers to the idea of applying a band-pass filter (and/or doing a multi-resolution analysis, MRA) to the verification set.  Typically, skill is assessed on a scale-by-scale basis.  However, other techniques are also applied.  For example, denoising the field before applying traditional statistics, using the variogram, or applying a statistical test based on the variogram (these last are less similar to the spirit of the scale separation idea, but are at least somewhat related).

There is functionality to do the wavelet methods proposed in Briggs and Levine (1997).  In particular, to simply denoise the fields before applying traditional verification statistics, use 

\code{wavePurifyVx}.

To apply verification statistics to detail fields (in either the wavelet or field space), use:

 \code{waverify2d} (dyadic fields) or \code{mowaverify2d} (non-dyadic or dyadic) fields.

The intensity-scale technique introduced in Casati et al. (2004) and the new developments of the technique proposed in Casati (2010) can be performed with

\code{waveIS}.

Although not strictly a \dQuote{scale separation} method, the structure function (for which the variogram is a special case) is in the same spirit in the sense that it analyzes the field for different separation distances, and these \dQuote{scales} are separate from each other (i.e., the score does not necessarily improve or decline as the scale increases).  This package contains essentially wrapper functions to the \code{vgram.matrix} and \code{plot.vgram.matrix} functions from the fields package, but there is also a function called

\code{variogram.matrix}

that is a modification of \code{vgram.matrix} that allows for missing values.  The primary function for doing this is called

\code{griddedVgram}, which has a \code{plot} method function associated with it.

There are also slight modifications of these functions (small modifications of the fields functions) to calculate the structure function of Harris et al. (2001).  These functions are called

 \code{structurogram} (for non-gridded fields) and \code{structurogram.matrix} (for gridded fields).

The latter allows for ignoring zero-valued grid points (as detailed in Harris et al., 2001) where the former does not (they must be removed prior to calling the function).

2. Displacement Methods:  

In Gilleland et al (2009), this category was broken into two main types as field deformation and features-based.  The former lumped together binary image measures/metrics with field deformation techniques because the binary image measures inform about the \dQuote{similarity} (or dissimilarity) between the spatial extent or pattern of two fields (across the entire field).  Here, they are broken down further into those that yield only a single (or small vector of) metric(s) or measure(s) (location measures), and those that have mechanisms for moving grid-point locations to match the fields better spatially (field deformation).

2a. Distance-based and Spatial-Alignment Summary Measures: 

Gilleland (2020) introduced a new spatial alignment summary measure that falls between zero and one, with one representing a perfect match and zero a bad match.  There is one user-selectable parameter/argument that determines the rate of decrease of the measure towards zero.  Another two summary measures also incorporate intensity information.  These summaries are available via \code{Gbeta}, \code{GbetaIL} and \code{G2IL}.

In addition to the above new measures, older well-known measures are inlcuded, including: the Hausdorff metric, partial-Hausdorff measure, FQI (Venugopal et al., 2005), Baddeley's delta metric (Baddeley, 1992; Gilleland, 2011; Schwedler et al., 2011), metrV (Zhu et al., 2011), as well as the localization performance measures described in Baddeley, 1992: mean error distance, mean square error distance, and Pratt's Figure of Merit (FOM).

\code{locmeasures2d}, \code{metrV}, \code{distob}, \code{locperf}

Image moments can give useful information about location errors, and are used within feature-based methods, particularly MODE, as they give the centroid of an image (or feature), as well as the orientation angle, among other useful properties.  See the \code{imomenter} function for more details.

2b. Field deformation:

Thanks to Caren Marzban for supplying his optical flow code for this package (it has been modified some).  These functions perform the analyses described in Marzban and Sandgathe (2010) and are based on the work of Lucas and Kanade (1981).  See the help file for

\code{OF}.

Rigid transformations can be estimated using the \code{rigider} function.  To simply rigidly transform a field (or feature) using specified parameters (x- and y- translations and/or rotations), the \code{rigidTransform} function can be used.  For these functions, which may result in transformations that do not perfectly fall onto grid points, the function \code{Fint2d} can be used to interpolate from nearest grid points.  Interpolation options include \dQuote{round} (simply take the nearest location value), \dQuote{bilinear} and \dQuote{bicubic}.

2c. Features-based methods: These methods are also sometimes called object-based methods (the term \dQuote{features} is used in this package in order to differentiate from an R object), and have many similarities to techniques used in Object-Based Image Analysis (OBIA), a relatively new research area that has emerged primarily as a result of advances in earth observations sensors and GIScience (Blaschke et al., 2008).  It is attempted to identify individual features within a field, and subsequently analyze the fields on a feature-by-feature basis.  This may involve intensity error information in addition to location-specific error information.  Additionally, contingency table verifcation statistics can be found using new definitions for hits, misses and false alarms (correct negatives are more difficult to asses, but can also be done).

Currently, there is functionality for performing the analyses introduced in Davis et al. (2006,2009), including the merge/match algorithm of Gilleland et al (2008), as well as the SAL technique of Wernli et al (2008, 2009).  Some functionality for composite analysis (Nachamkin, 2004) is provided by way of placing individual features onto a relative grid so that each shares the same centroid.  Shape analysis is partially supported by way of functions to identify boundary points (Micheas et al. 2007; Lack et al. 2010).  In particular, see:

Functions to identify features: \code{FeatureFinder}

Functions to match/merge features: \code{centmatch}, \code{deltamm}, \code{minboundmatch}

Functions to diagnose features and/or compare matched features:

    \code{FeatureAxis}, \code{FeatureComps}, \code{FeatureMatchAnalyzer}, \code{FeatureProps}, \code{FeatureTable}, \code{interester}

See \code{compositer} for setting up composited objects, and see \code{hiw} (along with \code{distill} and \code{summary} method functions) for some shape analysis functionality.

The cluster analysis methods of Marzban and Sandgathe (2006; 2008) have been added.  The former method was written from scratch by Eric Gilleland 

\code{clusterer}

and the latter variation was modified from code originally written by Hilary Lyons 

\code{CSIsamples}.

The Structure, Amplitude and Location (SAL) method can be performed with \code{saller}.

2d. Geometrical characterization measures:

Perhaps the measures in this sub-heading are best described as part-and-parcel of 2c.  They are certainly useful in that domain, but have been proposed also for entire fields by AghaKouchak et al. (2011); though similar measures have been applied in, e.g., MODE.  The measures introduced in AghaKouchak et al. (2011) available here are: connectivity index (Cindex), shape index (Sindex), and area index (Aindex):

\code{Cindex}, \code{Sindex}, \code{Aindex}

3. Statistical inferences for spatial (and/or spatiotemporal) fields:

In addition to the methods categorized in Gilleland et al. (2009), there are also functions for making comparisons between two spatial fields.  The field significance approach detailed in Elmore et al. (2006), which requires a temporal dimension as well, involves using a circular block bootstrap (CBB) algorithm (usually for the mean error) at each grid point (or location) individually to determine grid-point significance (null hypothesis that the mean error is zero), and then a semi-parametric Monte Carlo method viz. Livezey and Chen (1983) to determine field significance.

\code{spatbiasFS}, \code{LocSig}, \code{MCdof}

In addition, the spatial prediction comparison test (SPCT) introduced by Hering and Genton (2011) is included via the functions: \code{lossdiff}, \code{empiricalVG} and \code{flossdiff}.  Supporting functions for calculating the loss functions include: absolute error (\code{abserrloss}), square error (\code{sqerrloss}) and correlation skill (\code{corrskill}), as well as the distance map loss function (\code{distmaploss}) introduced in Gilleland (2013).

4. Other:

The bias corrected TS and ETS (or TS dHdA and ETS dHdA) introduced in Mesinger (2008) are now included within the \code{vxstats} function.

The 2-d Gaussian Mixture Model (GMM) approach introduced in Lakshmanan and Kain (2010) can be carried out using the

\code{gmm2d}

function (to estimate the GMM) and the associated \code{summary} function calculates the parameter comparisons.  Also available are \code{plot} and \code{predict} method functions, but it can be very slow to run.  The \code{gmm2d} employs an initialization function that takes the K largest object areas (connected components) and uses their centroids as initial estimates for the means, and uses the axes as initial guesses for the standard deviations.  However, the user may supply their own initial estimate function.

The S1 score and anomaly correlation (ACC) are available through the functions

\code{S1} and \code{ACC}.

See Brown et al. (2012) and Thompson and Carter (1972) for more on these statistics.

Also included is a function to do the geographic box-plot of Willmott et al. (2007).  Namely,

\code{GeoBoxPlot}.

Datasets:

All of the initial Spatial Forecast Verification Inter-Comparison Project (ICP, \url{https://ral.ucar.edu/projects/icp/}) data sets used in the special collection of the Weather and Forecasting journal are included.  See the help file for 

\code{obs0426},

which gives information on all of these datasets that are included, as well as two examples for plotting them: one that does not preserve projections, but plots the data without modification, and another that preserves the projections, but possibly with some interpolative smoothing.

Ebert (2008) provides a nice review of these methods.  Roberts and Lean (2008) describes one of the methods, as well as the primary boxcar kernel smoothing method used throughout this package.  Gilleland et al. (2009, 2010) provides an overview of most of the various recently proposed methods, and Ahijevych et al. (2009) describes the data sets included in this package.  Some of these have been applied to the ICP test cases in Ebert (2009).

Additionally, one of the NIMROD cases (as provided by the UK Met Office) analyzed in Casati et al (2004) (case 6) is included along with approximate lon/lat locations.  See the help file for UKobs6 more information.

A spatio-temporal verification dataset is also included for testing the method of Elmore et al. (2006).  See the help file for

\code{GFSNAMfcstEx}.

A simulated dataset similar to the one used in Marzban and Sandgathe (2010) is also available and is called 

\code{hump}.
}

\author{
Eric Gilleland
}
\references{
AghaKouchak, A., Nasrollahi, N. Li, J. Imam, B. and Sorooshian, S. (2011) Geometrical characterization of precipitation patterns.  \emph{J. Hydrometeorology}, \bold{12}, 274--285, doi:10.1175/2010JHM1298.1.

Ahijevych, D., Gilleland, E., Brown, B. G. and Ebert, E. E. (2009) Application of spatial verification methods to idealized and NWP gridded precipitation forecasts. \emph{Wea. Forecasting}, \bold{24} (6), 1485--1497.

Baddeley, A. J. (1992) An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, Forstner, W. and Ruwiedel, S. Eds., Wichmann, 59--78.

Blaschke, T., Lang, S. and  Hay, G. (Eds.) (2008)  Object-Based Image Analysis.  Berlin, Germany: Springer-Verlag, 818 pp.

Briggs, W. M. and Levine, R. A. (1997) Wavelets and field forecast verification. \emph{Mon. Wea. Rev.}, \bold{125}, 1329--1341.

Brown, B. G., Gilleland, E. and Ebert, E. E. (2012) Chapter 6: Forecasts of spatial fields. pp. 95--117, In \emph{Forecast Verification: A Practitioner's Guide in Atmospheric Science}, 2nd edition. Edts. Jolliffe, I. T. and Stephenson, D. B., Chichester, West Sussex, U.K.: Wiley, 274 pp.

Casati, B. (2010) New Developments of the Intensity-Scale Technique within the Spatial Verification Methods Inter-Comparison Project. \emph{Wea. Forecasting} \bold{25}, (1), 113--143, doi:10.1175/2009WAF2222257.1.

Casati, B., Ross, G. and Stephenson, D. B. (2004) A new intensity-scale approach for the verification of spatial precipitation forecasts. \emph{Meteorol. Appl.} \bold{11}, 141--154.

Davis, C. A., Brown, B. G. and Bullock, R. G. (2006) Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784.

Ebert, E. E. (2008) Fuzzy verification of high resolution gridded forecasts: A review and proposed framework.  \emph{Meteorol. Appl.}, \bold{15}, 51--64. DOI: 10.1002/met.25 

Ebert, E. E. (2009) Neighborhood verification: A strategy for rewarding close forecasts.  \emph{Wea. Forecasting}, \bold{24}, 1498--1510, doi:10.1175/2009WAF2222251.1.

Elmore, K. L., Baldwin, M. E. and Schultz, D. M. (2006) Field significance revisited: Spatial bias errors in forecasts as applied to the Eta model.  \emph{Mon. Wea. Rev.}, \bold{134}, 519--531.

Gilleland, E., 2020. Novel measures for summarizing high-resolution forecast performance. Submitted to \emph{Advances in Statistical Climatology, Meteorology and Oceanography} on 19 July 2020.

Gilleland, E. (2013) Testing competing precipitation forecasts accurately and efficiently: The spatial prediction comparison test.  \emph{Mon. Wea. Rev.}, \bold{141}, (1), 340--355.

Gilleland, E. (2011)  Spatial forecast verification: Baddeley's delta metric applied to the ICP test cases.  \emph{Wea. Forecasting}, \bold{26}, 409--415, doi:10.1175/WAF-D-10-05061.1.

Gilleland, E., Lee, T. C. M., Halley Gotway, J., Bullock, R. G. and Brown, B. G. (2008) Computationally efficient spatial forecast verification using Baddeley's delta image metric.  \emph{Mon. Wea. Rev.}, \bold{136}, 1747--1757.

Gilleland, E., Ahijevych, D., Brown, B. G., Casati, B. and Ebert, E. E. (2009) Intercomparison of Spatial Forecast Verification Methods. \emph{Wea. Forecasting}, \bold{24}, 1416--1430, doi:10.1175/2009WAF2222269.1.

Gilleland, E., Ahijevych, D. A., Brown, B. G. and Ebert, E. E. (2010) Verifying Forecasts Spatially. \emph{Bull. Amer. Meteor. Soc.}, October, 1365--1373.

Harris, D., Foufoula-Georgiou, E., Droegemeier, K. K. and Levit, J. J. (2001)  Multiscale statistical properties of a high-resolution precipitation forecast.  \emph{J. Hydrometeorol.}, \bold{2}, 406--418.

Hering, A. S. and Genton, M. G. (2011) Comparing spatial predictions.  \emph{Technometrics} \bold{53}, (4), 414--425.

Lack, S., Limpert, G. L. and Fox, N. I. (2010) An object-oriented multiscale verification scheme. \emph{Wea. Forecasting}, \bold{25}, 79--92, DOI: 10.1175/2009WAF2222245.1

Lakshmanan, V. and Kain, J. S. (2010) A Gaussian Mixture Model Approach to Forecast Verification. \emph{Wea. Forecasting}, \bold{25} (3), 908--920.

Livezey, R. E. and Chen, W. Y. (1983) Statistical field significance and its determination by Monte Carlo techniques.  \emph{Mon. Wea. Rev.}, \bold{111}, 46--59.

Lucas, B D. and Kanade, T. (1981)  An iterative image registration technique with an application to stereo vision.  \emph{Proc. Imaging Understanding Workshop}, DARPA, 121--130.

Marzban, C. and Sandgathe, S. (2006) Cluster analysis for verification of precipitation fields.  \emph{Wea. Forecasting}, \bold{21}, 824--838.

Marzban, C. and Sandgathe, S. (2008) Cluster Analysis for Object-Oriented Verification of Fields: A Variation. \emph{Mon. Wea. Rev.}, \bold{136}, (3), 1013--1025.

Marzban, C. and Sandgathe, S. (2009) Verification with variograms.  \emph{Wea. Forecasting}, \bold{24} (4), 1102--1120, doi: 10.1175/2009WAF2222122.1

Marzban, C. and Sandgathe, S. (2010) Optical flow for verification.  \emph{Wea. Forecasting}, \bold{25}, 1479--1494, doi:10.1175/2010WAF2222351.1.

Mesinger, F. (2008) Bias adjusted precipitation threat scores.  \emph{Adv. Geosci.}, \bold{16}, 137--142.

Micheas, A. C., Fox, N. I., Lack, S. A. and Wikle, C. K. (2007) Cell identification and verification of QPF ensembles using shape analysis techniques. \emph{J. of Hydrology}, \bold{343}, 105--116.

Nachamkin, J. E. (2004) Mesoscale verification using meteorological composites. \emph{Mon. Wea. Rev.}, \bold{132}, 941--955.

Roberts, N. M. and Lean, H. W. (2008) Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events.  \emph{Mon. Wea. Rev.}, \bold{136}, 78--97. doi:10.1175/2007MWR2123.1.

Schwedler, B. R. J. and Baldwin, M. E. (2011)  Diagnosing the sensitivity of binary image measures to bias, location, and event frequency within a forecast verification framework.  \emph{Wea. Forecasting}, \bold{26}, 1032--1044, doi:10.1175/WAF-D-11-00032.1.

Thompson, J. C. and Carter, G. M. (1972) On some characteristics of the S1 score.  \emph{J. Appl. Meteorol.}, \bold{11}, 1384--1385.

Venugopal, V., Basu, S. and Foufoula-Georgiou, E. (2005) A new metric for comparing precipitation patterns with an application to ensemble forecasts.  \emph{J. Geophys. Res.}, \bold{110}, D08111, doi:10.1029/2004JD005395, 11pp.

Wernli, H., Paulat, M., Hagen, M. and Frei, C. (2008)  SAL--A novel quality measure for the verification of quantitative precipitation forecasts.  \emph{Mon. Wea. Rev.}, \bold{136}, 4470--4487.

Wernli, H., Hofmann, C. and Zimmer, M. (2009)  Spatial forecast verification methods intercomparison project: Application of the SAL technique.  \emph{Wea. Forecasting}, \bold{24}, 1472--1484, doi:10.1175/2009WAF2222271.1

Willmott, C. J., Robeson, S. M. and Matsuura, K. (2007)  Geographic box plots.  \emph{Physical Geography}, \bold{28}, 331--344, DOI: 10.2747/0272-3646.28.4.331.

Zhu, M., Lakshmanan, V. Zhang, P. Hong, Y. Cheng, K. and Chen, S. (2011) Spatial verification using a true metric.  \emph{Atmos. Res.}, \bold{102}, 408--419, doi:10.1016/j.atmosres.2011.09.004.
}
\examples{
## See help files for above named functions and datasets
## for specific examples.
}
\keyword{ package}
