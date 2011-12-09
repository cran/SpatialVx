\name{SpatialVx-package}
\alias{SpatialVx-package}
\alias{SpatialVx}
\docType{package}
\title{
Spatial Forecast Verification
}
\description{
\pkg{SpatialVx} contains functions to perform several spatial forecast verification methods.  this release has functionality to do most of the traditional and neighborhood smoothing methods, as well as FQI and several binary image measures/metrics.  Some feature-based methods functions are also included.}

\details{

Primary functions include:

1. Filter Methods:

1a. Neighbordoohd Methods:

Neighborhood methods generally apply a convolution kernel smoother to one or both of the fields in the verificaiton set, and then apply the traditional scores.  Most of the methods reviewe in Ebert (2008, 2009) are included in this package.  The main functions are:

\code{hoods2dPrep}, \code{hoods2d}, \code{pphindcast2d}, \code{kernel2dsmooth}, and \code{plot.hoods2d}.

1b. Scale Separation Methods:

Scale separation refers to the idea of applying a band-pass filter (and/or doing a multi-resolution analysis, MRA) to the verification set.  Typically, skill is assessed on a scale-by-scale basis.  However, other techniques are also applied.  For example, denoising the field before applying traditional statistics, using the variogram, or applying a statistical test based on the variogram (these last are less similar to the spirit of the scale separation idea, but are at least somewhat related).

There is functionality to do the wavelet methods proposed in Briggs and Levine (1997).  In particular, to simply denoise the fields before applying traditional verification statistics, use \code{wavePurifyVx}.  To apply verification statistics to detail fields (in either the wavelet or field space), use \code{waverify2d} (dyadic fields) or \code{mowaverify2d} (non-dyadic or dyadic) fields.

The intensity-scale technique introduced in Casati et al. (2004) and the new developments of the technique proposed in Casati (2010) can be performed with \code{waveIS}.

2. Displacement Methods:  

In Gilleland et al (2009), this category was broken into two main types as field deformation and features-based.  The former lumped together binary image measures/metrics with field deformation techniques because the binary image measures inform about the "similarity" (or dissimilarity) between the spatial extent of two fields (across the entire field).  Here, they are broken down further into those that yield a only a single (or small vector of) metric(s) or measure(s) (location measures), and those that have mechanisms for moving grid-point locations to match the fields better spatially (field deformation).

2a. Location Measures: 

Included here are the well-known Hausdorff metric, the partial-Hausdorff measure, FQI (Venugopal et al., 2005), Baddeley's delta metric (Gilleland, 2011; Schwedler et al., 2011), metrV (Zhu et al., 2011), as well as the localization performance measures described in Baddeley, 1992: mean error distance, mean square error distance, and Pratt's Figure of Merit (FOM).

\code{locmeasures2dPrep}, \code{locmeasures2d}, \code{metrV}, \code{distob}, \code{locperf}

2b. Field deformation:

coming soon.  An image warping pacakge is in the works, and will be released as a separate package, but included as a dependency by this package potentially with some wrapper functions specific to spatial forecast verification needs.  Other field deformation aproaches are also in the hopper, and will be made available as soon as they are ready.

2c. Features-based methods: These methods are also sometimes called object-based methods, and have many similarities to techniques used in Object-Based Image Analysis (OBIA), a relatively new research area that has emerged primarily as a result of advances in earth observations sensors and GIScience (Blaschke et al., 2008).  It is attempted to identify individual features within a field, and subsequently analyze the fields on a feature-by-feature basis.  This may involve intensity error information in addition to location-specific error information.  Additionally, contingency table verifcation statistics can be found using new definitions for hits, misses and false alarms (correct negatives are more difficult to asses, but can also be done).

Currently, there is some functionality for performing the analyses introduced in Davis et al. (2006,2009), including the merge/match algorithm of Gilleland et al (2008), as well as the SAL technique of Wernli et al (2008, 2009).

\code{FeatureSuitePrep}, \code{FeatureFunPrep}, \code{FeatureSuite}, \code{saller}, \code{deltamm}, \code{convthresh}, \code{FeatureMatchAnalyzer}

3. Spatio-temporal inferences:

In addition to the methods categorized in Gilleland et al. (2009), there are also functions for making comparisons between two spatial fields.  Currently, this entails only one such method (more to come), whcih requires temporal information as well.  It is the field significance approach detailed in Elmore et al. (2006).  It involves using a circular block bootstrap (CBB) algorithm (usually for the mean error) at each grid point (or location) individually to determine grid-point significance (null hypothesis that the mean error is zero), and then a semi-parametric Monte Carlo method viz. Livezey and Chen (1983) to determine field significance.

\code{spatbiasFS}, \code{LocSig}, \code{MCdof}

4. Other:

Also included is a function to do the geographic box-plot of Willmott et al. (2007).  See the help file for \code{GeoBoxPlot}.

Datasets:

All of the initial Spatial Forecast Verification Inter-Comparison Project (ICP, http://www.ral.ucar.edu/projects/icp) data sets used in the special collection of the Weather and Forecasting journal are included.  See the help file for 'obs0426', which will give information on all of the datasets included, as well as an example for plotting the data.

Ebert (2008) provides a nice review of these methods.  Roberts and Lean (2008) describes one of the methods, as well as the primary boxcar kernel smoothing method used throughout this package.  Gilleland et al. (2009, 2010) provides an overview of most of the various recently proposed methods, and Ahijevych et al. (2009) describes the data sets included in this package.  Some of these have been applied to the ICP test cases in Ebert (2009).

Additionally, one of the NIMROD cases (as provided by the UK Met Office) analyzed in Casati et al (2004) (case 6) is included along with approximate lon/lat locations.  See the help file for UKobs6 more information.

A spatio-temporal verification dataset is also included for testing the method of Elmore et al. (2006).  See the help file for 'GFSNAMfcstEx' for more information.
}

\author{
Eric Gilleland
}
\references{
Ahijevych, D., E. Gilleland, B.G. Brown, and E.E. Ebert, 2009. Application of spatial verification methods to idealized and NWP gridded precipitation forecasts. \emph{Wea. Forecasting}, \bold{24} (6), 1485--1497.

Baddeley, A. J., 1992: An error metric for binary images.  In \emph{Robust Computer Vision Algorithms}, W. Forstner and S. Ruwiedel, Eds., Wichmann, 59--78.  Available at: http://www.google.com/url?sa=t&rct=j&q=an%20error%20metric%20for%20binary%20images&source=web&cd=1&ved=0CB4QFjAA&url=http%3A%2F%2Fciteseerx.ist.psu.edu%2Fviewdoc%2Fdownload%3Fdoi%3D10.1.1.52.3879%26rep%3Drep1%26type%3Dpdf&ei=GdXoTrSRJaSIsQLh1tiYCQ&usg=AFQjCNE2Ngehh_J2dXIsonPCMAvAL-fz_g&cad=rja

Blaschke, T., S. Lang, and  G. Hay (Eds.), 2008.  Object-Based Image Analysis.  Springer-Verlag, Berlin, Germany, 818 pp.

Briggs, W. M. and R. A. Levine, 1997. Wavelets and field forecast verification. \emph{Mon. Wea. Rev.}, \bold{125}, 1329--1341.

Casati, B., 2010: New Developments of the Intensity-Scale Technique within the Spatial Verification Methods Inter-Comparison Project. \emph{Wea. Forecasting} \bold{25}, (1), 113--143, DOI: 10.1175/2009WAF2222257.1.

Casati B, G Ross, and DB Stephenson, 2004. A new intensity-scale approach for the verification of spatial precipitation forecasts. \emph{Meteorol. Appl.} \bold{11}, 141--154.

Davis CA, BG Brown, and RG Bullock, 2006. Object-based verification of precipitation forecasts, Part I: Methodology and application to mesoscale rain areas. \emph{Mon. Wea. Rev.}, \bold{134}, 1772--1784

Ebert EE, 2008. Fuzzy verification of high resolution gridded forecasts: A review and proposed framework.  \emph{Meteorol. Appl.}, \bold{15}, 51--64. DOI: 10.1002/met.25 Available at http://www.ecmwf.int/newsevents/meetings/workshops/2007/jwgv/METspecialissueemail.pdf

Ebert, E. E., 2009: Neighborhood verification: A strategy for rewarding close forecasts.  \emph{Wea. Forecasting}, \bold{24}, 1498-1510, DOI: 10.1175/2009WAF2222251.1.

Elmore, K. L., M. E. Baldwin, and D. M. Schultz, 2006: Field significance revisited: Spatial bias errors in forecasts as applied to the Eta model.  \emph{Mon. Wea. Rev.}, \bold{134}, 519--531.

Gilleland, E., 2011.  Spatial forecast verification: Baddeley's delta metric applied to the ICP test cases.  \emph{Wea. Forecasting}, \bold{26}, 409--415, DOI: 10.1175/WAF-D-10-05061.1.

Gilleland, E., T. C. M. Lee, J. Halley Gotway, R. G. Bullock, and B. G. Brown, 2008: Computationally efficient spatial forecast verification using Baddeley's delta image metric.  \emph{Mon. Wea. Rev.}, \bold{136}, 1747--1757.

Gilleland, E., D. Ahijevych, B.G. Brown, B. Casati, and E.E. Ebert, 2009. Intercomparison of Spatial Forecast Verification Methods. \emph{Wea. Forecasting}, \bold{24}, 1416--1430, DOI: 10.1175/2009WAF2222269.1.

Gilleland, E., D.A. Ahijevych, B.G. Brown and E.E. Ebert, 2010: Verifying Forecasts Spatially. \emph{Bull. Amer. Meteor. Soc.}, October, 1365--1373.

Livezey, R. E. and W. Y. Chen, 1983: Statistical field significance and its determination by Monte Carlo techniques.  \emph{Mon. Wea. Rev.}, \bold{111}, 46--59.

Roberts, N. M. and H. W. Lean, 2008: Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events.  \emph{Mon. Wea. Rev.}, \bold{136}, 78--97. DOI: 10.1175/2007MWR2123.1.

Schwedler, B. R. J. and M. E. Baldwin, 2011.  Diagnosing the sensitivity of binary image measures to bias, location, and event frequency within a forecast verification framework.  Accepted to \emph{Wea. Forecasting}, doi: 10.1175/WAF-D-11-00032.1.

Venugopal, V., S. Basu, and E. Foufoula-Georgiou, 2005: A new metric for comparing precipitation patterns with an application to ensemble forecasts.  \emph{J. Geophys. Res.}, \bold{110}, D08111, doi:10.1029/2004JD005395, 11pp.

Wernli, H., M. Paulat, M. Hagen, and C. Frei, 2008.  SAL--A novel quality measure for the verification of quantitative precipitation forecasts.  \emph{Mon. Wea. Rev.}, \bold{136}, 4470--4487.

Wernli, H., C. Hofmann, and M. Zimmer, 2009.  Spatial forecast verification methods intercomparison project: Application of the SAL technique.  \emph{Wea. Forecasting}, \bold{24}, 1472--1484, DOI: 10.1175/2009WAF2222271.1

Willmott, C. J., S. M. Robeson, and K. Matsuura, 2007.  Geographic box plots.  \emph{Physical Geography}, \bold{28}, 331--344, DOI: 10.2747/0272-3646.28.4.331.

Zhu, M., V. Lakshmanan, P. Zhang, Y. Hong, K. Cheng, and S. Chen, 2011: Spatial verification using a true metric.  \emph{Atmos. Res.}, \bold{102}, 408--419, doi:10.1016/j.atmosres.2011.09.004.
}
\examples{
## The example below may take a few minutes to run,
## and therefore is commented out in order for R's
## package checks to run in a shorter amount of time.
## Compare the results to Ebert (2009) Fig. 5.
##
## Note that 'levels' refer to number of grid points,
## which for this example are ~4 km, so that the
## plot labels differ by a factor of 4 from Ebert (2009).
## Future versions of this package will allow for different
## labelling.
\dontrun{
data( pert004)
data(pert000)
hold <- hoods2dPrep( "pert004", "pert000", thresholds=c(1, 2, 5, 10, 20, 50), levels=c(1, 3, 5, 9, 17, 33, 65, 129, 257), units="mm/h")
look <- hoods2d( hold, which.methods=c("fss", "multi.event"), verbose=TRUE)
plot( look)
look2 <- pphindcast2d( hold, verbose=TRUE)
look2

# Location measures.
data( geom000)
data(geom001)
hold <- locmeasures2dPrep("geom001", "geom000", thresholds=c(0.1,50.1), k=c(4,0.975), alpha=c(0.1,0.9), units="in/100")
hold2 <- locmeasures2d( hold)
summary( hold2)

}
}
\keyword{ package}
