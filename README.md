rfUtilities 2.1.5

[![CRAN
status](http://www.r-pkg.org/badges/version/rfUtilities)](https://cran.r-project.org/package=rfUtilities)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/grand-total/rfUtilities)](https://cran.r-project.org/package=rfUtilities)

*R package for random forests model selection, class balance and validation*

Random Forests Model Selection, inference, fit and performance evaluation

rfUtilities 2.1.5 (GitHub development release)

Release of of "rfUtilities" 2.1-4 includes bug fixes and function enhancement  

Release of of "rfUtilities" 2.1-3 includes new functions for combining random forests ensembles, bivariate partial dependence plots and unsupervised classification using random forests.   

Release of of "rfUtilities" 2.1-2 includes new functions for calculating Log Loss performance evaluation a function implementing an Isotonic regression for calibration of the estimated posterior probabilities of a model. There is also a new function for deriving parameter effect size based on partial dependency (Cafri & Bailey, 2016). The statistics Mean Absolute Error (mae) and Mean Bias Error (mbe) were added to the rf.crossValidation function. 

Available functions in rfUtilities are:

          accuracy - A function, called by the rf.crossValidation function or independently, that provides validation statistics for    
                     binomial or regression models
          bivariate.partialDependence - Bivariate partial-dependency plot
          logLoss - Calculates Logarithmic loss (logLoss)
          multi.collinear - Multi-collinearity test with matrix permutation.
          occurrence.threshold - A statistical sensitivity test for occurrence probability thresholds
          probability.calibration - Isotonic probability calibration
          rf.class.sensitivity - Random Forests class-level sensitivity analysis
          rf.classBalance - Random Forests Class Balance (Zero Inflation Correction) Model
          rf.combine - Combine Random Forests Ensembles
          rf.crossValidation - Random Forests classification or regression cross-validation
          rf.effectSize - Random Forests parameter effect size
          rf.imp.freq - Random Forests variable selection frequency
          rf.modelSel - Random Forests Model Selection
          rf.partial.ci - Random Forests regression partial dependency plot with confidence intervals
          rf.partial.prob - Random Forest probability scaled partial dependency plots
          rf.regression.fit - Evaluates fit and overfit of random forests regression models
          rf.significance - Significance test for classification or regression random forests models
          rf.unsupervised - Unsupervised Random Forests
          
**Bugs**: Users are encouraged to report bugs here. Go to [issues](https://github.com/jeffreyevans/rfUtilities/issues) in the menu above, and press new issue to start a new bug report, documentation correction or feature request. You can direct questions to <jeffrey_evans@tnc.org>.

**To install `rfUtilities` in R use install.packages() to download curent stable release from CRAN** 

**or, for the development version, run the following (requires the remotes package):**
`remotes::install_github("jeffreyevans/rfUtilities")`

**Tutorial**: See (http://evansmurphy.wixsite.com/evansspatial/random-forest-sdm).
