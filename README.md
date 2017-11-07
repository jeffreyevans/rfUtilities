rfUtilities
===========

R package for random forests model selection, class balance and validation

New release of of "rfUtilities" 2.1-2 includes new functions for calculating Log Loss performance evaluation a function implementing an Isotonic regression for calibration of the estimated posterior probabilities of a model. There is also a new function for deriving parameter effect size based on partial dependency (Cafri & Bailey, 2016). The statistics Mean Absolute Error (mae) and Mean Bias Error (mbe) were added to the rf.crossValidation function.  

Available functions in rfUtilities are:

          accuracy - A function, called by the rf.crossValidation function or independently, that provides validation statistics for    
                     binomial or regression models
          logLoss - Calculates Logarithmic loss (logLoss)
          multi.collinear - Multi-collinearity test with matrix permutation.
          occurrence.threshold - A statistical sensitivity test for occurrence probability thresholds
          probability.calibration - Isotonic probability calibration
          rf.class.sensitivity - Random Forests class-level sensitivity analysis
          rf.classBalance - Random Forests Class Balance (Zero Inflation Correction) Model
          rf.crossValidation - Random Forests classification or regression cross-validation
          rf.effectSize - Random Forests parameter effect size
          rf.imp.freq - Random Forests variable selection frequency
          rf.modelSel - Random Forests Model Selection
          rf.partial.ci - Random Forests regression partial dependency plot with confidence intervals
          rf.partial.prob - Random Forest probability scaled partial dependency plots
          rf.regression.fit - Evaluates fit and overfit of random forests regression models
          rf.significance - Significance test for classification or regression random forests models
