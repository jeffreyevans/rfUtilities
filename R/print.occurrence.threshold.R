#' @title Print occurrence.threshold
#' @description Print method for occurrence.threshold objects
#'    
#' @param x    Object of class occurrence.threshold
#' @param ...  Ignored
#'
#' @method print occurrence.threshold
#'
#' @export
print.occurrence.threshold <- function(x, ...) {
  cat("Evaluation statistic:", x$statistic, "\n")
    cat("\n")  
    cat("Moments for thresholds:", "\n")
    cat("\n")  
      print( summary(x$thresholds) )
    cat("\n")  	
    if(x$statistic == "delta.ss") { 	
    message("Probability threshold with delta abs[sensitivity-specificity] = ", 
      names(which.min(x$thresholds)), "\n") 
	} else if (x$statistic == "sum.ss") {
    message("Probability threshold with cumulative [sensitivity+specificity] = ", 
      nnames(which.max(x$thresholds)), "\n")
    } else if (x$statistic == "kappa") {
    message("Probability threshold with maximum kappa = ", 
      names(which.max(x$thresholds)), "\n")		  
	} else if (x$statistic == "youden") {
    message("Probability threshold with of Youden's index [(sensitivity+specificity)-1] = ", 
      names(which.max(x$thresholds)), "\n")
	} else if (x$statistic == "logloss") {
    message("Probability threshold for minimization of logarithmic loss = ", 
      names(which.min(x$thresholds)), "\n")
    }
} 
