#' @title Print for combined random forests ensembles 
#' @description print method for combined random forests ensembles 
#' @param x    Object of class rf.ensembles
#' @param ...  Ignored
#'
#' @method print rf.ensembles
#'
#' @export
"print.rf.ensembles" <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n")
  cat("               Type of random forest: ", x$type, "\n", sep="")
  cat("               Number of random forests models: ", x$nrf, "\n", sep="")  
  cat("                     Number of trees: ", x$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", x$mtry, "\n\n", sep="")
  if(x$type == "classification") {
    if(!is.null(x$confusion)) {
      cat("        OOB estimate of  error rate: ",
        round(x$err.rate*100,2), "%\n", sep="")
      cat("Confusion matrix:\n")
        print(x$confusion)
      if(!is.null(x$test$err.rate)) {
        cat("        Test set error rate: ",
          round(x$test$err.rate*100,2), "%\n", sep="")
        cat("Confusion matrix:\n")
          print(x$test$confusion)
      }
    }
  }
  if(x$type == "regression") {
    if(!is.null(x$mse)) {
      cat("          Mean of squared residuals: ", x$mse,
          "\n", sep="")
		  
      cat("                    % Var explained: ",
          round(x$rsq, digits=2), "\n", sep="")
		  
      if(!is.null(x$test$mse)) {
        cat("                       Test set MSE: ",
            round(x$test$mse, digits=2), "\n", sep="")
        cat("                    % Var explained: ",
            round(x$test$rsq, digits=2), "\n", sep="")
      }      
    }
    if (!is.null(x$coefs)) {
      cat("  Bias correction applied:\n")
      cat("  Intercept: ", x$coefs[1], "\n")
      cat("      Slope: ", x$coefs[2], "\n")
    }
  }
}
