#' @title Summary for combined random forests ensembles 
#' @description summary method for combined random forests ensembles 
#' @param object     Object of class rf.ensembles
#' @param ...        Ignored
#'
#' @method summary rf.ensembles
#'
#' @export
"summary.rf.ensembles" <- function(object, ...) {
  cat("\nCall:\n", deparse(object$call), "\n")
  cat("               Type of random forest: ", object$type, "\n", sep="")
  cat("               Number of random forests models: ", object$nrf, "\n", sep="")  
  cat("                     Number of trees: ", object$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", object$mtry, "\n\n", sep="")
  if(object$type == "classification") {
    if(!is.null(object$confusion)) {
      cat("        OOB estimate of  error rate: ",
        round(object$err.rate*100,2), "%\n", sep="")
      cat("Confusion matrix:\n")
        print(object$confusion)
      if(!is.null(object$test$err.rate)) {
        cat("        Test set error rate: ",
          round(object$test$err.rate*100,2), "%\n", sep="")
        cat("Confusion matrix:\n")
          print(object$test$confusion)
      }
    }
  }
  if(object$type == "regression") {
    if(!is.null(object$mse)) {
      cat("          Mean of squared residuals: ", object$mse,
          "\n", sep="")
		  
      cat("                    % Var explained: ",
          round(object$rsq, digits=2), "\n", sep="")
		  
      if(!is.null(object$test$mse)) {
        cat("                       Test set MSE: ",
            round(object$test$mse, digits=2), "\n", sep="")
        cat("                    % Var explained: ",
            round(object$test$rsq, digits=2), "\n", sep="")
      }      
    }
    if (!is.null(object$coefs)) {
      cat("  Bias correction applied:\n")
      cat("  Intercept: ", object$coefs[1], "\n")
      cat("      Slope: ", object$coefs[2], "\n")
    }
  }
}
