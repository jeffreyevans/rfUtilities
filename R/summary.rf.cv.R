#' @title Summarizing cross-validation
#' @description Summarizing of the rf.crossValidation function 
#'
#' @param object  Object of class rf.cv
#' @param ... Ignored
#'
#' @method summary rf.cv
#'
#' @export
summary.rf.cv <- function(object, ...) {
  if(class(object)[2] == "classification") { 
    cat("Classification accuracy for cross-validation", "\n")
    cv <- data.frame()	
      cv <- rbind(cv,
              apply(object$cross.validation$cv.users.accuracy, MARGIN = 2, stats::median),
              apply(object$cross.validation$cv.producers.accuracy, MARGIN = 2, stats::median))
    	      row.names(cv) <- c("users.accuracy", "producers.accuracy")
                names(cv) <- names(object$cross.validation$cv.users.accuracy)
      cat("", "\n")
    print( cv )
    cat("", "\n")
    cat("Cross-validation Kappa", "=", stats::median(object$cross.validation$cv.oob[,"kappa"]), "\n")
    cat("Cross-validation OOB Error", "=", stats::median(object$cross.validation$cv.oob[,"OOB"]), "\n")
    cat("Cross-validation error variance", "=", stats::var(object$cross.validation$cv.oob[,"OOB"]), "\n")
      cat("", "\n")
    cat("", "\n")
    
    cat("Classification accuracy for model", "\n")
    mdl <- data.frame()	
      mdl <- rbind(mdl,
              apply(object$model$model.users.accuracy, MARGIN = 2, stats::median),
              apply(object$model$model.producers.accuracy, MARGIN = 2, stats::median))
    	      row.names(mdl) <- c("users.accuracy", "producers.accuracy")
                names(mdl) <- names(object$model$model.users.accuracy)
      cat("", "\n")
    print( mdl )
    cat("", "\n")
      cat("Model Kappa", "=", stats::median(object$model$model.oob[,"kappa"]), "\n")
    cat("Model OOB Error", "=", stats::median(object$model$model.oob[,"OOB"]), "\n")
    cat("Model error variance", "=", stats::var(object$model$model.oob[,"OOB"]), "\n")
  }
  if(class(object)[2] == "regression") {
    cat("Fit MSE", "=", object[["fit.mse"]], "\n")
    cat("Fit percent variance explained", "=", object[["fit.var.exp"]], "\n")
  	cat("Median permuted MSE", "=", stats::median(object[["model.mse"]]), "\n")
	cat("Median permuted percent variance explained", "=", stats::median(object[["model.varExp"]]), "\n")
	cat("Median permuted percent variance explained", "=", stats::median(object[["model.varExp"]]), "\n")	
	cat("Median cross-validation RMSE", "=", stats::median(object[["y.rmse"]]), "\n")
    cat("Median cross-validation MBE", "=", stats::median(object[["y.mbe"]]), "\n")
	cat("Median cross-validation MAE", "=", stats::median(object[["y.mae"]]), "\n")	
	cat("RMSE cross-validation error variance", "=", stats::var(object[["y.rmse"]]), "\n")
    cat("MBE cross-validation error variance", "=", stats::var(object[["y.mbe"]]), "\n")
	cat("MAE cross-validation error variance", "=", stats::var(object[["y.mae"]]), "\n")
  }
}
