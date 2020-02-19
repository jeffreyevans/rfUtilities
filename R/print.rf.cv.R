#' @title Print random forests cross-validation
#' @description Print method for rf.cv objects
#'    
#' @param x    Object of class rf.cv
#' @param ...  Ignored
#'
#' @method print rf.cv
#'
#' @export
print.rf.cv <- function(x, ...) {
  if(class(x)[2] == "classification") {
  cat("Classification accuracy for model", "\n")	
    mdl <- x$model.error[,1:2]
        rownames(mdl) <- paste0("Class ", names(x$cv.users.accuracy))
    cat("", "\n")
	print( t(mdl) )
	cat("", "\n")
    cat("Model Kappa", "=", x$model.error[,"model.kappa"][1], "\n")
	cat("Model OOB Error", "=", stats::median(x$model.error[,"model.oob"][1]), "\n") 
  cat("Classification accuracy for cross-validation", "\n")
  cv <- data.frame()	
    cv <- rbind(cv, apply(x$cv.users.accuracy, MARGIN = 2, stats::median, na.rm=TRUE),
                    apply(x$cv.producers.accuracy, MARGIN = 2, stats::median, na.rm=TRUE))
        row.names(cv) <- c("Cross-validation users.accuracy", "Cross-validation producers.accuracy")
      names(cv) <- paste0("Class ", names(x$cv.users.accuracy))
    cat("", "\n")
	print( cv )
	cat("", "\n")
	cat("Cross-validation Kappa", "=", stats::median(x$cv.oob[,"CV.kappa"],na.rm=TRUE), "\n")
	cat("Cross-validation Model OOB Error", "=", stats::median(x$cv.oob[,"Model.PCC"],na.rm=TRUE), "\n")
	cat("Cross-validation CV OOB Error", "=", stats::median(x$cv.oob[,"CV.PCC"],na.rm=TRUE), "\n")
	cat("Cross-validation error variance", "=", stats::var(x$cv.oob[,"CV.PCC"]), "\n")
    cat("", "\n")
  }
  if(class(x)[2] == "regression") {
    cat("Fit MSE", "=", x[["fit.mse"]], "\n")
    cat("Fit percent variance explained", "=", x[["fit.var.exp"]], "\n")
  	cat("Median permuted MSE", "=", stats::median(x[["model.mse"]]), "\n")
	cat("Median permuted percent variance explained", "=", stats::median(x[["model.varExp"]]), "\n")
	cat("Median cross-validation RMSE", "=", stats::median(x[["y.rmse"]]), "\n")
    cat("Median cross-validation MBE", "=", stats::median(x[["y.mbe"]]), "\n")
	cat("Median cross-validation MAE", "=", stats::median(x[["y.mae"]]), "\n")
    cat("Range of ks p-values", "=", range(x[["p.val"]]), "\n")
    cat("Range of ks D statistic", "=", range(x[["D"]]), "\n")	
	cat("RMSE cross-validation error variance", "=", stats::var(x[["y.rmse"]]), "\n")
    cat("MBE cross-validation error variance", "=", stats::var(x[["y.mbe"]]), "\n")
	cat("MAE cross-validation error variance", "=", stats::var(x[["y.mae"]]), "\n")
  }  
}
