#' @title Print random forests model selection
#' @description Print method for rf.modelSel objects
#'    
#' @param x    Object of class rf.modelSel
#' @param ...  Ignored
#'
#' @method print rf.modelSel
#'
#' @export
print.rf.modelSel <- function(x, ...) {
  cat("Selected variables:", "\n")
    cat("\t", x$selvars, "\n")
    cat("\n")
	
  for(i in 1:nrow(x$parameters)) {
    cat("Variables in parameter set", i, "\n")
      cat("\t", stats::na.omit(as.character(x$parameters[i,])), "\n")
  	cat("\n")
  }
  
  cat("Variable importance for all parameters:", "\n")
    cat("\n")
      print(x$importance)  
  	cat("\n")
	
  cat("Summary of competing parameter sets (models):", "\n")
    cat("\n")
      print(x$test)	
    cat("\n")	
	
  if( "rf.final" %in% ls(x) ) {
    cat("##################################", "\n")
    cat("Selected random forests model:", "\n")
    print(x$rf.final)
  }  
}
 