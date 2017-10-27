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
  for(i in 1:length(x$parameters)) {
    cat("Variables in parameter set", i, "\n")
      cat("\t", x$parameters[[i]], "\n")
  	cat("\n")
  }
  cat("Variable importance for selected parameters:", "\n")
    cat("\n")
    print(x$importance)  
  	cat("\n")
	
  cat("Variable importance test for selected parameters:", "\n")
    cat("\n")
	
    imp <- data.frame(var=rownames(x$importance), imp = x$importance$imp) 
    imp <- imp[order(-imp$imp),]	
	p <- as.vector(t(x$parameters[[as.numeric(1)]]))
    parameters <- as.data.frame(t(p[match(imp$var, p)]))
      for(i in 2:length(x$parameters)) {
        p <- as.vector(t(x$parameters[[as.numeric(i)]]))
        p <- as.data.frame(t(p[match(imp$var, p)]))
  	    parameters <- merge(parameters, p, all=TRUE)
	  }
	parameters <- parameters[ as.numeric(rownames(x$test)) ,]
      names(parameters) <- paste0("x", 1:ncol(parameters))
	print( data.frame(x$test, parameters) )  
    cat("\n")	
	
  if( "rf.final" %in% ls(x) ) {
    cat("##################################", "\n")
    cat("Selected random forests model:", "\n")
    print(x$rf.final)
  }  
}
 