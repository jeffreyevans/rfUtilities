#' @title Plot random forests model selection
#' @description Dot plot function for rf.modelSel importance values 
#'
#' @param  x      A rf.modelSel object
#' @param  imp    Plot selected ("sel") or all ("all") importance used in model selection
#' @param  ...    Additional arguments passed to plot
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'                    
#' @method plot rf.modelSel
#'
#' @export    	     
plot.rf.modelSel <- function(x, imp = c("all", "sel"),  ...) {
  plot.ms <- function(x, ...) {	  
    dots <- as.list(match.call(expand.dots = TRUE)[-1])
	  n <- x$parameter 
        x <- as.matrix( x$importance )
        rownames(x) <- n		
      ord <- rev(order(x[,1], decreasing=TRUE)[1:nrow(x)])
    dots[["x"]] <- x[ord,1]
    if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <-  lable
    if (is.null(dots[["pch"]]) & "pch" %in% names(dots) == FALSE) dots[["pch"]] <-  20
  	do.call("dotchart", dots)
  }
    if( imp[1] == "sel" ) { 
	  imp = x$importance[which(x$importance$parameter %in% x$selvars),]
	} else if (imp[1] == "all") { 
	  imp = x$importance 
	}
      if (x$scaling[1]=="mir") {lable = "Row Standardization Variable Importance"} 	
    if (x$scaling[1]=="se") {lable = "Standardized Error Variable Importance"}
  plot.ms(imp, ...) 
}
