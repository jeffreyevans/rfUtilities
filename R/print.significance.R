#' @title Print significance
#' @description print method for class "significance"
#' @param x    Object of class significance
#' @param x    P-value to accept or reject significance
#' @param ...  Ignored
#'
#' @method print significance
#'
#' @export
print.significance <- function(x, p=0.05, ...) {
  cat("Number of permutations: ", x[["nPerm"]], "\n")
  cat("Randomization fraction: ", x[["rand.frac"]], "\n")
    if(tolower(x[["rf.type"]]) == "classification") {
      cat("\t", "p-value for OOB randomization: ", x[["oob.p.value"]], "\n")
	  cat("\t", "Use Kappa for OOB error: ", x[["use.kappa"]], "\n")
      cat("\t", "p-value forMcnemar's test test: ", x[["p.value"]], "\n")	
      cat("\t", "Model OOB error: ", x[["test.OOB"]], "\n")	    
      cat("\t", "Random OOB error: ", stats::median(x[["RandOOB"]]), "\n")
      cat("\t", "min random global error:", min(x[["RandOOB"]]), "\n")
      cat("\t", "max random global error: ",  max(x[["RandOOB"]]), "\n")
      cat("\t", "min random within class error:", min(x[["RandMaxError"]]), "\n")
      cat("\t", "max random within class error: ", max(x[["RandMaxError"]]), "\n")
      if(x[["oob.p.value"]] <= p) { accept = TRUE } else { accept = FALSE }		
    } else if(tolower(x[["rf.type"]]) == "regression") {
      cat("\t", "p-value for R-square randomization: ", x[["p.value"]], "\n")
      cat("\t", "Kolmogorov-Smirnov p-value under equal hypothesis: ", 
	      as.numeric(x$ks.p.value[1,][1]), "\n")	
      cat("\t", "Kolmogorov-Smirnov p-value under less-than hypothesis: ",
	      as.numeric(x$ks.p.value[1,][2]), "\n")
      cat("\t", "Kolmogorov-Smirnov p-value under greater-than hypothesis: ", 
	      as.numeric(x$ks.p.value[1,][3]), "\n")	  
      cat("\t", "Model R-square: ", x[["R.square"]], "\n")
      cat("\t", "Random R-square: ", stats::median(x[["RandR.square"]]), "\n")
      cat("\t", "Randomization variance: ", stats::var(x[["RandR.square"]]), "\n")
      if(x[["p.value"]] <= p) { accept = TRUE } else { accept = FALSE }
    }
    if(accept == TRUE) { 
	  cat("\t", "Model significant at p = ", round(p,5), "\n")
    } else if(accept == FALSE) { 
	  cat("\t", "Model not significant at p = ", round(p,5), "\n")
	}
}	
