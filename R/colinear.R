#' @title Collinearity test
#' @description Test for collinearity in data 
#'
#' @param x    A symmetric correlation matrix
#' @param p    The correlation cutoff (default is 0.85)
#'
#' @return Messages and a vector of correlated variables 
#'
#' @author Jeffrey S. Evans  <jeffrey_evans<at>tnc.org>
#'
#' @note Evaluation of the pairwise correlated variables to remove is accomplished through 
#'       calculating the mean correlations of each variable and selecting the 
#'       variable with higher mean.  
#' 
#' @examples 
#' dat <- data.frame(v1=seq(0.1, 5, length=100), 
#'                   v2=seq(0.1, 5, length=100), 
#'                   v3=dnorm(runif(100)), 
#' 				     v4=dnorm(runif(100)))
#' 
#' # Evaluate collinearity
#' ( cor.vars <- collinear(cor(dat), p = 0.80) )
#'  			       
#' # Remove identified variable(s)
#' head( dat[,-which(names(dat) %in% cor.vars)] )
#' 
#' @export collinear
collinear <- function (x, p = 0.85) {
  if(!class(x) == "matrix")
    stop("x does not appear to be a matrix")    
  diag(x) <- 0
  if (!isTRUE(all.equal(x, t(x)))) 
    stop("correlation matrix is not symmetric")  
  if(is.null(colnames(x))) { 
    colnames(x) <- paste0("X", 1:ncol(x))
	rownames(x) <- paste0("X", 1:nrow(x))
  }	
  if (is.null(rownames(x))) 
    stop("'x' must have row names")		
  if (!any(x[!is.na(x)] > p))
    stop("All correlations are <=", p, "\n")
  if (dim(x)[1] < 2) 
    stop("There is only one variable")
      x2 <- abs(x)
        originalOrder <- 1:ncol(x2)
          averageCorr <- function(x) mean(x, na.rm = TRUE)
            tmp <- x2
              diag(tmp) <- NA
                maxAbsCorOrder <- order(apply(tmp, 2, averageCorr), 
				                        decreasing = TRUE)
                x2 <- x2[maxAbsCorOrder, maxAbsCorOrder]
              newOrder <- originalOrder[maxAbsCorOrder]
            rm(tmp)
      diag(x2) <- NA 
	  
	  combine.vars <- expand.grid(rownames(x2),colnames(x2))
	    combine.vars <- combine.vars[-which(combine.vars[,1] == combine.vars[,2]),]
	      deletecol <- rep(FALSE,length(colnames(x2)))
		    names(deletecol) <- colnames(x)		
    for (i in 1:nrow(x2)) {
      i.name = rownames(x2)[i]
      if( deletecol[grep(i.name, names(deletecol))][1] == TRUE ) {
        next
      }
      for(j in (1:ncol(x2))[-i]){
        idx.names = c(i.name,colnames(x2)[j])
        if( deletecol[grep(colnames(x2)[j], names(deletecol))][1] == TRUE ) {
          next
        }  	
          if(x2[i, j] > p) {
            rc1 <- mean(x2[i, ], na.rm = TRUE)
            cc2 <- mean(x2[-j, ], na.rm = TRUE)
    	      deletecol[grep(idx.names[which.max(c(rc1,cc2))], names(deletecol))] <- TRUE
              cat("Collinearity between ", idx.names[1], "and ", 
                idx.names[2], "correlation = ", round(x2[i,j], 4), "\n")
              cat("  Correlation means: ", round(rc1, 3), "vs", round(cc2, 3), "\n")	
              cat("   recommend dropping", idx.names[which.max(c(rc1,cc2))], "\n", "\n")			
            } 
        } 
    }  
  return(unique(names(deletecol[deletecol==TRUE])))
}
