#' @title Multi-collinearity test
#' @description Test for multi-collinearity in data using qr-matrix decomposition
#'
#' @param x             data.frame or matrix object
#' @param perm         (FALSE/TRUE) Should a permutation be applied
#' @param leave.out    (FALSE/TRUE) Should a variable be left out at each permutation
#' @param n             Number of permutations
#' @param p             multi-collinearity threshold
#' @param na.rm        (FALSE/TRUE) Remove NA values
#'
#' @return If perm == TRUE a data.frame of indicating the frequency that a variable was collinear and, if leave.out = TRUE 
#'         the number of times it was omitted. Otherwise, a vector of collinear variables is returned. If no colinear variables are identified a NULL is returned. 
#'
#' @author Jeffrey S. Evans  <jeffrey_evans<at>tnc.org>
#'
#' @references
#'  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole. 
#'  Dongarra, J. J., Bunch, J. R., Moler, C. B. and Stewart, G. W. (1978) LINPACK Users Guide. Philadelphia: SIAM Publications. 
#'
#' @note
#' A permutation approach is not available where, at each iteration, the columns are randomly rearranged and a parameter dropped. The frequency that a variable is identified as collinear is accumulated. 
#' The multi-collinearity threshold needs to be adjusted based on number of parameters. For small number(s) of variables (<20) use ~1e-07 and for larger ~0.05  
#'
#' @examples  
#' test <- data.frame(v1=seq(0.1, 5, length=100), v2=seq(0.1, 5, length=100), 
#'                    v3=dnorm(runif(100)), v4=dnorm(runif(100)))
#'						
#' # Single test					
#'   ( cl <- multi.collinear(test) )
#'
#' # Permutated test with leave out	
#' ( cl.test <- multi.collinear(test, perm = TRUE, leave.out = TRUE, n = 999) )
#'     cl.test[cl.test$frequency > 0,]$variables
#' 			       
#'  # Remove identified variable(s)
#'  head( test[,-which(names(test) %in% cl.test[cl.test$frequency > 0,]$variables)] )
#'
#' @export
multi.collinear <- function(x, perm = FALSE, leave.out = FALSE, n = 99, p = 1e-07, na.rm = FALSE) {
  if (!inherits(x, "data.frame") & !inherits(x, "matrix")) stop("x must be a data.frame or matrix")
    if ( (dim(x)[2] < 2) == TRUE) stop("Need at least two parameters to test")
      if(!inherits(x, "matrix")) x <- as.data.frame(x)
	    if(na.rm) { x <- stats::na.omit(x) }
    qrd <- function(x) {
      x <- as.matrix(x)
      n <- ncol(x)
      m <- nrow(x)
      q <- matrix(0, m, n)
      r <- matrix(0, n, n)  
        for (j in 1:n) {
          v = x[,j] 
          if (j > 1) {
            for (i in 1:(j-1)) {
              r[i,j] <- t(q[,i]) %*% x[,j]
              v <- v - r[i,j] * q[,i] 
            }      
          }
          r[j,j] <- sqrt(sum(v^2))
          q[,j] <- v / r[j,j]
        }  
      return( list('qr'=q, 'rank'=r) )
    }	  
	mc.test <- function(v, p) {
	  qrx <- qr(v, tol = p)
	  if (length(names(v)[qrx$pivot[1:qrx$rank]]) != length(v) ) {  
        keep <- names(v)[qrx$pivot[1:qrx$rank]]
      return(paste(setdiff(names(v), keep)))
	  }
	  return(NULL)
	}
    if(perm == TRUE) { 
	  freq <- data.frame(variables = names(x), frequency = 0)
	  if(leave.out == TRUE) freq <- data.frame(freq, leave.out = 0)
	    for(i in 1:n) {
	      x.data <- x[,sample(1:length(names(x)))]
	        if(leave.out == TRUE) { 
		      x.data <- x.data[,-sample(1:length(names(x.data)),1)]
		  	  lo.idx <- which(freq$variables %in% setdiff(names(x), names(x.data))) 
		  	  freq[lo.idx,]$leave.out <- freq[lo.idx,]$leave.out + 1
		    }			
		  mc <- mc.test(x.data, p = p) 
		    if(!is.null(mc)) {
              idx <- which(freq$variables %in% mc) 			
			  freq[idx,]$frequency <- freq[idx,]$frequency + 1
			}
		}
      return( freq )
    } else {
	  mc <- mc.test(x, p = p)
      if( length(mc) > 0) { 
	    return( mc )
      } else {
        return( NULL )
      }		
    }
}
