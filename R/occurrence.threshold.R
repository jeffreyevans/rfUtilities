#' @title Test occurrence probability thresholds
#' @description A statistical sensitivity test for occurrence probability thresholds 
#'       
#' @param x             A classification randomForest model object  
#' @param class         What class to test, quoted
#' @param p             Vector of probability thresholds 
#' @param type          What statistic to use in evaluation ("delta.ss", "sum.ss", 
#'                      "kappa") 
#'
#' @return An "occurrence.threshold" class object containing a "thresholds" vector 
#'         object with evaluation statistic and probability thresholds as names.  
#'
#' @details
#' Available threshold evaluation statistics:
#' \itemize{
#' \item kappa - [The Kappa statistic is maximized]
#' \item sum.ss - [The sum of sensitivity and specificity is maximized]
#' \item delta.ss - [The absolute value of the difference between sensitivity and 
#'                   specificity is minimized]
#'  }
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references 
#' Jimenez-Valverde, A., & J.M. Lobo (2007). Threshold criteria for 
#'   conversion of probability of species presence to either-or 
#'   presence-absence. Acta Oecologica 31(3):361-369 
#' @references 
#' Liu, C., P.M. Berry, T.P. Dawson, R.G. Pearson (2005). Selecting 
#'   thresholds of occurrence in the prediction of species distributions. 
#'   Ecography 28:385-393.
#'
#' @examples
#' library(randomForest)
#' 
#'  data(imports85)
#'   imp85 <- imports85[,-2] 
#'   imp85 <- imp85[complete.cases(imp85), ]
#'   imp85[] <- lapply(imp85, function(x) if (is.factor(x)) x[, drop=TRUE] else x)
#'   y <- ifelse( imp85$numOfDoors != "four", "0", "1")   
#'   
#' ( rf.mdl <- randomForest(y = as.factor(y), x = imp85[,-5]) )
#'   ( delta.ss.t <- occurrence.threshold(rf.mdl, class = "1") )
#'   ( sum.ss.t <- occurrence.threshold(rf.mdl, class = "1", type = "sum.ss") ) 
#'   ( kappa.ss.t <- occurrence.threshold(rf.mdl, class = "1", type = "kappa") )
#'   
#'   par(mfrow=c(2,2))
#'     plot(sum.ss.t)
#'     plot(delta.ss.t)
#'     plot(kappa.ss.t)   
#' 
#' # using ranger
#' library(ranger) 
#' ( rf.mdl <- ranger(y = as.factor(y), x = imp85[,-5], probability = TRUE) )
#'  ( kappa.ss.t <- occurrence.threshold(rf.mdl, class = "1", type = "kappa") )
#'   plot(kappa.ss.t) 
#'
#' @exportClass occurrence.threshold 
#' @export  
occurrence.threshold <- function(x, class, p = seq(0.10, 0.9, 0.05), 
                                 type = "delta.ss") {
	if(missing(x)) 
	  stop( "x (randomForest object) is a required argument")
	if(missing(class)) 
	  stop( "class is a required argument") 	  
    if (!inherits(x, c("randomForest", "ranger"))) 
      stop("x is not randomForest class object")
    if(length(grep("~", x$call[[2]])) > 0)
      stop("does not support a formula interface, please use x, y arguments")
    # formating call and pulling data
    a <- as.list(x$call)[-1] 
      xdata = eval(a[[which(names(a) == "x")]]) 
      y = eval(a[[which(names(a) == "y")]]) 

	if(inherits(x, "randomForest")) {
	  if(tolower(x$type) != "classification")	
	    stop("Regression or unsupervised not supported \n")	
	  probs <- stats::predict(x, xdata, type="prob")
    } else if(inherits(x, "ranger")) {
	  if(tolower(x$treetype) != "probability estimation")	
	    stop("Not a probability forests, please use;  probability = TRUE  \n")
	  probs <- stats::predict(x,data = xdata)$predictions
	}	
    probs <- probs[,which(colnames(probs) %in% class)]
	   y <- ifelse( y == class, 1, 0) 
	    test.p <- vector()
    for(i in p) {
        y.p <- ifelse( probs >= i, 1, 0)
      if( type == "kappa") {	
        test.p <- append(test.p, accuracy(y.p, y)$kappa) 
	    }
	  else if (type == "delta.ss") {	
        test.p <- append(test.p, abs(accuracy(y.p, y)$sensitivity - 
		                 accuracy(y.p, y)$specificity)) 
	    }
	  else if(type == "sum.ss") {
	    test.p <- append(test.p, (accuracy(y.p, y)$sensitivity + 
		                 accuracy(y.p, y)$specificity))
      }
	}
    names(test.p) <- p
      s <- list(thresholds = test.p, statistic = type)
      class(s) <- c("occurrence.threshold")	  
    return( s )
  }
