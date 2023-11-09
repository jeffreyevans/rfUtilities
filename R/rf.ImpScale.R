#' @title Scaling of Random Forests importance values
#' @description Various scaling approaches for scaling permuted parameter 
#'              importance values
#'              
#' @param x         randomForest or ranger object with a classification or 
#'                  regression instance 
#' @param scaling   Type of importance scaling, options are ("mir","se", "p")
#' @param n         If scaling = "p" how many permutations
#' @param sort      Sort output by importance
#'
#' @return A data.frame with:
#' \itemize{  
#'   \item   {"parameter"} {Name of the parameter}
#'   \item   {"importance"} {scaled importance values}
#'   \item   {"pvalue"} {Optional p-value for ranger objects}
#' }
#'
#' @details
#' The "mir" scale option performs a row standardization and the "se" option 
#' performs normalization using the "standard errors" of the permutation-based 
#' importance measure. Both options result in a 0-1 range but, "se" sums to 1.
#'
#' The scaled importance measures are calculated as: 
#'   mir = i/max(i) and se = (i / se) / ( sum(i) / se).
#' 
#' The "p" scale option calculates a p-value using the Janitza et al, (2018)
#'   method where, assuming that noise variables vary randomly around zero,
#'   uses negative values to build a null distribution.  
#' 
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @references 
#' Altmann, A., Tolosi, L., Sander, O. & Lengauer, T. (2010). Permutation importance: 
#'   a corrected feature importance measure, Bioinformatics 26:1340-1347. 
#' @references 
#' Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas 
#'  connectivity in Yellowstone National Park with landscape genetics. 
#'  Ecology 91:252-261
#' @references 
#' Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species 
#'  distribution and change using Random Forests CH.8 in Predictive Modeling 
#'  in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer
#' @references 
#' Janitza, S., E. Celik, and A.-L. Boulesteix (2018). A computationally fast 
#'  variable  importance test for random forests for high-dimensional data. 
#'  Advances in Data Analysis and Classification 12(4):885-915
#'
#' @examples
#' library(randomForest)
#' library(ranger)
#' 
#' #### Classification
#' data(iris)
#' iris$Species <- as.factor(iris$Species)
#' 
#' # Random Forests using x, y arguments with index
#' ( rf.class <- randomForest(x = iris[,1:4], y = iris[,"Species"], 
#'                            importance = TRUE) )
#'
#'   rf.ImpScale(rf.class, scaling="mir")
#'   rf.ImpScale(rf.class, scaling="se")
#' 
#' # ranger using formula interface 
#' ( rgr.class <- ranger(Species ~., data=iris, 
#'                       importance = "permutation") )
#'					   
#'   rf.ImpScale(rgr.class, scaling="mir")
#'   rf.ImpScale(rgr.class, scaling="p", n=10)
#' 
#' #### Regression
#' data(airquality)
#'   airquality <- na.omit(airquality)
#'   
#' ( rf.reg <- randomForest(x=airquality[,2:6], y=airquality[,1],
#'                          importance = TRUE) )
#'   rf.ImpScale(rf.reg, scaling="mir")
#'   rf.ImpScale(rf.reg, scaling="se")
#' 
#' ( rgr.reg <- ranger(x=airquality[,2:6], y=airquality[,1],
#'                     importance = "permutation") )
#'   rf.ImpScale(rgr.reg, scaling="mir")
#'   rf.ImpScale(rgr.reg, scaling="p", n=10)
#' 
#' @seealso \code{\link[ranger]{importance_pvalues}} details on Altmann's p-value method
#'
#' @export rf.ImpScale
rf.ImpScale <- function (x, scaling=c("mir","se", "p"), n=99, sort = FALSE) { 
  if (!any(sort(class(x))[1] == "randomForest") && !any(sort(class(x))[1] == "ranger")) 
    stop(deparse(substitute(x)), " Must be a randomForest or ranger object") 
  if(any(sort(class(x))[1] == "ranger")) {
    #if(length(grep("~", x$call[[2]])) > 0)
    #  stop("This package does not support a formula interface, 
	#        please use x, y arguments")
	if(scaling[1] == "se")
	  stop("Ranger does not support standard error scaling")
  if(x$importance.mode != "permutation")
     warning("It is highly recommend that you use the permuted importance")
	   mtype = tolower(x$treetype)  
  } else if(any(sort(class(x))[1] == "randomForest")) {
    mtype = tolower(x$type)
    #if( length(grep("~", x$call[[2]])) > 0 )
    #  stop("This package does not support a formula interface, please use x, y arguments")
	if(scaling[1] == "p")
	  stop("randomForest does not support importance p-values")  
  }
  #**** extract importance
  # regression 
  if (mtype == "regression") {   
    if(sort(class(x))[1] == "ranger") {  
      if (any(x$variable.importance == "none"))
        stop("ranger object does not contain importance") 
	  rf.imp <- x$variable.importance 
    } else if(sort(class(x))[1] == "randomForest") {	
      if (!"%IncMSE" %in% colnames(x$importance))
        stop("randomForest object does not contain importance") 
      rf.imp <- x$importance[,"%IncMSE"]
    }  	
  # classification 	
  } else if(any( mtype %in% c("classification", "unsupervised", 
                              "probability estimation"))) {
    if(sort(class(x))[1] == "ranger") {  
      if(any(x$importance == "none"))
        stop("ranger object does not contain importance") 
	  rf.imp <- x$variable.importance 
    } else if(sort(class(x))[1] == "randomForest") {	
      if (!"MeanDecreaseAccuracy" %in% colnames(x$importance))
        stop("randomForest object does not contain importance") 
      rf.imp <- x$importance[,"MeanDecreaseAccuracy"]
    } 	
  }  
  #**** importance scaling  
  if (scaling[1] == "mir") {
    i <- rf.imp / max(rf.imp) 
  } else if(scaling[1] == "se") {
      if(mtype == "regression") {
        rf.impSD <- x$importanceSD
	  } else if(mtype == "classification") {
	    rf.impSD <- x$importanceSD[,"MeanDecreaseAccuracy"]
	  }	
        rf.impSD[rf.impSD == 0] <- 0.000000001
      i <- ( rf.imp / rf.impSD ) / sum(rf.imp / rf.impSD, na.rm=TRUE)
  } else if(scaling[1] == "p") {
	a <- as.list(x$call)[-1]
    if(!"y" %in% names(a)) {
	  xy <- data.frame(eval(a[[which(names(a) == "data")]]))
	  fml = stats::formula(a[[1]], data = xy)	
	i <- ranger::importance_pvalues(x, num.permutations = n,
               method = "altmann", formula = fml, data = xy)
    } else {	
    xy <- data.frame(y=eval(a[[which(names(a) == "y")]]), 
                       eval(a[[which(names(a) == "x")]]))
      fml = stats::formula(y ~ ., data=xy)
	i <- ranger::importance_pvalues(x, num.permutations = n,
               method = "altmann", formula = fml, data = xy)
    }				   
  }    
  #**** format results  
  if(scaling[1] == "p") {
    i[,"importance"] <- i[,"importance"] / max(i[,"importance"])
      i <- data.frame(parameter = rownames(i), i)
        rownames(i) <- 1:nrow(i)	  
  } else {
    i <- data.frame(parameter = names(i),importance=i) 
	  rownames(i) <- 1:nrow(i)
  }
  return( if(sort) { i[order(i$importance),] } else { i } )
}
