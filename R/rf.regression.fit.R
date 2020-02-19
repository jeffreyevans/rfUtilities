#' @title Random Forest fit statistics  
#' @description Evaluates fit and overfit of random forests regression 
#' 
#' @param x    randomForest or ranger regression object
#' @param y    Vector of dependent variable used in model
#'
#' @return     A list and rf.fit class object with "fit" matrix of fit statistics 
#'             and "message" indicating overfit risk. 
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @examples
#'   library(randomForest)
#'   set.seed(131)
#'   data(airquality)
#'   airquality <- na.omit(airquality)
#'   
#'   ( rf.aq <- randomForest(airquality[,1:3], airquality[,"Ozone"]) )
#'   rf.regression.fit(rf.aq, airquality[,"Ozone"])
#'
#' @export
rf.regression.fit <- function(x, y) {
  if(!any(class(x) %in% c("randomForest","randomForest.formula","ranger")))
    stop("x is not a randomForest object")
  if(missing(y) | missing(x))
    stop("Both x and x arguments must be provided")
  if(any(class(x) %in% c("randomForest","randomForest.formula"))) {	
    if(!x$type == "regression") 
      stop("Classification models not supported")
	if(length(y) != length(x$predicted))
      stop("y does not match model")	 
        rmse <- sqrt(mean((x$predicted - y)^2))
          r2 <- mean(x$rsq)
            f2 <- r2 / (1 - r2)
          n <- length(y)
        k <- nrow(x$importance) 
      fit.Actuals.pred <- cbind(x$predicted, y)
  } else if(any(class(x) %in% c("ranger"))) {
    if(!x$treetype == "Regression") 
      stop("Classification models not supported") 
	if(length(y) != length(x$predictions))
      stop("y does not match model")		  
        rmse <- sqrt(mean((x$predictions - y)^2))
          r2 <- x$r.squared
            f2 <- r2 / (1 - r2)
          n <- x$num.samples
        k <- x$num.independent.variables
      fit.Actuals.pred <- cbind(x$predictions , y) 
  }  
  overfit.ratio <- round(n / k)  # overfitting ratio, ideal is > 10
  accuracy <- stats::median(apply(fit.Actuals.pred, 1, min) / 
                            apply(fit.Actuals.pred, 1, max))
    if(overfit.ratio > 9) { 
	  msg <- "Model is not overfit" 
       } else { 
	  msg <- "Model may be overfit" 
	}
	fit.stats <- t(data.frame(RMSE = round(rmse,3), 
	                          R.squared = round(r2,3), 
						      Cohen.f2 = round(f2,3),
                     	      Accuracy = round(accuracy,3), 
						      Overfitting.ratio = round(overfit.ratio,3)))
							  fit <- list(fit = fit.stats, message = msg)
    class(fit) <- c("rf.fit","list")	
  return(fit)
}
