#' @title Logarithmic loss (logLoss)
#' @description Evaluation of estimate quality in binomial models using cross-entropy or log likelihood loss 
#'    
#' @param p            vector of predicted probabilities {0-1} 
#' @param y            vector of observed binomial values {0,1}
#' @param likelihood   (FALSE/TRUE) return log likelihood loss, default is (FALSE) for log loss  
#' @param global       (TRUE/FALSE) return local or global log loss values, if FALSE local values are returned
#' @param eps          epsilon scaling factor to avoid NaN values
#'
#' @return If likelihood TRUE the log likelihood loss will be returned. If global FALSAE, a list with observed (y), probability (p) and log loss (log.loss) otherwise, a vector of global log loss value 
#' 
#' @note The log loss metric, based on cross-entropy, measures the quality of predictions rather than the accuracy.
#' @note Effectively, the log loss is a measure that gauges additional error comming the estimates as opposed to the true values.
#' @note As the estimated probability diverges from its observed value the log loss increases with an expected of [0-1] where 0 would be a perferct model.  
#' @note For a single sample with true value yt in {0,1} and estimated probability yp that yt = 1, the log loss is derived as: -log P(yt | yp) = -(yt log(yp) + (1 - yt) log(1 - yp))
#' @note eps is used where log loss is undefined for p=0 or p=1, so probabilities are clipped to: max(eps, min(1 - eps, p))
#' @note If likelihood is output, the eps and local arguments are ignored.  
#'       
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references C.M. Bishop (2006). Pattern Recognition and Machine Learning. Springer, p. 209.
#'
#' @examples 
#'
#'   require(randomForest)
#'     data(iris)
#'     iris$Species <- ifelse( iris$Species == "versicolor", 1, 0 ) 
#'     # Add some noise
#'       idx1 <- which(iris$Species %in% 1)
#'       idx0 <- which( iris$Species %in% 0)
#'       iris$Species[sample(idx1, 2)] <- 0
#'       iris$Species[sample(idx0, 2)] <- 1
#'     
#'  ( mdl <- randomForest(x=iris[,1:4], y=as.factor(iris[,"Species"])) )
#' 	
#'   # Global log loss	
#'     logLoss(y = iris$Species, p = predict(mdl, iris[,1:4], type="prob")[,2]) 
#' 			   
#'   # Local log loss
#'     ( ll <- logLoss(y = iris$Species, p = predict(mdl, iris[,1:4], 
#'                    type="prob")[,2], global = FALSE) )
#' 
#'   # Log likelihood loss
#'     logLoss(y = iris$Species, p = predict(mdl, iris[,1:4], 
#' 			    type="prob")[,2], likelihood = TRUE) 
#' 				   
#' @export 
logLoss <- function(y, p, likelihood = FALSE, global = TRUE, eps = 1e-15) {
  if( length(p) != length(y) ) stop("Error: y and p are not the same length")
    if( min(p) < 0 | max(p) > 1) stop("p is out of {0-1} probability range")
      if( !is.numeric(y) ) y <- as.numeric(as.character(y))
        if( min(y) < 0 | max(y) > 1) stop("y is not binomial {0,1}")
    loglike.loss <- function(y, x) {
      if(class(y) == "factor" | class(y) == "character")
        y <- as.numeric(as.character(y))  
      x <- 1/(1 + exp(-x))
      grad <- x - y
      hess <- x * (y - x)
      return(hess)
    }
  if( likelihood == FALSE) {	
    if( global == TRUE) { 
      p <- matrix(sapply(p, function(x) max(eps, x)), nrow = length(p))
        p <- matrix(sapply(p, function(x) min(1 - eps, x)), nrow = length(p))
        ll <- sum(y * log(p) + (1 - y) * log(1 - p)) 
      ll <- ll * -1 / length(y)
	  class(ll) <- "global.log.loss"
    } else {
    ll <- vector()  
      for( i in 1:length(p)) {  
        ll <- append(ll, pmin(pmax(p[i], eps), 1 - eps) - (sum(y[i] * log(p[i]) + 
	                (1 - y[i]) * log(1 - p[i]))) / length(y[i]) )
      } 
          ll[is.nan(ll)] <- 0
	    ll <- list(y = y, p = p, log.loss = ll)
      class(ll) <- "local.log.loss"	
    } 
  } else {
    ll <- loglike.loss(y, p)
  }
  return(ll)  
}
