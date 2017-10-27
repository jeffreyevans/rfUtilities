#' @title Isotonic probability calibration
#' @description Performs an isotonic regression calibration of posterior probability to minimize log loss.
#'    
#' @param y                Binominal response varible used to fit model   
#' @param p                Estimated probabilities from fit model
#' @param regularization   (FALSE/TRUE) should regularization be performed on the probabilities? (see notes) 
#'
#' @return a vector of calibrated probabilities 
#' 
#' @note Isotonic calibration can correct for monotonic distortions. 
#' @note regularization defines new minimum and maximum bound for the probabilities using:
#' @note   pmax = ( n1 + 1) / (n1 + 2), pmin = 1 / ( n0 + 2); where n1 = number of prevalence values and n0 = number of null values 
#'       
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references Platt, J. (1999) Probabilistic outputs for support vector machines and comparison to regularized likelihood methods. Advances in Large Margin Classifiers (pp 61-74).
#' @references Niculescu-Mizil, A., & R. Caruana (2005) Obtaining calibrated probabilities from boosting. Proc. 21th Conference on Uncertainty in Artificial Intelligence (UAI 2005). AUAI Press.
#'
#' @examples 
#'  library(randomForest)
#'    data(iris)
#'    iris$Species <- ifelse( iris$Species == "versicolor", 1, 0 ) 
#'    
#'    # Add some noise
#'    idx1 <- which(iris$Species %in% 1)
#'    idx0 <- which( iris$Species %in% 0)
#'    iris$Species[sample(idx1, 2)] <- 0
#'    iris$Species[sample(idx0, 2)] <- 1
#'  
#'  # Specify model  
#'  y = iris[,"Species"] 
#'  x = iris[,1:4]
#'  set.seed(4364)  
#'  ( rf.mdl <- randomForest(x=x, y=factor(y)) )
#'  y.hat <- predict(rf.mdl, iris[,1:4], type="prob")[,2] 
#'  
#'  # Calibrate probabilities
#'  calibrated.y.hat <- probability.calibration(y, y.hat, regularization = TRUE) 
#'
#'  # Plot calibrated aganist origianl probability estimate
#'  plot(density(y.hat), col="red", xlim=c(0,1), ylab="Density", xlab="probabilities",
#'       main="Calibrated probabilities" )
#'         lines(density(calibrated.y.hat), col="blue")
#'           legend("topright", legend=c("original","calibrated"), 
#'  	            lty = c(1,1), col=c("red","blue"))
#'   
#' @export 
probability.calibration <- function(y, p, regularization = FALSE) {
  if( length(p) != length(y)) stop("Vectors do not match")
  if(!is.numeric(y)) if(is.factor(y)) { y <- as.numeric(as.character(y)) } else { 
     stop("y is not valid binomial vector") }
  if(length(unique(y)) > 2) stop("y is not a valid binomial vector") 
  if(!min(unique(y)) == 0) stop("y is not a valid binomial vector")
  if(!max(unique(y)) == 1) stop("y is not a valid binomial vector")  
  if(!is.numeric(p)) stop("p arguments must be numeric")    
    if(regularization == TRUE) {
      p.max <- (length(y[y == 1]) + 1) / (length(y[y == 1]) + 2)
      p.min <- 1 / (length(y[y == 0]) + 2)
	  p <- ifelse( p < p.min, p.min, p)
	  p <- ifelse( p > p.max, p.max, p)
    }
  idx <- duplicated(p)
    idx <- which( idx == TRUE)
      p.unique <- p[-idx]
      y.unique <- y[-idx]  
  isotonic.calibration <- function(iso, x0) {
    o = iso$o
      if (is.null(o))
        o = 1:length(x)
        x = iso$x[o]
        y = iso$yf
      ind = cut(x0, breaks = x, labels = FALSE, include.lowest = TRUE)
      min.x <- min(x)
      max.x <- max(x)
      adjusted.knots <- iso$iKnots[c(1, which(iso$yf[iso$iKnots] > 0))]
      fits = sapply(seq(along = x0), function(i) {
      j = ind[i]
    if (is.na(j)) {
      if (x0[i] > max.x) j <- length(x)
        else if (x0[i] < min.x) j <- 1
      }
    upper.step.n <- min(which(adjusted.knots > j))
    upper.step <- adjusted.knots[upper.step.n]
    lower.step <- ifelse(upper.step.n==1, 1, adjusted.knots[upper.step.n -1] )
    denom <- x[upper.step] - x[lower.step]
    denom <- ifelse(denom == 0, 1, denom)
    val <- y[lower.step] + (y[upper.step] - y[lower.step]) * (x0[i] - x[lower.step]) / (denom)
    val <- ifelse(val > 1, max.x, val)
      val <- ifelse(val < 0, min.x, val)
        val <- ifelse(is.na(val), max.x, val)
        val
      })
      return( fits )
    } 
    iso.mdl <- stats::isoreg(p.unique, y.unique)
  return( isotonic.calibration(iso.mdl, p) )
}	
