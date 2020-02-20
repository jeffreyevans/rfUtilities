#' @title Random Forests class-level sensitivity analysis 
#' @description Performs a sensitivity analysis on a specified class 
#'              in a random forests model 
#' 
#' @param x         randomForest or ranger class object
#' @param class     Which class to perturb
#' @param p         Proportion of class to be randomized
#' @param nperm     Number of permutations
#' @param plot      Plot results (TRUE/FALSE)
#' @param seed      Random seed value 
#'
#' @return List object with following components:
#' \itemize{ 
#' \item  mean.error - [Mean of RMSE]
#' \item  sd.error - [Standard deviation of RMSE]
#' \item  rmse - [Root mean squared error (RMSE) for each perturbed probability]
#' \item  probs - [data.frame with "true" estimate in first column and 
#'                 perturbed probabilities in subsequent columns.]        
#' }
#'
#' @details
#' Wildlife survey data likely decreases the proportion of imperfect detection 
#' (false absences or presences) but can still be a source of error. Because of 
#' this it is often necessary to test the model sensitivity of a given class 
#' (eg., used verses available habitat). Model sensitivity of false absences is 
#' evaluated by randomly assigning a proportion of the specified positive class 
#' to the other, refitting the model and estimating the probabilities. Each 
#' perturbed estimate is compared against the "true" estimate. Currently only 
#' supports binomial models.
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @references
#' Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species 
#'   distribution and change using Random Forests CH.8 in Predictive Modeling 
#'   in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer
#' @references
#' Gardner, R.H., R.V. O'Neill, M.G. Turner, and V.H. Dale (1989). Quantifying 
#'   scale-dependent effects of animal movements with simple percolation models. 
#'   Landscape Ecology 3:217-227.
#'
#' @examples
#' library(randomForest)
#' library(ranger) 
#' data(iris)
#'   y <- as.factor(ifelse(iris$Species == "setosa" | 
#'                  iris$Species == "virginica", 1, 0) )
#'     xdata <- iris[,1:4] 
#' 
#' rf.mdl <- randomForest(xdata, y, ntree=501) 
#'   ua <- rf.class.sensitivity(rf.mdl, nperm=20) 
#'   
#' rf.mdl <- randomForest(xdata, y, probability=TRUE) 
#'   ua <- rf.class.sensitivity(rf.mdl, nperm=20) 
#'            
#' @export        
rf.class.sensitivity <- function(x, class = "1", p = 0.05, nperm = 999, 
                                 plot=TRUE, seed=NULL) {	 
 if (!inherits(x, c("randomForest", "ranger"))) 
    stop("x is not randomForest class object")
  if(length(grep("~", x$call[[2]])) > 0)
    stop("This package does not support a formula interface, please use x, y arguments")
  if(inherits(x, "randomForest")) {
    if (x$type != "classification") 
      stop("Unsupervised classification or regression not supported")
  } else if(inherits(x, "ranger")) {
    if(tolower(x$treetype) != "probability estimation")	
	  stop("Not a probability forests, please use;  probability = TRUE  \n")
  }
  if(!is.null(seed)) { set.seed(seed) } else { set.seed(.Random.seed[1]) }
  
  # formating call and pulling data
  a <- as.list(x$call)[-1] 
    xdata = eval(a[[which(names(a) == "x")]]) 
    ydata = eval(a[[which(names(a) == "y")]]) 
 
    rmse <- function(o,p) sqrt( mean( (o - p )^2 ) )	
	values <- as.character(unique(y))	
	if(inherits(x, "randomForest")) {
	  pred <- stats::predict(x, xdata, type="prob")
    } else if(inherits(x, "ranger")) {
	  pred <- stats::predict(x, data = xdata)$predictions
	}	
	nc <- which(colnames(pred) == class)
	  mpred <- data.frame(obs = pred[,nc])
	    names(mpred) <- "obs" 
    sample.size = round( (length(ydata) * p) / length(values), digits=0)
	if(plot == TRUE) 
	  graphics::plot(stats::density(mpred[,1]), main="Perturbed model probabilities", type="n")
	    for(i in 1:nperm) {
		  y <- ydata
            samp <- sample(which( y %in% class ), sample.size)
	          y[samp] <- values[values != class] 
		        a[["y"]] <- y
	        if(inherits(x, "randomForest")) {
	          pmdl <- do.call(randomForest::randomForest, a)			
	    	  mpred <- data.frame(mpred, stats::predict(pmdl, xdata, type="prob")[,nc])
            } else if(inherits(x, "ranger")) {
			  pmdl <- do.call(ranger::ranger, a)
	          mpred <- data.frame(mpred, stats::predict(x, data = xdata)$predictions[,nc])
	        }	
	      if(plot == TRUE) {
		    pden <- stats::density(mpred[,i+1])
		    graphics::lines(pden$x, pden$y, col="grey")
		  }
	    }
	    if(plot == TRUE) {
          d <- stats::density(mpred[,1])
	      graphics::lines(d$x, d$y, col="red", lwd=2)
	      graphics::legend("topleft", legend=c("observed", "perturbed"), 
	                       col=c("red", "grey"), lwd=c(2,1), bg="white")
        }		
          names(mpred)[2:ncol(mpred)] <- paste0("sim", seq(1:(ncol(mpred)-1))) 
    	    error <- vector()
        for(i in 2:ncol(mpred)) error <- append(error, rmse(mpred[,1],mpred[,i]))
      cat("Mean error: ", mean(error), "\n")	
    cat("Standard deviation of error: ", stats::sd(error), "\n")	 
  list(mean.error = mean(error), sd.error = stats::sd(error), rmse = error, probs = mpred)
  }
