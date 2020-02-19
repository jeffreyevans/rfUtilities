#' @title Random Forest model significance test 
#' @description Performs significance test for classification and regression 
#'              Random Forests models. 
#'
#' @param x              randomForest class object
#' @param p              p-value to test for significance in regression models
#' @param q              Quantile threshold to test classification models
#' @param nperm          Number of permutations
#' @param randomization  Fraction (0.01-1) of randomization, 
#'                       default is 1 (total randomization) 
#'
#' @return A list class object with the following components:
#' For Regression problems:    
#' * RandR.square Vector of random R-square values 
#' * R.square The R-square of the "true" model
#' * p.value p-values of randomizations of R-square
#' * ks.p.value p-value(s) evaluation of Kolmogorov-Smirnov test
#' * nPerm number of permutations 
#' * rf.type Type of Random Forests
#' * rand.frac Amortization fraction
#'
#' For Classification problems:
#' * RandOOB Vector of random out-of-bag (OOB) values 
#' * RandMaxError Maximum error of randomizations 
#' * test.OOB Error OOB error of the "true" model
#' * test.MaxError maximum class OOB error of the "true" model
#' * p.value p-value based on Mcnemar's test
#' * oop.p.value p-value based on permutation of OOB error
#' * nPerm Number of permutations 
#' * rf.type Type of Random Forests
#' * rand.frac Amortization fraction
#' @md
#'
#' @notes
#' Please note that previous versions of this function required xdata and "..."  
#' arguments that are no longer necessary. The model object is now used in 
#' obtaining the data and arguments used in the original model  
#'
#' @details
#' If the p-value is small, it suggests a near certainty 
#' that the difference between the two populations is 
#' significant. alternative = c("two.sided", "less", "greater")
#'
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo 
#'               boreas connectivity in Yellowstone National Park with landscape 
#'               genetics. Ecology 91:252-261
#' @references Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling 
#'               species distribution and change using Random Forests CH.8 in  
#'               PredictiveModeling in Landscape Ecology eds Drew, CA, 
#'               Huettmann F, Wiersma Y. Springer 
#' 
#' @examples 
#' \dontrun{
#' #### Regression
#' library(randomForest)
#' library(ranger)
#' 
#'   set.seed(1234)	
#'   data(airquality)
#'   airquality <- na.omit(airquality)
#'  
#'  # randomForest 
#'  ( rf.mdl <- randomForest(x=airquality[,2:6], y=airquality[,1]) )
#'    ( rf.perm <- rf.significance(rf.mdl, nperm=99) )
#'
#'  # ranger
#'  ( rf.mdl <- ranger(x=airquality[,2:6], y=airquality[,1]) )
#'    ( rf.perm <- rf.significance(rf.mdl, nperm=99) )
#'
#'  
#' #### Classification
#' ydata = as.factor(ifelse(airquality[,1] < 40, 0, 1))
#' ( rf.mdl <- ranger(x = airquality[,2:6], y = ydata) )
#'       ( rf.perm <- rf.significance(rf.mdl, nperm=99) )
#'	   
#' ( rf.mdl <- randomForest(x = airquality[,2:6], y = ydata) )
#'      ( rf.perm <- rf.significance(rf.mdl, nperm=99) )
#' 
#' 
#'   set.seed(1234)	
#'     data(iris)
#'       iris$Species <- as.factor(iris$Species) 
#'	
#'  	
#'  ( rf.mdl <- randomForest(x=iris[,1:4], y=iris[,"Species"]) )
#'    ( rf.perm <- rf.significance(rf.mdl, nperm=99) )
#'
#'  ( rf.mdl <- ranger(x=iris[,1:4], y=iris[,"Species"]) )
#'    ( rf.perm <- rf.significance(rf.mdl, nperm=99) )
#'
#' }
#'
#' @exportClass significance 
#' @export
rf.significance <- function(x, nperm = 999, randomization = 1, 
                            kappa = FALSE) { 
  if (!inherits(x, c("randomForest", "ranger"))) 
    stop("x is not randomForest class object")
  if(length(grep("~", x$call[[2]])) > 0)
      stop("This package does not support a formula interface, please use x, y arguments")
  if(randomization <= 0.01 | randomization > 1)
    stop("randomization fraction must be between 0.01 - 1")
  if(inherits(x, "randomForest")) {
    mtype <- tolower(x$type)
  } else if(inherits(x, "ranger")) {
    mtype <- tolower(x$treetype)
  }
    calculate.p <- function(x, test, nperm) { 
      if ( length( x[x >= test] ) < 1 ) { 
	    return( 0 )
	  } else { 
        return( round(length(x[x >= test]) / nperm, 5) )
	  }	
    }
  options(warn=-1)	
  #*********************************	
  if(mtype == "classification") {
  #*********************************
    ptest = TRUE
	mdl.call <- as.list(x$call)[-1]  
	  obs = eval(mdl.call[[which(names(mdl.call) == "y")]])
	    if(length(unique(obs)) > 2) ptest = FALSE
	if(inherits(x, "ranger")) { 
	  cm <- accuracy(x$confusion.matrix)
	  if(ptest == FALSE) {
	    test.p <- NA
	  } else {
        test.p <- mcnemar.test(table(x$predictions, obs))$p.value
      }	  
	} else if(inherits(x,"randomForest")) {
	  cm <- accuracy(x$confusion[,1:(ncol(x$confusion)-1)])
	  if(ptest == FALSE) {
	    test.p <- NA
	  } else {
        test.p <- mcnemar.test(table(x$predicted, obs))$p.value	  
	  }
	}
  if(kappa) { test.oob <- cm$kappa }else{ test.oob <- cm$PCC }	
    test.max <- max(cm$producers.accuracy, na.rm=TRUE) 
	rand.oob <- vector()
	rand.max <- vector()
    rand.p <- vector()	
      for(i in 1:nperm) {
        a <- as.list(x$call)[-1]  
        ydata = eval(a[[which(names(a) == "y")]])
		  if(kappa == TRUE) {
            if(length(unique(ydata)) > 2) {
			  kappa = FALSE 
			  kappa.warn = TRUE
			}
          } 
          ns <- round(length(ydata) * randomization,0)
            yidx <- sample(1:length(ydata), ns) 
			  ridx <- sample(yidx,replace=FALSE)
			  y.new <- ydata 
		    for(r in 1:length(yidx)) { y.new[yidx] <- y.new[ridx]  }
		  a[["y"]] <- y.new 	
		if(inherits(x, "ranger")) {    
          rfr <- do.call(ranger::ranger, a)
		    cm <- accuracy(rfr$confusion.matrix)
			if(ptest == FALSE) {
	          test.p <- NA
	        } else {
			  rand.p[i] <- mcnemar.test(table(rfr$predictions, ydata))$p.value
			}  
        } else if(inherits(x,"randomForest")) {
          rfr <- do.call(randomForest::randomForest, a)
		    cm <- accuracy(rfr$confusion[,1:(ncol(rfr$confusion)-1)])
	        if(ptest == FALSE) {
	          test.p <- NA
	        } else {			
              rand.p[i] <- mcnemar.test(table(rfr$predicted, ydata))$p.value
			}
        }
	    if(kappa) { rand.oob[i] <- cm$kappa }else{ rand.oob[i] <- cm$PCC }
	    rand.max[i] <- max(cm$producers.accuracy, na.rm=TRUE)		
      }
    if(exists("kappa.warn")){
      if(kappa.warn == TRUE)
        warning("Cannot run Kappa on more than 2 classes, reverting to PCC")
    }   
	if(ptest == FALSE) {
	  p.value <- NA
	} else {
	  p.value <- class.p(x=rand.p, test=test.p, nperm=nperm)	 
	}
	oob.p.value = class.p(x=rand.oob, test=test.oob, nperm=nperm)
    sig <- list(RandOOB = 1-(rand.oob/100), RandMaxError = 1-(rand.max/100), 
                test.OOB = 1-(test.oob/100), test.MaxError = 1-(test.max/100), 
                p.value = p.value, oob.p.value = oob.p.value, nPerm = nperm, 
				rf.type = mtype, rand.frac = randomization, use.kappa=kappa)		   
    #*********************************	
    } else if(mtype == "regression") {
    #*********************************	
	mdl.call <- as.list(x$call)[-1]  
	obs = eval(mdl.call[[which(names(mdl.call) == "y")]])
	  if(inherits(x, "ranger")) { 
	    test.rsq <- x$r.squared
	    test.mse <- x$prediction.error
		p.equal.test <- stats::ks.test(x$predictions, stats::ecdf(obs), 
	                  alternative ="two.sided", exact=TRUE)$p.value
	    p.lt.test <- stats::ks.test(x$predictions, stats::ecdf(obs), 
	                  alternative ="less", exact=TRUE)$p.value
        p.gt.test <- stats::ks.test(x$predictions, stats::ecdf(obs), 
	                  alternative ="greater", exact=TRUE)$p.value
	  } else if(inherits(x,"randomForest")) {
	    test.rsq <- stats::median(x$rsq)
	    test.mse <- stats::median(x$mse)
		p.equal.test <- stats::ks.test(x$predicted, stats::ecdf(obs), 
	                  alternative ="two.sided", exact=TRUE)$p.value
	    p.lt.test <- stats::ks.test(x$predicted, stats::ecdf(obs), 
	                  alternative ="less", exact=TRUE)$p.value
        p.gt.test <- stats::ks.test(x$predicted, stats::ecdf(obs), 
	                  alternative ="greater", exact=TRUE)$p.value
	  }
	  rand.rsq <- vector()
	  rand.mse <- vector()
      p.equal <- vector()	
      p.lt <- vector() 
      p.gt <- vector() 
      for(i in 1:nperm) {
        a <- as.list(x$call)[-1]  
        ydata = eval(a[[which(names(a) == "y")]])
          ns <- round(length(ydata) * randomization,0)
            yidx <- sample(1:length(ydata), ns) 
			ridx <- sample(yidx,replace=FALSE)
		  y.new <- ydata 
			for(r in 1:length(yidx)) { y.new[yidx] <- y.new[ridx]  }
			a[["y"]] <- y.new 
		if(inherits(x, "ranger")) {    
          rfr <- do.call(ranger::ranger, a)
		    rand.rsq[i] <- rfr$r.squared  
		    rand.mse[i] <- rfr$prediction.error
			p.equal[i] <- stats::ks.test(rfr$predictions, stats::ecdf(obs), 
	                          alternative ="two.sided", exact=TRUE)$p.value
	        p.lt[i] <- stats::ks.test(rfr$predictions, stats::ecdf(obs), 
	                      alternative ="less", exact=TRUE)$p.value
            p.gt[i] <- stats::ks.test(rfr$predictions, stats::ecdf(obs), 
	                      alternative ="greater", exact=TRUE)$p.value
        } else if(inherits(x,"randomForest")) {
          rfr <- do.call(randomForest::randomForest, a)
		    rand.rsq[i] <- stats::median(rfr$rsq)
		    rand.mse[i] <- stats::median(rfr$mse)
			p.equal[i] <- stats::ks.test(rfr$predicted, stats::ecdf(obs), 
	                         alternative ="two.sided", exact=TRUE)$p.value
	        p.lt[i] <- stats::ks.test(rfr$predicted, stats::ecdf(obs), 
	                      alternative ="less", exact=TRUE)$p.value
            p.gt[i] <- stats::ks.test(rfr$predicted, stats::ecdf(obs), 
	                      alternative ="greater", exact=TRUE)$p.value
        }
      }
    p.equal[is.nan(p.equal)] <- 0	
	  p.lt[is.nan(p.lt)] <- 0  
	    p.gt[is.nan(p.gt)] <- 0  
	ks.pValue = data.frame(
	     ks.equal = regress.p(x = p.equal, test = p.equal.test, nperm = nperm),
	     ks.lt = regress.p(x = p.lt, test = p.lt.test, nperm = nperm),
	     ks.gt = regress.p(x = p.gt, test = p.gt.test, nperm = nperm))
	rsq.pValue = round(regress.p(x = rand.rsq, test = test.rsq, nperm = nperm), digits=5)
    sig <- list(RandR.square = rand.rsq, R.square = test.rsq, p.value = rsq.pValue,  
                ks.p.value = ks.pValue, nPerm = nperm, rf.type = mtype,  
				rand.frac = randomization)
    }
	  options(warn=0)
    class(sig) <- c("significance","list")
  return( sig )			   
} 
