#' @title Random Forest Classification or Regression Model Cross-validation 
#' @description Implements a permutation test cross-validation for Random Forests models
#'    
#' @param x                 random forest object
#' @param xdata             x data used in model
#' @param ydata             optional y data used in model, default is to use x$y from model object
#' @param p                 Proportion data withhold (default p=0.10)
#' @param n                 Number of cross validations (default n=99)
#' @param seed              Sets random seed in R global environment
#' @param normalize         (FALSE/TRUE) For regression, should rmse, mbe and mae be normalized using (max(y) - min(y))
#' @param bootstrap         (FALSE/TRUE) Should a bootstrap sampling be applied. If FALSE, an n-th percent withold will be conducted
#' @param trace             Print iterations
#' @param ...               Additional arguments passed to Random Forests 
#'
#' @return  For classification a "rf.cv"", "classification" class object with the following components:
#' \itemize{ 
#' \item  cross.validation$cv.users.accuracy        Class-level users accuracy for the subset cross validation data   
#' \item  cross.validation$cv.producers.accuracy    Class-level producers accuracy for the subset cross validation data    
#' \item  cross.validation$cv.oob                   Global and class-level OOB error for the subset cross validation data    
#' \item  model$model.users.accuracy                Class-level users accuracy for the model 
#' \item  model$model.producers.accuracy            Class-level producers accuracy for the model 
#' \item  model$model.oob                           Global and class-level OOB error for the model  
#'  }
#'
#' @return  For regression a "rf.cv", "regression" class object with the following components:
#' \itemize{ 
#' \item  fit.var.exp    Percent variance explained from specified fit model 
#' \item  fit.mse        Mean Squared Error from specified fit model   
#' \item  y.rmse         Root Mean Squared Error (observed vs. predicted) from each Bootstrap iteration (cross-validation)    
#' \item  y.mbe          Mean Bias Error from each Bootstrapped model
#' \item  y.mae          Mean Absolute Error from each Bootstrapped model
#' \item  D              Test statistic from Kolmogorov-Smirnov distribution Test (y and estimate)
#' \item  p.val          p-value for Kolmogorov-Smirnov distribution Test (y and estimate)
#' \item  model.mse      Mean Squared Error from each Bootstrapped model
#' \item  model.varExp   Percent variance explained from each Bootstrapped model   
#'  }
#'
#' @details
#' For classification problems, the cross-validation statistics are based on the prediction error on the withheld data: 
#' Total observed accuracy represents the percent correctly classified (aka, PCC) and is considered as a naive measure of agreement. 
#' The diagonal of the confusion matrix represents correctly classified observations where off-diagonals represent cross-classification error. The primary issue with this evaluation is that does not reveal if error was evenly distributed between classes.
#' To represent the balance of error one can use omission and commission statistics such as estimates of users and producers accuracy. User's accuracy corresponds to error of commission (inclusion), observations being erroneously included in a given class.
#' The commission errors are represented by row sums of the matrix. Producer's accuracy corresponds to error of omission (exclusion), observations being erroneously excluded from a given class. The omission errors are represented by column sums of the matrix.
#' None of the previous statistics account for random agreement influencing the accuracy measure. The kappa statistic is a chance corrected metric that reflects the difference between observed agreement and agreement expected by random chance.
#' A kappa of k=0.85 would indicate that there is 85% better agreement than by chance alone.  
#'   \itemize{ 
#'   \item   pcc = [Number of correct observations / total number of observations] 
#'   \item   pcc = [Number of correct observations / total number of observations] 
#'   \item   producers accuracy =  [Number of correct / total number of correct and omission errors] 
#'   \item   k = (observed accuracy - chance agreement) / (1 - chance agreement) where; change agreement = sum[product of row and column totals for each class] 
#'    }
#' For regression problems, a Bootstrap is constructed and the subset models MSE and percent variance explained is reported. Additional, the RMSE between the withheld response variable (y) and the predicted subset model   
#'
#' @author Jeffrey S. Evans <jeffrey_evans<at>tnc.org>
#'
#' @references Evans, J.S. and S.A. Cushman (2009) Gradient Modeling of Conifer Species Using Random Forest. Landscape Ecology 5:673-683.
#' @references Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas connectivity in Yellowstone National Park with landscape genetics. Ecology 91:252-261
#' @references Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer 
#' 
#' @examples 
#' \dontrun{
#' library(randomForest)
#'
#' # For classification
#'   data(iris)
#'     iris$Species <- as.factor(iris$Species)    	
#'       set.seed(1234)	
#' ( rf.mdl <- randomForest(iris[,1:4], iris[,"Species"], ntree=501) )
#'   ( rf.cv <- rf.crossValidation(rf.mdl, iris[,1:4], p=0.10, n=99, ntree=501) )
#'
#'    # Plot cross validation versus model producers accuracy
#'    par(mfrow=c(1,2)) 
#'      plot(rf.cv, type = "cv", main = "CV producers accuracy")
#'      plot(rf.cv, type = "model", main = "Model producers accuracy")
#'
#'    # Plot cross validation versus model oob
#'    par(mfrow=c(1,2)) 
#'      plot(rf.cv, type = "cv", stat = "oob", main = "CV oob error")
#'      plot(rf.cv, type = "model", stat = "oob", main = "Model oob error")	  
#'
#' # For regression
#' data(airquality)
#' airquality <- na.omit(airquality) 
#' rf.mdl <- randomForest(y=airquality[,"Ozone"], x=airquality[,2:4])
#' ( rf.cv <- rf.crossValidation(rf.mdl, airquality[,2:4], 
#'                               p=0.10, n=99, ntree=501) )
#'  par(mfrow=c(2,2))
#'    plot(rf.cv)  
#'    plot(rf.cv, stat = "mse")
#'    plot(rf.cv, stat = "var.exp")
#'	plot(rf.cv, stat = "mae")
#' }	 
#'	  
#' @seealso \code{\link[randomForest]{randomForest}} for randomForest ... options 
#'
#' @exportClass rf.cv
#' @export 	
rf.crossValidation <- function(x, xdata, ydata=NULL, p=0.10, n=99, seed=NULL, normalize = FALSE, 
                               bootstrap = FALSE, trace = FALSE, ...) {
  if(!any(class(x) %in% c("randomForest","list"))) stop("x is not a randomForest object")
  # add check length of y and dim of x
    if(!is.null(seed)) { set.seed(seed) }
	  if(bootstrap) cat("Bootstrap sampling is being applied,", paste("p",p,sep="="), "argument is ignored", "\n")
    
  if (x$type == "unsupervised") { 
    stop("Unsupervised classification not supported")

  ###########################################
  #### Start regression cross-validation ####	   
  } else if (x$type == "regression") {
    cat("running:", x$type, "cross-validation", "with", n, "iterations", "\n")
	# Validation statistics (RMSE, MBE, MAE)
        # Root Mean Square Error (RMSE) 
        rmse <- function(y, x, norm = FALSE){
          if( length(y[is.na(y)]) > 0) stop("NA values present in y data")
          if( length(x[is.na(x)]) > 0) stop("NA values present in x data")
            e <- sqrt(mean((y - x)^2)) 
        	  if( norm ) e <- e / diff(range(y, na.rm = TRUE))
            return( e )		   
        }          
        # Mean Bias Error (MBE) 
        mbe <- function(y, x, norm = FALSE){
          if( length(y[is.na(y)]) > 0) stop("NA values present in y data")
          if( length(x[is.na(x)]) > 0) stop("NA values present in x data")
            e <- mean(x - y)
            # e <- mean(x - y) / mean(y) * 100 
        	  if( norm ) e <- e / diff(range(y, na.rm = TRUE))
            return( e )		   
        }     
        # Mean Absolute Error (MAE) 
        mae <- function(y, x, norm = FALSE){
          if( length(y[is.na(y)]) > 0) stop("NA values present in y data")
          if( length(x[is.na(x)]) > 0) stop("NA values present in x data") 
            e <- mean(abs(y - x))
              if( norm ) e <- e / diff(range(y))
            return( e )		   
        }
	    # Kolmogorov-Smirnov Test (1=D, 2=p.value)
        ks <- function(y, x, s = c(1,2)) {
		  stats::ks.test(x, stats::ecdf(y))[s[1]]
		}
    # Define validation vectors
    y.rmse <- rep(NA, n)  
 	y.mae <- rep(NA, n)
	y.mbe <- rep(NA, n)
	model.varExp <- rep(NA, n)
	model.mse <- rep(NA, n)
	ks.p <- rep(NA, n)
	ks.d <- rep(NA, n)
	if( is.null(ydata) ) { ydata <- x$y }
    if(bootstrap) boot.sample.size <- rep(NA, n)	
      sample.size = round( (length(ydata) * p), digits=0) #population sample size 
	  for(i in 1:n) {
	    if(trace) cat("running iteration:", i, "\n")
        dat <- data.frame(y=ydata, xdata)
        # Draw random sample		
	    if(!bootstrap) {
	      sidx <- sample(1:nrow(dat), sample.size)   
	        dat.sub <- dat[-sidx,]                   
              dat.cv <- dat[sidx,]		            
	    } else {	
	      dat.sub <- dat[sample(1:nrow(dat), replace=TRUE),]              
          dat.cv <- dat[which(!rownames(dat) %in% rownames(dat.sub)),]	   
        }	
	     rf.fit <- randomForest::randomForest(y=dat.sub[,"y"], 
		                    x=dat.sub[,2:ncol(dat.sub)], ...)
 		   model.mse[i] <- rf.fit$mse[length(rf.fit$mse)]
	       model.varExp[i] <- round(100*rf.fit$rsq[length(rf.fit$rsq)], digits=2)          
		   y.rmse[i] <- rmse(dat.cv[,"y"], stats::predict(rf.fit, 
		                    newdata = dat.cv[,2:ncol(dat.cv)]), 
							norm=normalize) 
           y.mbe[i] <- mbe(dat.cv[,"y"], stats::predict(rf.fit, 
		                   newdata = dat.cv[,2:ncol(dat.cv)]), 
						   norm=normalize) 
		   y.mae[i] <- mae(dat.cv[,"y"], stats::predict(rf.fit, 
		                   newdata = dat.cv[,2:ncol(dat.cv)]), 
						   norm=normalize) 
           ks.p[i] <- as.numeric(ks(dat.cv[,"y"], stats::predict(rf.fit, 
		                 newdata = dat.cv[,2:ncol(dat.cv)]), s=2))  
           ks.d[i] <- as.numeric(ks(dat.cv[,"y"], stats::predict(rf.fit, 
		                 newdata = dat.cv[,2:ncol(dat.cv)]), s=1))
		 if(bootstrap) boot.sample.size[i] <- nrow(dat.cv) 
	  }		 
	 r.cv <- list(fit.var.exp=round(100*x$rsq[length(x$rsq)], digits=2), 
	              fit.mse=stats::median(x$mse), y.rmse = y.rmse, 
				  y.mbe = y.mbe, y.mae = y.mae, D = ks.d, 
				  p.val = ks.p, model.mse = model.mse, 
				  model.varExp = model.varExp )			
      class(r.cv) <- c("rf.cv", "regression", "list")
    return( r.cv )
	 
   ###############################################
   #### Start classification cross-validation ####
   } else if (x$type == "classification") {
      cat("running:", x$type, "cross-validation", "with", n, "iterations", "\n")
	    if( is.null(ydata) ) { ydata <- x$y }
    if(bootstrap) boot.sample.size <- vector()	
    classes <- as.vector(levels( ydata ))	
      sample.size = round( (length(ydata) * p) / length(x$classes), digits=0) 
	    cv.ua <- as.data.frame(array(0, dim=c(0,length(classes))))
	      cv.pa <- as.data.frame(array(0, dim=c(0,length(classes))))
	      mdl.ua <- as.data.frame(array(0, dim=c(0,length(classes)))) 
	    mdl.pa <- as.data.frame(array(0, dim=c(0,length(classes))))
      mdl.oob <- as.data.frame(array(0, dim=c(0,2+length(classes))))
	cv.oob <- as.data.frame(array(0, dim=c(0,2+length(classes))))	
	
    # Create class-level sample sizes
	nclass <- length(unique(ydata)) 
	p.class = p / nclass
	sample.sizes <- vector()
	  for(s in 1:nclass) { 
	    sample.sizes[s] <- round(length(ydata[ydata == unique(ydata)[s]]) * 
		                         p.class, digits=0)    
	  } 
      for(i in 1:n) {
	    if(trace) cat("running iteration:", i, "\n")
        # Draw random sample		
	    if(!bootstrap) {
          sidx <- list()				
            for(s in 1:nclass) {
              sidx[[s]] <- sample(which(ydata %in% unique(ydata)[s]), sample.sizes[s])
            }
			sidx <- unlist(sidx)  
              tx <- xdata[sidx,]
              ty <- ydata[sidx]
              mx <- xdata[-sidx,]
              my <- ydata[-sidx]	           
	    } else {	
          dat <- data.frame(y=ydata, xdata) 	
            tx <- dat[sample(1:nrow(dat), replace=TRUE),]
              ty <- tx$y 
                tx <- tx[,2:ncol(tx)]
              mx <- dat[-which(!rownames(dat) %in% rownames(tx)),]
            my <- mx$y
          mx <- mx[,2:ncol(mx)] 
       }
      rf.fit <- randomForest::randomForest(y=as.factor(my), x = mx, ytest=as.factor(ty), 
	                                       xtest=tx, ...)        
      options(warn=-1)
      cv.acc <- accuracy(rf.fit$test$predicted,  ty)
	    if(bootstrap) boot.sample.size <- append(boot.sample.size, length(my))
 		cv.ua <- rbind(cv.ua, cv.acc$users.accuracy) 
		  cv.pa <- rbind(cv.pa, cv.acc$producers.accuracy) 
	        cv.oob <- rbind(mdl.oob, c(apply(rf.fit$test$err.rate, MARGIN = 2, stats::median),
		                    cv.acc$kappa))	  
		mdl.acc <- accuracy(rf.fit$predicted, my) 
		  if(!length(classes) == length(unique(my))) {
	      #### add code to check for presence of classes and set to NA if absent classes occur   
		  }
	    mdl.ua <- rbind(mdl.ua, mdl.acc$users.accuracy) 
		  mdl.pa <- rbind(mdl.pa, mdl.acc$producers.accuracy) 
            mdl.oob <- rbind(mdl.oob, c(apply(rf.fit$err.rate, MARGIN = 2, stats::median),
		                     mdl.acc$kappa))
	options(warn=0)
	}  
    names(cv.ua) <- c(classes) 
      names(cv.pa)  <- c(classes) 
        names(mdl.ua) <- c(classes) 	  
          names(mdl.pa) <- c(classes) 	  
            names(mdl.oob) <- c("OOB", classes, "kappa") 
              names(cv.oob) <- c("OOB", classes, "kappa")	  
    acc <- list( cross.validation = list(cv.users.accuracy = cv.ua, 
	            cv.producers.accuracy=cv.pa, cv.oob = cv.oob), 
	            model = list(model.users.accuracy = mdl.ua, 
				             model.producers.accuracy=mdl.pa, 
			    			 model.oob=mdl.oob) )
	    if(bootstrap) acc[["boot.sample.size"]] <- boot.sample.size						 
      class(acc) <- c("rf.cv",  "classification", "list")
    return( acc )  							 
    }   
}  
