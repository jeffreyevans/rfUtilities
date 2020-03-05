#' @title Cross-validation for Random Forest models 
#' @description Implements a permutation based cross-validation test  
#'              for classification or regression Random Forests models
#'    
#' @param x                 A randomForest or ranger object
#' @param p                 Proportion data withhold (default p=0.10)
#' @param n                 Number of cross validations (default n=99)
#' @param seed              Sets random seed in R global environment
#' @param normalize         (FALSE/TRUE) For regression, should rmse, mbe and mae be normalized 
#'                          using (max(y) - min(y))
#' @param bootstrap         (FALSE/TRUE) Should a bootstrap sampling be applied. If FALSE, 
#'                          an n-th percent withhold will be conducted
#' @param p.threshold       If ranger probability forest, threshold to use in validation 
#' @param trace             Print iterations
#'
#' @return  For classification a "rf.cv"", "classification" class object with the following 
#'          components:
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
#' \item  y.rmse         Root Mean Squared Error (observed vs. predicted) from each Bootstrap 
#'                       iteration (cross-validation)    
#' \item  y.mbe          Mean Bias Error from each Bootstrapped model
#' \item  y.mae          Mean Absolute Error from each Bootstrapped model
#' \item  D              Test statistic from Kolmogorov-Smirnov distribution Test (y and 
#'                       estimate)
#' \item  p.val          p-value for Kolmogorov-Smirnov distribution Test (y and estimate)
#' \item  model.mse      Mean Squared Error from each Bootstrapped model
#' \item  model.varExp   Percent variance explained from each Bootstrapped model   
#'  }
#'
#' @note
#' Please note that previous versions of this function required ydata, xdata and   
#' "..." arguments that are no longer necessary. The model object is now used in 
#' obtaining the data and arguments used in the original model  
#'
#' @details
#' For classification problems, the cross-validation statistics are based on the 
#' prediction error on the withheld data: Total observed accuracy represents the 
#' percent correctly classified (aka, PCC) and is considered as a naive measure 
#' of agreement.
#' 
#' The diagonal of the confusion matrix represents correctly classified observations 
#' where off-diagonals represent cross-classification error. The primary issue with 
#' this evaluation is that does not reveal if error was evenly distributed between 
#' classes.
#'
#' To represent the balance of error one can use omission and commission statistics such 
#' as estimates of users and producers accuracy. User's accuracy corresponds to error of 
#' commission (inclusion), observations being erroneously included in a given class.
#'
#' The commission errors are represented by row sums of the matrix. Producer's accuracy 
#' corresponds to error of omission (exclusion), observations being erroneously excluded 
#' from a given class. The omission errors are represented by column sums of the matrix.
#'
#' None of the previous statistics account for random agreement influencing the accuracy 
#' measure. The kappa statistic is a chance corrected metric that reflects the difference 
#' between observed agreement and agreement expected by random chance. A kappa of k=0.85 
#' would indicate that there is 85% better agreement than by chance alone.  
#'   \itemize{ 
#'   \item   pcc = [Number of correct observations / total number of observations] 
#'   \item   pcc = [Number of correct observations / total number of observations] 
#'   \item   producers accuracy =  [Number of correct / total number of correct and 
#'                                 omission errors] 
#'   \item   k = (observed accuracy - chance agreement) / (1 - chance agreement) where; change 
#'                agreement = sum[product of row and column totals for each class] 
#'    }
#' For regression problems, a Bootstrap is constructed and the subset models MSE and percent 
#' variance explained is reported. Additional, the RMSE between the withheld response variable 
#' (y) and the predicted subset model   
#'
#' @author Jeffrey S. Evans <jeffrey_evans<at>tnc.org>
#'
#' @references 
#'   Evans, J.S. and S.A. Cushman (2009) Gradient Modeling of Conifer Species Using Random 
#'     Forest. Landscape Ecology 5:673-683.
#' @references 
#'   Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas connectivity in 
#'     Yellowstone National Park with landscape genetics. Ecology 91:252-261
#' @references 
#'   Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution 
#'     and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds 
#'     Drew, CA, Huettmann F, Wiersma Y. Springer 
#' 
#' @examples 
#' \dontrun{
#' library(randomForest)
#' library(ranger)
#'
#' data(airquality)
#' airquality <- na.omit(airquality)
#' yclass = as.factor(ifelse(airquality[,1] < 40, 0, 1))
#' 
#' # regression with ranger
#' rf.mdl <- ranger(x = airquality[,2:6], y = airquality[,1])
#'   ( rf.cv <- rf.crossValidation(rf.mdl, p=0.10) )
#' 
#'   # plot results
#'   par(mfrow=c(2,2))
#'     plot(rf.cv)  
#'     plot(rf.cv, stat = "mse")
#'     plot(rf.cv, stat = "var.exp")
#'     plot(rf.cv, stat = "mae")
#' 
#' # regression with randomForest
#' rf.mdl <- randomForest(airquality[,2:6], airquality[,1])
#'   ( rf.cv <- rf.crossValidation(rf.mdl, p=0.10) )
#' 
#' # classification with ranger
#' rf.mdl <- ranger(x = airquality[,2:6], y = yclass)
#'   ( rf.cv <- rf.crossValidation(rf.mdl, p=0.10) )
#' 
#'     # Plot cross validation versus model producers accuracy
#'     par(mfrow=c(1,2)) 
#'       plot(rf.cv, type = "cv", main = "CV producers accuracy")
#'       plot(rf.cv, type = "model", main = "Model producers accuracy")
#'     
#'     # Plot cross validation versus model oob
#'     par(mfrow=c(1,2)) 
#'       plot(rf.cv, type = "cv", stat = "oob", main = "CV oob error")
#'       plot(rf.cv, type = "model", stat = "oob", main = "Model oob error")	 
#' 
#' # classification with randomForest
#' rf.mdl <- randomForest(x = airquality[,2:6], y = yclass)
#'   ( rf.cv <- rf.crossValidation(rf.mdl, p=0.10) )
#' 
#' # multi-class classification
#' data(iris)
#'   iris$Species <- as.factor(iris$Species)    	
#' ( rf.mdl <- randomForest(iris[,1:4], iris[,"Species"], ntree=501) )
#'   ( rf.cv <- rf.crossValidation(rf.mdl) )
#'
#' }	 
#'	  
#' @seealso \code{\link[randomForest]{randomForest}} for randomForest details 
#' @seealso \code{\link[ranger]{ranger}} for ranger details 
#'
#' @exportClass rf.cv
#' @export 	
rf.crossValidation <- function(x, p=0.10, n=99, seed=NULL, normalize = FALSE, 
                               bootstrap = FALSE, p.threshold = 0.60, 
							   trace = FALSE) {
  if (!inherits(x, c("randomForest", "ranger"))) 
    stop("x is not randomForest class object")
  if(length(grep("~", x$call[[2]])) > 0)
      stop("This package does not support a formula interface, please use x, y arguments")
  if(inherits(x, "randomForest")) {
    if (x$type == "unsupervised") 
      stop("Unsupervised classification not supported")
    mtype <- tolower(x$type)
  } else if(inherits(x, "ranger")) {
    mtype <- tolower(x$treetype)
  }
  if(!is.null(seed)) { set.seed(seed) }
  if(bootstrap) 
    cat("Bootstrap sampling is being applied,", paste("p",p,sep="="), 
	    "argument is ignored", "\n")
 
  # formating call and pulling data
  a <- as.list(x$call)[-1] 
    xdata = eval(a[[which(names(a) == "x")]]) 
    ydata = eval(a[[which(names(a) == "y")]]) 
    dat <- data.frame(y=ydata, xdata) 
	
  #*********************************	
  if(mtype == "regression") {
  #*********************************
  cat("running:", mtype, "cross-validation", "with", n, "iterations", "\n")
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
    if(bootstrap) boot.sample.size <- rep(NA, n)	
      sample.size = round( (length(ydata) * p), digits=0)
      #**************************
      # cross-validation for loop
	  #**************************
	  for(i in 1:n) {
	    if(trace) cat("running iteration:", i, "of", n, "\n")
		# Draw random sample		
	    if(!bootstrap) {
	      sidx <- sample(1:nrow(dat), sample.size)   
	        dat.sub <- dat[-sidx,]                   
            dat.cv <- dat[sidx,]		            
	    } else {	
	      dat.sub <- dat[sample(1:nrow(dat), replace=TRUE),]              
          dat.cv <- dat[which(!rownames(dat) %in% rownames(dat.sub)),]	   
        }		
		a[["y"]] <- dat.sub[,"y"]
		a[["x"]] <- dat.sub[,2:ncol(dat.sub)]
		if(inherits(x, "ranger")) {    
          rf.fit <- do.call(ranger::ranger, a)
 		    model.mse[i] <- rf.fit$prediction.error
	        model.varExp[i] <- rf.fit$r.squared
            prd <- stats::predict(rf.fit,data = dat.cv[,2:ncol(dat.cv)])$predictions 			

        } else if(inherits(x,"randomForest")) {
          rf.fit <- do.call(randomForest::randomForest, a)
 		    model.mse[i] <- rf.fit$mse[length(rf.fit$mse)]
	        model.varExp[i] <- round(100*rf.fit$rsq[length(rf.fit$rsq)], digits=2) 
            prd <- stats::predict(rf.fit, newdata = dat.cv[,2:ncol(dat.cv)]) 			
        }
	      y.rmse[i] <- rmse(dat.cv[,"y"], prd, norm=normalize) 
          y.mbe[i] <- mbe(dat.cv[,"y"], prd, norm=normalize) 
		  y.mae[i] <- mae(dat.cv[,"y"], prd, norm=normalize) 
          ks.p[i] <- as.numeric(ks(dat.cv[,"y"], prd, s=2))  
          ks.d[i] <- as.numeric(ks(dat.cv[,"y"], prd, s=1))
		if(bootstrap) boot.sample.size[i] <- nrow(dat.cv) 
	  }		   
    if(inherits(x, "ranger")) {
      fit.var.exp = x$r.squared
	  fit.mse = x$prediction.error
	} else if(inherits(x,"randomForest")) {
	  fit.var.exp = round(100*x$rsq[length(x$rsq)], digits=2) 
	  fit.mse = stats::median(x$mse)  
    }
	r.cv <- list(fit.var.exp=fit.var.exp, 
	             fit.mse=fit.mse, y.rmse = y.rmse, 
	    		 y.mbe = y.mbe, y.mae = y.mae, D = ks.d, 
				 p.val = ks.p, model.mse = model.mse, 
				 model.varExp = model.varExp )	
      class(r.cv) <- c("rf.cv", "regression")
	  
    #*********************************	
    } else if(any(mtype %in% c("classification","probability estimation"))) {
    #*********************************		
    cat("running:", mtype, "cross-validation", "with", n, "iterations", "\n")
		if(inherits(x, "ranger")) { 
		  if("probability estimation" == mtype) {  
		    ppred <- ifelse(x$predictions[,2] < p.threshold, 0, 1)
	      } else {
		    ppred <- x$predictions 
		  }
		  cm <- accuracy(ydata, ppred) 
        } else if(inherits(x,"randomForest")) {
		  cm <- accuracy(ydata, x$predicted)
		}
      mdl.oob <- cm$PCC	
	  mdl.kappa <- cm$kappa	
	  mdl.ua <- cm$users.accuracy 
	  mdl.pa <- cm$producers.accuracy 
	  model.error = data.frame(model.users.accuracy = mdl.ua, 
	                           model.producers.accuracy=mdl.pa, 
			                   model.oob=mdl.oob, model.kappa=mdl.kappa) 
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
	  cv.oob <- c(0,0,0)
      #**************************
      # cross-validation for loop
	  #**************************	  
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
          tx <- dat[sample(1:nrow(dat), replace=TRUE),]
            ty <- tx$y 
              tx <- tx[,2:ncol(tx)]
              mx <- dat[-which(!rownames(dat) %in% rownames(tx)),]
            my <- mx$y
          mx <- mx[,2:ncol(mx)] 
       }
	   
	   	dat.sub <- data.frame(y=my, mx)
        dat.cv <- data.frame(y=ty, tx) 
		a[["y"]] <- dat.sub[,"y"]
		a[["x"]] <- dat.sub[,2:ncol(dat.sub)]
		if(inherits(x, "ranger")) {    
          rf.fit <- do.call(ranger::ranger, a)
		    prd <- stats::predict(rf.fit,data = dat.cv[,2:ncol(dat.cv)])$predictions
		      if("probability estimation" == mtype) {	
                prd <- ifelse(prd[,2] < p.threshold, 0, 1)
			  } 
			cv.acc <- accuracy(prd,  ty)
	        mdl.oob <- rf.fit$prediction.error
        } else if(inherits(x,"randomForest")) {
          rf.fit <- do.call(randomForest::randomForest, a)
            prd <- stats::predict(rf.fit, newdata = dat.cv[,2:ncol(dat.cv)])
			cv.acc <- accuracy(prd,  ty)
			mdl.oob <- c(apply(rf.fit$err.rate, MARGIN = 2, stats::median))[1] 
        }
	    if(bootstrap) boot.sample.size <- append(boot.sample.size, length(my))
		cv.oob <- rbind(cv.oob, c(mdl.oob, cv.acc$PCC, cv.acc$kappa))	
 		  cv.ua <- rbind(cv.ua, cv.acc$users.accuracy) 
		    cv.pa <- rbind(cv.pa, cv.acc$producers.accuracy) 
		  #if(!length(classes) == length(unique(my))) {
	      #### add code to check for presence of classes and set 
		  ####  to NA if absent classes occur   
		  #}
	}  
    names(cv.ua) <- c(classes) 
      names(cv.pa)  <- c(classes)
        cv.oob <- as.data.frame(cv.oob)	  
          names(cv.oob) <- c("Model.PCC", "CV.PCC", "CV.kappa")
		    cv.oob <- cv.oob[-1,]
              rownames(cv.oob) <- 1:nrow(cv.oob)		  
    r.cv <- list(cv.users.accuracy = cv.ua, cv.producers.accuracy=cv.pa, 
	             cv.oob = cv.oob, model.error = model.error)
	    if(bootstrap) r.cv[["boot.sample.size"]] <- boot.sample.size
      class(r.cv) <- c("rf.cv",  "classification")
  }   
  return( r.cv )
}  
