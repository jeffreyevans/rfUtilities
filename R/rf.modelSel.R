#' @title Random Forest Model Selection
#' @description Implements Murphy et al., (2010) Random Forests model selection 
#'              approach. 
#' 
#' @param ydata                 Y Data for model
#' @param xdata X               Data for model
#' @param imp.scale             Type of scaling for importance values (mir or se), 
#'                              default is mir
#' @param r                     Vector of importance percentiles to test i.e., 
#'                                seq(0,1,0.2)[2:5]
#' @param final.model           Run final model with selected variables (TRUE/FALSE)
#' @param seed                  Sets random seed in the R global environment. 
#'                              This is highly suggested.
#' @param parsimony             Threshold for competing model (0-1)
#' @param kappa                 Use the chance corrected kappa statistic rather 
#'                              than PCC
#' @param method                Use the fast C++ ranger implementation "Wright" or
#'                              original "Breiman" Fortran code  
#' @param pvalue                Calculate a p-value and filter parameters with this 
#'                              threshold

#' @param nperm                 Number of permutations to calculate p-value  
#' @param ...                   Additional arguments to pass to randomForest or ranger 
#'                              (e.g., ntree=1000, replace=TRUE, proximity=TRUE)
#'
#' @return \strong{A rf.modelSel class object with the following components:} 
#' \itemize{  
#'   \item   {"rf.final"} {Final selected model, if final = TRUE(randomForest 
#'                         model object)}
#'   \item   {"sel.vars"} {Final selected variables (vector)}
#'   \item   {"test"} {Validation parameters used on model selection (data.frame)}
#'   \item   {"sel.importance"} {Importance values for selected model (data.frame)}
#'   \item   {"importance"} {Importance values for all models (data.frame)}
#'   \item   {"parameters"} {Variables used in each tested model (list)}
#'   \item   {"scaling"} {Type of scaling used for importance}
#' }
#' 
#' @details
#' If you want to run classification, make sure that y is a factor, otherwise the 
#' randomForest model runs in regression mode For classification problems the model 
#' selection criteria is: smallest OOB error, smallest maximum within class error, 
#' and fewest parameters. For regression problems, the model selection criteria is 
#' largest percent variation explained, smallest MSE and fewest parameters.
#' 
#' @details
#' The "mir" scale option performs a row standardization and the "se" option 
#' performs normalization using the "standard errors" of the permutation-based 
#' importance measure. Both options result in a 0-1 range but, "se" sums to 1.
#' The scaled importance measures are calculated as: 
#'   mir = i/max(i) and se = (i / se) / ( sum(i) / se).
#' @details
#' The parsimony argument is the percent of allowable error surrounding 
#' competing models. For example, if there are two competing models, 
#' a selected model with 5 parameters and a competing model with 3 parameters, 
#' and parsimony = 0.05, if there is +/- 5% error in the fewer 
#' parameter model it will be selected at the final model. 
#' @details
#' If you specify the pvalue and nperm arguments then a permutation test is 
#' applied and parameters that do not meet the specified significance are removed 
#' before the model selection process. Please note that the p-value will be a  
#' function of the number of permutations. So a pvlaue=0.10 would be adequate for
#' nperm=99. 
#' @details
#' Using the kappa = TRUE argument will base error optimization on the kappa
#' rather than percent correctly classified (PCC). This will correct the PCC
#' for random agreement. The method = "Breiman" specifies the
#' use of the original Breiman Fortran code whereas "Wright" uses the C++ 
#' implementation from the ranger package (which exhibits a considerable 
#' improvement in speed).    
#' 
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @references 
#' Evans, J.S. and S.A. Cushman (2009) Gradient Modeling of Conifer Species 
#'   Using Random Forest. Landscape Ecology 5:673-683.
#' @references 
#' Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas 
#'   connectivity in Yellowstone National Park with landscape genetics. 
#'   Ecology 91:252-261
#' @references 
#' Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species 
#'   distribution and change using Random Forests CH.8 in Predictive Modeling 
#'   in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer
#'
#' @examples
#' require(randomForest)
#'   data(airquality)
#'   airquality <- na.omit(airquality)
#'
#'   xdata = airquality[,2:6]
#'   ydata = airquality[,1]
#'
#'  #### Regression example
#'  
#'  #### Using Breiman's original Fortran code from randomForest package
#'  ( rf.regress <- rf.modelSel(airquality[,2:6], airquality[,1], 
#'                              imp.scale="se") )
#'  
#'  #### Using Wright's C++ code from ranger package
#'  ( rf.regress <- rf.modelSel(airquality[,2:6], airquality[,1], 
#'                              method="Wright") )
#'
#'  #### Classification example
#'  ydata = as.factor(ifelse(ydata < 40, 0, 1))
#'  
#'   #### Using Breiman's original Fortran code from randomForest package
#'   ( rf.class <- rf.modelSel(xdata, ydata, ntree=1000) )
#'   
#'      # Use selected variables (same as final.model = TRUE
#'	  vars <- rf.class$selvars
#'      ( rf.fit <- randomForest(x=iris[,vars], y=iris[,"Species"]) )
#'   
#'      # Use results to select competing model
#'	  vars <- na.omit(as.character(rf.class$parameters[2,]))
#'      ( rf.fit <- randomForest(x=xdata[,vars], y=ydata) )   
#'   
#'   #### Using Wright's C++ code from ranger package
#'   ( rf.class <- rf.modelSel(xdata, ydata, method="Wright") )	
#'   	
#' \dontrun{
#'    # Using ranger package, filter p-values for classification
#'    ( rf.class <- rf.modelSel(xdata, ydata, method="Wright", 
#'                              pvalue=0.1, nperm=99, num.trees=1000) )
#'
#'   # Using ranger package, filter p-values for regression
#'   ( rf.class <- rf.modelSel(airquality[,1], ydata, method="Wright", 
#'                             pvalue=0.1, num.trees=1000) )
#' }
#'
#' @seealso \code{\link[randomForest]{randomForest}} for randomForest ... model options when method = "Breiman" 
#' @seealso \code{\link[ranger]{ranger}} for ranger ... model options when method = "Wright"
#' @seealso \code{\link[rfUtilities]{rf.ImpScale}} details on p-values
#'
#' @exportClass rf.modelSel
#' @export

rf.modelSel <- function(xdata, ydata, imp.scale = c("mir", "se"), r = c(0.25, 0.50, 0.75),  
                        final.model = FALSE, seed = NULL, parsimony = NULL, kappa = FALSE, 
						method = c("Breiman", "Wright"), pvalue=NULL, nperm = 99, ...) {    
  if(missing(ydata) | missing(xdata))
    stop("Both ydata and xdata arguments must be defined")
  if(!is.null(seed)) { set.seed(seed) }
  if(method ==  "Breiman" && !is.null(pvalue))
    warning("randomForest does not support importance p-values, argument will be ignored") 
  if(method ==  "Wright" && imp.scale == "se"){
    warning("ranger does not support standard error importance, defaulting to mir") 
      imp.scale == "mir"
   }
  r <- unique(c(0,r,1)) 
    expected <- seq(0,1,0.01)
      range.idx <- findInterval(expected, r, all.inside = TRUE) 
  			  
  # format ellipse ... 
  dots <- as.list(match.call(expand.dots = TRUE)[-1])
  rm.idx <- names(dots) %in% c("expand.dots", "xdata", "ydata", "imp.scale",   
    "r", "final.model", "seed", "parsimony", "kappa","method", "pvalue", "nperm")   
    if(length(rm.idx > 0)) dots <- dots[!rm.idx]
    if (!"y" %in% names(dots)) 
      dots[["y"]] <- ydata
    if (!"x" %in% names(dots))	  
	  dots[["x"]] <- xdata
    if(method[1] == "Breiman") {
      if (!"importance" %in% names(dots))
        dots[["importance"]] <- TRUE 
    } else if(method[1] == "Wright") {
      if (!"importance" %in% names(dots))  
       dots[["importance"]] <- "permutation"
    }
  
  RFtype <- is.factor(ydata) 

  ## classification ##
  if (RFtype == TRUE) {
    model.vars <- list()
    ln <- 0
    if(method[1] == "Breiman") {
	  print(names(dots))
	  rf.all <- do.call(randomForest::randomForest, dots)
	    model.vars[[ln <- ln + 1]] <- rownames(rf.all$importance)
          cm <- accuracy(rf.all$confusion[,1:(ncol(rf.all$confusion)-1)])
    } else if(method[1] == "Wright") {
      rf.all <- do.call(ranger::ranger, dots)
	    if(!is.null(pvalue)) {
		  imp_pval <- rf.ImpScale(rf.all, scaling="p", n=nperm)
		    dropped <- imp_pval[which(imp_pval$pvalue > pvalue),]$parameter
			if((length(dropped) > 0) == TRUE) {
			  cat("\n", "The following parameter(s) did not meet p-value threshold: ", 
			      dropped, "\n")
		          xdata <- xdata[,which(!names(xdata) %in% dropped)]
		       dots[["x"]] <- xdata
             rf.all <- do.call(ranger::ranger, dots)
             }			  
		}
	  model.vars[[ln <- ln + 1]] <- names(rf.all$variable.importance)
		if("probability" %in% names(dots)) {
		  cmat <- table(y, ifelse(rf.all$predictions[,2] < 0.50, 0, 1)) 
		} else {  
          cmat <- rf.all$confusion.matrix
		}
      cm <- accuracy(cmat)	  
    }

    if(kappa) {
      class.errors <- data.frame(1-t(c(cm$kappa, (cm$producers.accuracy/100))))
        names(class.errors)[1] <- "kappa"
    } else {		
      class.errors <- data.frame(1-t(c(cm$PCC, cm$producers.accuracy)/100))
        names(class.errors)[1] <- "pcc"
	}    

    # error, class.error, threshold, nparameters 
    errors <- data.frame(error = class.errors[,1], 
                         class.error = max(class.errors[,2:ncol(class.errors)]),
    					 threshold = 1, 
    					 nparameters = ncol(xdata))	
      nan.idx <- which(is.nan(errors))
	    if(length(nan.idx) > 0) errors[1,][nan.idx] <- 1
      inf.idx <- which(is.infinite(errors))
	    if(length(inf.idx) > 0) errors[1,][inf.idx] <- 0
		
	# build importance 	
	  imp <- rf.ImpScale(rf.all, scaling = imp.scale)
        if(imp.scale == "se") imp[,2] <- imp[,2] / max(imp[,2])
          imp <- imp[order(imp$importance),]
            imp$p <- findInterval(imp[,2], r, all.inside = TRUE)   
		miss.idx <- which(!sort(unique(findInterval(expected, r, all.inside=TRUE))) %in% 
		                  sort(unique(imp$p)))
		  if(miss.idx > 0) {
		    miss.range <- tapply(expected, findInterval(expected, r, 
		                         all.inside = TRUE), range)[miss.idx] 
		  					   
		    miss.range <- paste0("Missing importance ranges; ", unlist(lapply(miss.range, 
		                         FUN=function(x) paste0(x[1],"-",x[2]))))					   
              for(i in miss.range) message(i)
            }
	  }	
      for (p in unique(imp$p)) {
		sel.imp <- imp[imp$p == p,]
		  if(!exists("sel.vars")) {
		    sel.vars <- imp$parameter
		  }    
		   if(!any((sel.vars %in% sel.imp$parameter) == FALSE)) {
		     cat("\n", "Parameters match last threshold, skipping threshold set", p, "\n")
		     next
		   }		  
        sel.vars <- sel.imp$parameter
		np = length(sel.vars)
      if (length(sel.vars) > 1) {
        dots[["x"]] <- xdata[,sel.vars]	  
        if(method[1] == "Breiman") {
          rf.model <- do.call(randomForest::randomForest, dots)
            model.vars[[ln <- ln + 1]] <- rownames(rf.model$importance)
              cm <- accuracy(rf.model$confusion[,1:(ncol(rf.model$confusion)-1)])
          } else if(method[1] == "Wright") {
            rf.model <- do.call(ranger::ranger, dots)
              model.vars[[ln <- ln + 1]] <- names(rf.model$variable.importance)
	        if("probability" %in% names(dots)) {
	        	  cmat <- table(y, ifelse(rf.all$predictions[,2] < 0.50, 0, 1)) 
	        	} else {  
                  cmat <- rf.model$confusion.matrix
	        	}
	        	cm <- accuracy(cmat)	
          }
            if(kappa) {
              e <- data.frame(1-t(c(cm$kappa, (cm$producers.accuracy/100))))
                names(e)[1] <- "kappa"
            } else {
			  e <- data.frame(1-t(c(cm$PCC, cm$producers.accuracy)/100))
         	    names(e) <- c("pcc",names(cm$producers.accuracy))
	        } 
		   e <- data.frame(error = e[,1], 
                           class.error = max(e[,2:ncol(e)]),
    					   threshold = thres, 
    					   nparameters = np)
            nan.idx <- which(is.nan(e))
	          if(length(nan.idx) > 0) e[1,][nan.idx] <- 1
            inf.idx <- which(is.infinite(e))
	          if(length(inf.idx) > 0) e[1,][inf.idx] <- 0
            errors <- rbind(errors, e)
        } else {
          warning(paste0("The ", min(expected[which(range.idx %in% p)]), "-",  
		          max(expected[which(range.idx %in% p)]), 
				  " threshold has <= 1 parameter and cannot be evaluated") )
		}		 
      } 
	  
      n <- max(unlist(lapply(model.vars,FUN=length)))
        for(l in 1:length(model.vars)){
          x <- model.vars[[l]]
            length(x) <- n
          model.vars[[l]] <- x
        }
	model.vars <- as.data.frame(do.call(rbind, model.vars))
	  names(model.vars) <- paste0("parameter", 1:ncol(model.vars))
      errors <- data.frame(errors, model.vars)
        rownames(errors) <- 1:nrow(errors) 	  
          errors <- errors[order(errors$class.error, errors$error, errors$nparameters),]

	  if(!is.null(parsimony)) { 
  	    if(parsimony < 0.00000001 | parsimony > 0.9) stop( "parsimony must range 0-1")
          oob <- "TRUE"
          for(i in 2:nrow(errors)) {
            if( abs((errors[i,][1] - errors[1,][1]) / errors[1,][2]) <= parsimony  &
                abs( (errors[i,][2] - errors[1,][2] ) / errors[1,][3] ) <= parsimony ) 
			{ oob <- append(oob, "TRUE") } else { oob <- append(oob, "FALSE") }
              final <- errors[which( oob == "TRUE" ),]
          	  final <- min(final[final$nparameters == min(final$nparameters) ,]$threshold)
            }   
		} else {		
          final <- as.vector(errors[,"threshold"])[1]
        }	
	sel.post <- which( errors[,"threshold"] == final)	
	  sel.imp <- errors[sel.post,]
	    sel.vars <- as.character(sel.imp[,5:ncol(sel.imp)])
	      sel.vars <- sel.vars[!is.na(sel.vars)]
    errors <- rbind(errors[sel.post,],errors[-sel.post,])

  ## regression ##	
  } else if(RFtype == "FALSE") {
    model.vars <- list()
    ln <- 0
    if(method[1] == "Breiman") {
	  rf.all <- randomForest::randomForest(x=xdata, y=ydata, importance=TRUE, ...)
	    model.vars[[ln <- ln + 1]] <- rownames(rf.all$importance)
          e <- c(varexp=(stats::median(rf.all$rsq)), mse=mean(rf.all$mse))
    } else if(method[1] == "Wright") {
      rf.all <- ranger::ranger(x=xdata, y=ydata, importance="permutation", ...)
	    if(!is.null(pvalue)) {
		  imp_pval <- rf.ImpScale(rf.all, scaling="p", n=nperm)
		    dropped <- imp_pval[which(imp_pval$pvalue > pvalue),]$parameter
			if(length(dropped) > 0) {
			  cat("\n", "The following parameter(s) did not meet p-value threshold: ", dropped, "\n")
		      xdata <- xdata[,imp_pval[which(imp_pval$pvalue <= pvalue),]$parameter]
		      rf.all <- ranger::ranger(x=xdata, y=ydata, importance="permutation", ...)
             }			  
		}
	    model.vars[[ln <- ln + 1]] <- names(rf.all$variable.importance)
	  e <- c(varexp=rf.all$r.squared, mse=rf.all$prediction.error)
    }
	
    # varexp, mse, threshold, nparameters 
    errors <- data.frame(varexp = e[1], 
                         mse = e[2],
    					 threshold = 1, 
    					 nparameters = ncol(xdata))	
					
	# build importance 	
	imp <- rf.ImpScale(rf.all, scaling = imp.scale)
      if(imp.scale == "se") imp[,2] <- imp[,2] / max(imp[,2])
        imp <- imp[order(imp$importance),]
          imp$p <- findInterval(imp[,2], r, all.inside = TRUE)   
	    miss.idx <- which(!sort(unique(findInterval(seq(0,1,0.01), r, all.inside=TRUE))) %in% 
	                      sort(unique(imp$p)))
	    if(miss.idx > 0) {
	      expected <- seq(0,1,0.01)
	      miss.range <- tapply(expected, findInterval(expected, r, 
	                           all.inside = TRUE), range)[miss.idx] 
	  					   
	      miss.range <- paste0("Missing importance ranges; ", unlist(lapply(miss.range, 
	                           FUN=function(x) paste0(x[1],"-",x[2]))))
          for(i in miss.range) message(i)
        }

      for (p in unique(imp$p)) {
		sel.imp <- imp[imp$p == p,]
		  if(!exists("sel.vars")) {
		    sel.vars <- imp$parameter
		  }    
		   if(!any((sel.vars %in% sel.imp$parameter) == FALSE)) {
		     cat("\n", "Parameters match last threshold, skipping threshold", r[p] , "\n")
		     next
		   }		  
        sel.vars <- sel.imp$parameter
		np = length(sel.vars)
      if (length(sel.vars) > 1) {    
        if(method[1] == "Breiman") {
          rf.model <- randomForest::randomForest(x=xdata[,sel.vars], y=ydata, 
		                                         importance=TRUE, ...)
              model.vars[[ln <- ln + 1]] <- rownames(rf.model$importance)
			e <- c(varexp=(stats::median(rf.all$rsq)), mse=mean(rf.all$mse)) 
          } else if(method[1] == "Wright") {
            rf.model <- ranger::ranger(x=xdata[,sel.vars], y=ydata, 
			                           importance="permutation", ...)
              model.vars[[ln <- ln + 1]] <- names(rf.model$variable.importance) 
            e <- c(varexp=rf.all$r.squared, mse=rf.all$prediction.error)
          } 
		    e <- data.frame(varexp = e[1], 
                            mse = e[2],
    					    threshold = thres, 
    					    nparameters = np)
            errors <- rbind(errors, e)
        } else {
          warning(paste0("The ", min(expected[which(range.idx %in% p)]), "-",  
		          max(expected[which(range.idx %in% p)]), 
				  " threshold has <= 1 parameter and cannot be evaluated") )  
        }		
      } 
      n <- max(unlist(lapply(model.vars,FUN=length)))
        for(l in 1:length(model.vars)){
          x <- model.vars[[l]]
            length(x) <- n
          model.vars[[l]] <- x
        }
	model.vars <- as.data.frame(do.call(rbind, model.vars))
	  names(model.vars) <- paste0("parameter", 1:ncol(model.vars))
      errors <- data.frame(errors, model.vars)
        rownames(errors) <- 1:nrow(errors) 	 	  
          errors <- errors[order(errors$varexp, errors$mse, errors$nparameters),] 
	  if(!is.null(parsimony)) { 
  	    if(parsimony < 0.00000001 | parsimony > 0.9) stop( "parsimony must range 0-1")
          oob <- "TRUE"
          for(i in 2:nrow(errors)) {
            if( abs((errors[i,][1] - errors[1,][1]) / errors[1,][2]) <= parsimony  &
                abs( (errors[i,][2] - errors[1,][2] ) / errors[1,][3] ) <= parsimony ) 
			{ oob <- append(oob, "TRUE") } else { oob <- append(oob, "FALSE") }
              final <- errors[which( oob == "TRUE" ),]
          	  final <- min(final[final$nparameters == min(final$nparameters) ,]$threshold)
            }   
		} else {		
          final <- as.vector(errors[,"threshold"])[1]
        }	
	sel.post <- which( errors[,"threshold"] == final)	
	  sel.imp <- errors[sel.post,]
	    sel.vars <- as.character(sel.imp[,5:ncol(sel.imp)])
	      sel.vars <- sel.vars[!is.na(sel.vars)]
    errors <- rbind(errors[sel.post,],errors[-sel.post,])
  } # end of regression

    if (final.model == TRUE) {
	  dots[["x"]] <- xdata[,sel.vars] 
      if(method[1] == "Breiman") {
	    print(names(dots))
        rf.final <- do.call(randomForest::randomForest, dots)
      } else if(method[1] == "Wright") {
	    print(names(dots))
	    rf.final <- do.call(ranger::ranger, dots)
	  }
      mdl.sel <- list(rf.final = rf.final, selvars = sel.vars, test = errors, importance = imp, 
	                  sel.importance = sel.imp, parameters = model.vars, scaling = imp.scale)      
    } else {
      mdl.sel <- list(selvars = sel.vars, test = errors, importance = imp, sel.importance = sel.imp, 
	                  scaling = imp.scale, parameters = model.vars) 
    }
    class( mdl.sel ) <- c("rf.modelSel", "list")	
  return( mdl.sel )	
}
