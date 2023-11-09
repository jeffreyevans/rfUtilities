#' @title Spatial uncertainty
#' @description Calculates the standard errors and +/- 95% confidence
#'    
#' @param x           randomForest or ranger model object
#' @param y           terra raster object representing prediction
#' @param vars        multiband terra raster object used to predict y
#' @param b           number of Bootstrap replicates, defaults to model  
#' @param seed        Optional seed
#' @param out.raster  Optional raster to disk
#' @param ...         Additional arguments passed to writeRaster
#'
#' @return terra raster object with standard errors, lower and upper CI's.  
#'
#' @details 
#' Uses an Infinitesimal Jackknife (Wager et al., 2014) to calculate standard 
#' errors, lower and upper 95% confidence interval(s). The y argument is the 
#' spatial prediction raster from the model vars are a raster object, containing 
#' independent variables, used for the spatial estimates. The names of the
#' parameters used in the model must match vars and the extent of y and vars must 
#' be the same. The spatial prediction should  represent a continuous process
#' or probability and not a nominal prediction. 
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @references
#' Wager, S., T. Hastie, & B. Efron (2014) Confidence Intervals for Random 
#'   Forests: The Jackknife and the Infinitesimal Jackknife. J Mach Learn 
#'   Res 15:1625-1651.
#'
#' @examples
#' \dontrun{
#' library(randomForest)
#' library(terra)
#' library(ranger)
#'  
#' r <- rast(system.file("ex/logo.tif", package="terra"))   
#' 
#' # known presence and absence points
#' p <- matrix(c(48, 48, 48, 53, 50, 46, 54, 70, 84, 85, 74, 84, 95, 85, 
#'    66, 42, 26, 4, 19, 17, 7, 14, 26, 29, 39, 45, 51, 56, 46, 38, 31, 
#'    22, 34, 60, 70, 73, 63, 46, 43, 28), ncol=2)
#' a <- matrix(c(22, 33, 64, 85, 92, 94, 59, 27, 30, 64, 60, 33, 31, 9,
#'    99, 67, 15, 5, 4, 30, 8, 37, 42, 27, 19, 69, 60, 73, 3, 5, 21,
#'    37, 52, 70, 74, 9, 13, 4, 17, 47), ncol=2)
#' 
#' # extract values for points
#' xy <- rbind(cbind(1, p), cbind(0, a))
#' dat <- data.frame(cbind(pa=xy[,1], extract(r, xy[,2:3])))
#' 
#' # randomForest example
#' ( rfm <- randomForest(x=dat[,2:ncol(dat)], y=as.factor(dat$pa)) )
#'     e <- predict(r, rfm, type="prob")[[2]]
#' 	   names(e) <- "probs"
#' 
#' ( ci <- spatial.uncertainty(rfm, e, r) )
#' 
#' # Ranger example
#' ( rfm <- ranger(x=dat[,2:ncol(dat)], y=as.factor(dat$pa), 
#'                 probability = TRUE, num.trees = 501, 
#' 				   importance="permutation", write.forest = TRUE,  
#' 				   keep.inbag = TRUE) )
#'   rf.predict <- function(model, data) {
#'     as.numeric(ranger:::predict.ranger(model, data = data,
#'       type = "response")$predictions[,2])
#'   }				
#'   e <- predict(r, rfm, fun=rf.predict)
#'     names(e) <- "probs"
#' 	 
#' ( ci <- spatial.uncertainty(rfm, e, r) )
#'   plot(c(e, ci))
#'
#' } 
#'
#' @export spatial.uncertainty
spatial.uncertainty <- function(x, y, vars, b = NULL, seed=NULL, 
                                out.raster=NULL, ...) {
  if(missing(x))
    stop("x argument must be provided")
  if(missing(y))
    stop("y argument must be provided")
  if(missing(vars))
    stop("vars argument must be provided")	
  if (!inherits(x, c("randomForest", "ranger"))) 
    stop(deparse(substitute(y)), "is not randomForest or ranger class object")	
  if(inherits(x, "randomForest")) {
    if (x$type == "unsupervised") 
      stop("Unsupervised classification not supported")
    mtype <- tolower(x$type)
	  message("Model will be re-fit using ranger using same parameters")
  } else if(inherits(x, "ranger")) {
    mtype <- tolower(x$treetype)
  }		
  if (!inherits(y, "SpatRaster"))  
    stop(deparse(substitute(y)), " must be a terra SpatRaster class object") 
  if (!inherits(vars, "SpatRaster"))  
    stop(deparse(substitute(vars)), " must be a terra SpatRaster class object") 
  if(!terra::compareGeom(y, vars[[1]], res=TRUE))
    stop(deparse(substitute(y)), "and", deparse(substitute(vars)), "do not match")   
  if(!is.null(seed)) { set.seed(seed) }
  a <- as.list(x$call)[-1] 
    if(length(grep("~", a[[1]])) > 0) {
      dv <- all.vars(a[[1]])[1]
      xdata <- get(a$data) 	  
      ydata = xdata[,which(names(xdata) == dv[1])]	  
      xdata = xdata[,-which(names(xdata) == dv[1])]
    } else {  
      xdata = eval(a[[which(names(a) == "x")]]) 
      ydata = eval(a[[which(names(a) == "y")]]) 
    }
    if(is.null(b)) {
      mb <- c("ntree", "num.trees") 
      if(any(mb %in% names(a))) {
     b = a[[mb[which(mb %in% names(a))]]]
      } else {
        b = 501
      }	 
    }

	
  if(any(!names(xdata) %in% names(vars)))
    stop("Names in vars do not match model parameters")
  if(mtype == "regression") {
    mdl <- ranger::ranger(x=xdata, y=ydata, num.trees = b, 
                  importance="permutation", write.forest = TRUE, 
				  keep.inbag = TRUE)
  } else if(any(mtype %in% c("classification","probability estimation"))) {
    mdl <- ranger::ranger(x=xdata, y=ydata, probability = TRUE, 
                  num.trees = b, importance="permutation", 
				  write.forest = TRUE, keep.inbag = TRUE)
  }
  predict.se <- function(model, data) {
    as.numeric(ranger:::predict.ranger(model, data = data,
      type = "se", se.method = "infjack")$se[,2])
  }
  se <- terra::predict(vars[[names(xdata)]], mdl, fun=predict.se)
    se <- terra::mask(terra::crop(se, y),y)  
      ci95 <- c(se, y - (se * 1.96), y + (se * 1.96) )
        names(ci95) <- c("std.err", "lower.ci", "upper.ci")
	if(!is.null(out.raster)) {
      message("writing uncertainty raster to: ", out.raster)	
	    terra::writeRaster(ci95, out.raster, ...)
	}  
  return(ci95)		
}  

#( rfm <- ranger(pa ~ red + green + blue, data=dat, 
#                 probability = TRUE, num.trees = 501, 
# 				  importance="permutation", write.forest = TRUE,  
# 				  keep.inbag = TRUE) )

#( rfm <- ranger(pa ~ ., data=dat, probability = TRUE, 
#                num.trees = 501, importance="permutation", 
#				 write.forest = TRUE,  
# 				 keep.inbag = TRUE) )