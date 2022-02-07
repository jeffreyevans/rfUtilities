#' @title Spatial uncertainty
#' @description Calculates the standard errors and +/- 95% confidence
#'    
#' @param x           randomForest or ranger model object
#' @param y           raster object representing prediction
#' @param vars        raster stack/brick used to predict y
#' @param b           number of jackknife replicates 
#' @param seed        Optional seed
#' @param out.raster  Optional raster to disk
#' @param ...         Additional arguments passed to writeRaster
#'
#' @return raster stack object with standard errors, lower and upper CI's.  
#'
#' @details 
#' Uses an Infinitesimal Jackknife to calculate standard errors,
#' lower and upper 95% confidence interval(s). The y argument is the spatial
#' prediction and vars are a raster stack/brick object used for the spatial
#' estimates. The names of the parameters used in the model must match vars
#' and the extent of y and vars must be the same. The spatial prediction should
#' represent a continuous process or probability and not a nominal prediction. 
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
#' library(raster)
#' library(ranger)
#' 
#' r <- brick(system.file("external/rlogo.grd", package="raster"))
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
#'   e <- predict(r, rfm, type="prob", index=2)
#' 
#' ( ci <- spatial.uncertainty(rfm, e, r) )
#' 
#' # Ranger example
#' ( rfm <- ranger(x=dat[,2:ncol(dat)], y=as.factor(dat$pa), 
#'                 probability = TRUE, num.trees = 501, 
#' 				importance="permutation", write.forest = TRUE,  
#' 				keep.inbag = TRUE) )
#'   rf.predict <- function(model, data) {
#'     as.numeric(ranger:::predict.ranger(model, data = data,
#'       type = "response")$predictions[,2])
#'   }				
#'   e <- predict(vars, rfm, fun=rf.predict, progress="window")
#' 
#' ( ci <- spatial.uncertainty(rfm, e, r) )
#' } 
#'
#' @export spatial.uncertainty
spatial.uncertainty <- function(x, y, vars, b = 99, seed=NULL, 
                                out.raster=NULL, ...) {
  if(!any(class(y)[1] %in% c("RasterLayer"))) 
    stop("y is not a raster stack or brick object")
  if(!any(class(vars)[1] %in% c("RasterBrick","RasterStack")))
    stop("vars is not a raster stack or brick object")
  if(!raster::compareRaster(y,vars[[1]]))	
    stop("y and vars rasters do not match")	
  if (!inherits(x, c("randomForest", "ranger"))) 
    stop("x is not randomForest or ranger class object")	
  if(length(grep("~", x$call[[2]])) > 0)
      stop("This functions does not support a formula interface, 
	        please use x, y arguments")			
  if(inherits(x, "randomForest")) {
    if (x$type == "unsupervised") 
      stop("Unsupervised classification not supported")
    mtype <- tolower(x$type)
  } else if(inherits(x, "ranger")) {
    mtype <- tolower(x$treetype)
  }
  if(!is.null(seed)) { set.seed(seed) }
  a <- as.list(x$call)[-1] 
    xdata = eval(a[[which(names(a) == "x")]]) 
    ydata = eval(a[[which(names(a) == "y")]]) 
  if(any(!names(xdata) %in% names(vars)))
    stop("Names in vars do not match model parameters")
  if(mtype == "regression") {
    mdl <- ranger(x=xdata, y=ydata, num.trees = b, 
                  importance="permutation", write.forest = TRUE, 
				  keep.inbag = TRUE)
  } else if(any(mtype %in% c("classification","probability estimation"))) {
    mdl <- ranger(x=xdata, y=ydata, probability = TRUE, 
                  num.trees = b, importance="permutation", 
				  write.forest = TRUE, keep.inbag = TRUE)
  }
  predict.se <- function(model, data) {
    as.numeric(ranger:::predict.ranger(model, data = data,
      type = "se", se.method = "infjack")$se[,2])
  }
  se <- predict(vars[[names(xdata)]], mdl, fun=predict.se, progress="window")
    se <- mask(crop(se, y),y)  
      ci95 <- stack(se, y - (se * 1.96), y + (se * 1.96) )
        names(ci95) <- c("std.err", "lower.ci", "upper.ci") 
  if(!is.null(out.raster))
    raster::writeRaster(ci95, out.raster, ...)  
  return(ci95)		
}  
