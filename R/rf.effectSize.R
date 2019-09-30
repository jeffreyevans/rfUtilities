#' @title Random Forest effect size
#' @description Parameter effect size based on partial dependency (Cafri & Bailey, 2016)
#'       
#' @param x      A randomForest model object 
#' @param y      A vector represent the independent variable of intrests
#' @param ...    Arguments passed to the partial dependency function, requires x.var, pred.data,  
#
#' @return A vector (single value) of the parameter effect size 
#'
#' @note
#' Effect size based on partial dependency and parameter-weighted OLS (does not support factoral or dichotomous variables)  
#' The algorithm follows:
#'   1) Grow a forest
#'   2) Estimate partial dependence (for a single variable).
#'        a. Create datasets for all observation in the dataset only let them take on one value for the 
#'             variable of interest while keeping values of all other variables unchanged.
#'        b. Pass the dataset through each tree and average the predictions over the trees in the forest.
#'   3) Construct a point estimate of the proposed effect size by fitting a weighted least squares model 
#'        with response based on the tree-averaged predicted values obtained in Step 2, the explanatory variable 
#'         corresponding to the value used to generate each tree-averaged prediction, and weight based on the 
#'         frequency each value the explanatory variable takes on in the original data.
#'   4) For confidence intervals, repeat Steps 1-3 for as many bootstrap samples as desired
#' Modified partialPlot function uses distinct X values to construct partial dependence for non-factor variables
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references Cafri, G., B.A. Bailey (2016) Understanding Variable Effects from Black Box Prediction: Quantifying Effects in Tree Ensembles Using Partial Dependence. Journal of Data Science 14:67-96
#'
#' @examples 
#'  library(randomForest)
#'    data(airquality)
#'    airquality <- na.omit(airquality)
#'      fit.reg <- randomForest(Ozone ~ ., data=airquality)
#'  
#'  # Parameter effect sizes	
#'  rf.effectSize(fit.reg, y = airquality$Solar.R, pred.data = airquality, x.var = Solar.R)
#'  rf.effectSize(fit.reg, y = airquality$Wind, pred.data = airquality, x.var = Wind)
#'  rf.effectSize(fit.reg, y = airquality$Temp, pred.data = airquality, x.var = Temp)
#'  rf.effectSize(fit.reg, y = airquality$Month, pred.data = airquality, x.var = Month)
#'  rf.effectSize(fit.reg, y = airquality$Day, pred.data = airquality, x.var = Day)
#'  	
#' \dontrun{
#'  # Bootstrap of effect size for Wind and Temp
#'  B = 999
#'  n = nrow(airquality)
#'  es.boot.wind <- vector()
#'  es.boot.temp <- vector()
#'    for(i in 1:B) {
#'      boot.samples <- airquality[sample(1:nrow(airquality), n, replace = TRUE),]
#'        fmla <- stats::as.formula(paste(paste("Ozone", "~", sep=""), paste(".", collapse= "")))   
#'          fit <- randomForest(fmla, data = boot.samples)
#'      es.boot.wind <- append(es.boot.wind, rf.effectSize(fit, y = boot.samples$Wind, 
#'  	                       pred.data = boot.samples, x.var = Wind))
#'      es.boot.temp <- append(es.boot.temp, rf.effectSize(fit, y = boot.samples$Temp, 
#'  	                       pred.data = boot.samples,x.var = Temp))
#'    }        
#'   se <- function(x) sqrt(var(x, na.rm = TRUE) / length(na.omit(x)))
#'     cat("Bootstrap variance for Wind:", var(es.boot.wind), "\n")
#'       cat("Bootstrap standard error for Wind:", se(es.boot.wind), "\n","\n")
#'     cat("Bootstrap variance for Temp:", var(es.boot.temp), "\n")
#'   cat("Bootstrap standard error for Temp:", se(es.boot.temp), "\n")
#' 
#'  # Confidence intervals of Bootstrap of effect size for Wind
#'  p=0.95
#'  y <- sort(es.boot.wind) 
#'  x <- 1:length(y)
#'    plx <- stats::predict(stats::loess(y ~ x), se=TRUE)
#'    lci = plx$fit - stats::qt(p, plx$df) * plx$se.fit
#'    uci = plx$fit + stats::qt(p, plx$df) * plx$se.fit
#'        graphics::plot(x, y, type="n", main="Effect size Bootstrap CI for Wind", 
#'  	    sub=paste("confidence intervals at", p))
#'        graphics::polygon(c(x,rev(x)), c(lci, rev(uci)), col="grey86")
#'        graphics::points(x, y, pch=20, cex=0.70)
#'        graphics::lines(x, plx[["fit"]], lty=3)
#'  
#'  # Confidence intervals of Bootstrap of effect size for Temp
#'  p=0.95
#'  y <- sort(es.boot.temp) 
#'  x <- 1:length(y)
#'    plx <- stats::predict(stats::loess(y ~ x), se=TRUE)
#'    lci = plx$fit - stats::qt(p, plx$df) * plx$se.fit
#'    uci = plx$fit + stats::qt(p, plx$df) * plx$se.fit
#'        graphics::plot(x, y, type="n", main="Effect size Bootstrap CI for Temp", 
#'  	    sub=paste("confidence intervals at", p))
#'        graphics::polygon(c(x,rev(x)), c(lci, rev(uci)), col="grey86")
#'        graphics::points(x, y, pch=20, cex=0.70)
#'        graphics::lines(x, plx[["fit"]], lty=3)
#'  	  
#'  # Plot bootstrap of wind effect size
#'  pdf <- density(es.boot.wind)
#'  plot(pdf, type="n", main="Bootstrap of effect size wind (n=99)", 
#'       ylab="p", xlab="effect size")
#'    polygon(pdf, col="grey")
#'    abline(v=mean(es.boot.wind))
#'      abline(v=mean(es.boot.wind)-sd(es.boot.wind), col="blue", lty=3)
#'      abline(v=mean(es.boot.wind)+sd(es.boot.wind), col="blue", lty=3)	  
#'  
#'  # Plot bootstrap of temp effect size	
#'  pdf <- density(es.boot.temp)
#'  plot(pdf, type="n", main="Bootstrap of effect size temp (n=99)", 
#'       ylab="p", xlab="effect size")
#'    polygon(pdf, col="grey")	
#'    abline(v=mean(es.boot.temp))
#'      abline(v=mean(es.boot.temp)-sd(es.boot.temp), col="blue", lty=3)
#'      abline(v=mean(es.boot.temp)+sd(es.boot.temp), col="blue", lty=3)
#' }
#'
#' @seealso \code{\link[randomForest]{partialPlot}} for ... options
#' @seealso \code{\link[randomForest]{randomForest}} for randomForest details
#'
#' @export
rf.effectSize <- function(x, y, ...) {
  if(!any(class(x) %in% c("randomForest","list"))) stop("x is not a randomForest object")
    if(is.factor(y) == TRUE | is.character(y) == TRUE) stop("factorial or dichotomous variables not supported")
  dots <- as.list(match.call(expand.dots = TRUE)[-1])    
    dots[["x"]] <- x
      if (is.null(dots[["x.var"]]) & "x.var" %in% names(dots) == FALSE) stop("x.var must be defined") 
      if (is.null(dots[["pred.data"]]) & "main" %in% names(dots) == FALSE) stop("pred.data must be defined")
      if(x$type == "classification") {	  
        if (is.null(dots[["which.class"]]) & "which.class" %in% names(dots) == FALSE) dots[["which.class"]] <- 1
	  }
      if (is.null(dots[["plot"]]) & "plot" %in% names(dots) == FALSE) dots[["plot"]] <-  FALSE  
    partialPlot.es <- function (x, pred.data, x.var, which.class, w, plot = FALSE, add = FALSE,
                                n.pt = min(length(unique(pred.data[, xname])), 51), rug = TRUE,
                                xlab = deparse(substitute(x.var)), ylab = "", main = paste("Partial Dependence on",
                                deparse(substitute(x.var))), ...) {
      classRF <- x$type != "regression"
        if (is.null(x$forest)) stop("The randomForest object must contain the forest.\n")
          x.var <- substitute(x.var)
            xname <- if (is.character(x.var)) x.var
        else {
          if (is.name(x.var))
          deparse(x.var)
        else {
          eval(x.var)
          }
        }
          xv <- pred.data[, xname]
          n <- nrow(pred.data)
      if (missing(w)) w <- rep(1, n)
      if (classRF) {
        if (missing(which.class)) {
          focus <- 1
        }
        else {
          focus <- charmatch(which.class, colnames(x$votes))
        if (is.na(focus))
          stop(which.class, "is not one of the class labels.")
        }
      }
      if (is.factor(xv) && !is.ordered(xv)) {
        x.pt <- levels(xv)
        y.pt <- numeric(length(x.pt))
      for (i in seq(along = x.pt)) {
        x.data <- pred.data
        x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
          if (classRF) {
            pr <- stats::predict(x, x.data, type = "prob")
            y.pt[i] <- stats::weighted.mean(log(ifelse(pr[, focus] >
                    0, pr[, focus], .Machine$double.eps)) - rowMeans(log(ifelse(pr >
                    0, pr, .Machine$double.eps))), w, na.rm = TRUE)
          }
        else y.pt[i] <- stats::weighted.mean(stats::predict(x, x.data), w, na.rm = TRUE)
      }
    
      if (add) {
        graphics::points(1:length(x.pt), y.pt, type = "h", lwd = 2, ...)
      }
      else {
        if (plot)
          graphics::barplot(y.pt, width = rep(1, length(y.pt)), col = "blue",
                  xlab = xlab, ylab = ylab, main = main, names.arg = x.pt, ...)
        }
      }
      else {
        if (is.ordered(xv))
          xv <- as.numeric(xv)
          x.pt <- sort(unique(xv))
          y.pt <- numeric(length(x.pt))
      for (i in seq(along = x.pt)) {
        x.data <- pred.data
        x.data[, xname] <- rep(x.pt[i], n)
          if (classRF) {
            pr <- stats::predict(x, x.data, type = "prob")
            y.pt[i] <- stats::weighted.mean(log(ifelse(pr[, focus] ==
                       0, .Machine$double.eps, pr[, focus])) - rowMeans(log(ifelse(pr ==
                       0, .Machine$double.eps, pr))), w, na.rm = TRUE)
          }
        else {
         y.pt[i] <- stats::weighted.mean(stats::predict(x, x.data),
         w, na.rm = TRUE)
        }
      }
        if (add) { graphics::lines(x.pt, y.pt, ...) } else {
      if (plot)
        plot(x.pt, y.pt, type = "l", xlab = xlab, ylab = ylab,
             main = main, ...)
      }
        if (rug && plot) {
          if (n.pt > 10) {
            rug(stats::quantile(xv, seq(0.1, 0.9, by = 0.1)), side = 1)
        }
       else {
          rug(unique(xv, side = 1))
          }
        }
      }
    invisible(list(x = x.pt, y = y.pt))
    }
    partial.dep <- do.call("partialPlot.es", dots)
	  d <- dots[["pred.data"]]
      xname <- dots[["x.var"]]  	  
      w <- table( y )								
        effect.size <- stats::lm(y ~ x, data = partial.dep, weights = w )
        a <- stats::coefficients(effect.size)[2]
	names(a) <- paste("Effect size for", dots[["x.var"]], sep=" ")
  return( a )
}
