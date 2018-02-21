#' @title Random Forest probability scaled partial dependency plots
#' @description Produces partial dependency plots with probability distribution based on scaled margin distances.
#'       
#' @param x              Object of class randomForest 
#' @param pred.data      Training data.frame used for constructing the plot, 
#' @param xname          Name of the variable for calculating partial dependence 
#' @param which.class    The class to focus on
#' @param w              Weights to be used in averaging (if not supplied, mean is not weighted) 
#' @param prob           Scale distances to probabilities
#' @param plot           (TRUE/FALSE) Plot results 
#' @param smooth         c(spline, loess) Apply spline.smooth or loess to 
#' @param conf           (TRUE/FALSE) Should confidence intervals be calculated for smoothing
#' @param smooth.parm    An appropriate smoothing parameter passed to loess or smooth.spline
#' @param pts            (FALSE/TRUE) Add raw points
#' @param raw.line       (FALSE/TRUE) Plot raw line (non-smoothed)
#' @param rug            Draw hash marks on plot representing deciles of x
#' @param n.pt           Number of points on the grid for evaluating partial dependence.
#' @param xlab           x-axis plot label
#' @param ylab           y-axis plot label
#' @param main           Plot label for main
#' @param ...            Additional graphical parameters passed to plot
#
#' @return A list class object with fit x,y. If smooth=c("spline","loess") y represents smoothed scaled margin distance values 
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer 
#' @references Baruch-Mordo, S., J.S. Evans, J. Severson, J. D. Naugle, J. Kiesecker, J. Maestas, & M.J. Falkowski (2013) Saving sage-grouse from the trees: A proactive solution to reducing a key threat to a candidate species Biological Conservation 167:233-241  
#'
#' @examples 
#'  require(randomForest)
#'    data(iris)
#'    iris.rf <- randomForest(iris[,1:4], iris[,5])		
#'  	 
#'  # plot all parameters	 
#'  par(mfrow=c(2,2))
#'    for(i in names(iris)[1:4]) {     
#'      rf.partial.prob(iris.rf, iris, i, "setosa", smooth="spline", 
#'                      n.pt=70, smooth.parm = 0.5)
#'     }
#'
#'  # Plot spline and loess smoothing for one parameter, with raw points and line
#'  par(mfrow=c(1,2))	 
#'    rf.partial.prob(x = iris.rf, pred.data = iris, xname = "Sepal.Length", 
#'                    which.class = "setosa", smooth = "spline", smooth.parm = 0.5,
#'    				  n.pt = 70, pts = TRUE, raw.line = TRUE, rug = TRUE)
#'    				
#'    rf.partial.prob(x = iris.rf, pred.data = iris, xname = "Sepal.Length", 
#'                    which.class = "setosa", smooth = "loess", smooth.parm = 0.20,
#'    				  n.pt = 70, pts = TRUE, raw.line = TRUE, rug = TRUE)
#'
#' @seealso \code{\link[stats]{smooth.spline}} for smooth.spline details on spar smoothing argument
#' @seealso \code{\link[stats]{loess}} for loess details of span smoothing argument 
#'
#' @export
rf.partial.prob <- function(x, pred.data, xname, which.class, w, prob=TRUE, plot=TRUE,
                            smooth, conf = TRUE, smooth.parm = NULL, pts = FALSE, 
							raw.line = FALSE, rug=FALSE, n.pt, xlab, ylab, main, ...) {  
    if(!any(class(x) %in% c("randomForest","list"))) stop("x is not a randomForest object")
	  if (is.null(x$forest)) stop("Object does not contain an ensemble \n")
	    if(!x$type != "regression")	stop("Regression not supported \n")	   
	      if (missing(which.class)) stop("Class name missing \n")
	    if (missing(xname)) stop("X Variable name missing \n")
	  if (missing(x)) stop("randomForest object missing \n")
    if (missing(pred.data)) stop("New data missing \n")
	focus <- charmatch(which.class, colnames(x$votes))
    if (is.na(focus)) stop(which.class, "is not one of the class labels")
    xv <- pred.data[, xname]
	n <- nrow(pred.data)
	  if(missing(n.pt)) n.pt <- min(length(unique(pred.data[, xname])), 51)
	  if (missing(w)) w <- rep(1, n)
	scale.dist <- function(d) { return( (exp(d) - min(exp(d))) / (max(exp(d)) - min(exp(d))) ) }
 
    if(is.factor(xv) && !is.ordered(xv)) {
      x.pt <- levels(xv)
      y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
          x.data <- pred.data
          x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
          pr <- stats::predict(x, x.data, type = "prob")
          y.pt[i] <- stats::weighted.mean(log(ifelse(pr[, focus] > 0, pr[, focus],
		              .Machine$double.eps)) - rowMeans(log(ifelse(pr > 0, pr, 
					  .Machine$double.eps))), 
					  w, na.rm=TRUE)
      }  
      if( prob == TRUE) { y.pt <- scale.dist(y.pt) } 
      if (plot) {
        if(missing(xlab)) xlab=xname
	    if(missing(ylab)) ylab=which.class
	    if(missing(main)) main="Partial Dependency Plot"
	      graphics::barplot(y.pt, width=rep(1,length(y.pt)), col="blue",
                            xlab=xlab, ylab=ylab, main=main,
                            names.arg=x.pt, ...)
      }				  
    }
	
    if (is.numeric(xv)) {
      x.pt <- seq(min(xv), max(xv), length=n.pt)
      y.pt <- numeric(length(x.pt))
        for (i in seq(along=x.pt)) {
          x.data <- pred.data
          x.data[, xname] <- rep(x.pt[i], n)
          pr <- stats::predict(x, x.data, type = "prob")
          y.pt[i] <- stats::weighted.mean(log(ifelse(pr[, focus] == 0, .Machine$double.eps,
               pr[, focus])) - rowMeans(log(ifelse(pr == 0, .Machine$double.eps, pr))),
               w, na.rm=TRUE)
        }
		
    if(prob == TRUE) { y.pt <- scale.dist(y.pt) }
	if (plot == TRUE) {
      if(missing(xlab)) xlab = xname
	    if(missing(ylab)) ylab = "probability"
	      if(missing(main)) main = paste("Partial Dependency Plot for Class",
		                                 which.class, sep=" - ")
      if(!missing(smooth)) {
	    if(smooth == "spline") {	
	      fit <- stats::smooth.spline(x.pt, y.pt, spar = smooth.parm, all.knots = TRUE)
          fit.y <- fit$y
            if(prob == TRUE) {
              fit.y[fit.y < 0] <- 0		
		      fit.y[fit.y > 1] <- 1
            }			  
            if(conf == TRUE) {
	# Using jackknifed residuals for upper and lower 95% confidence interval 
                res <- (fit$yin - fit$y)/(1-fit$lev)    
                  sigma <- sqrt(stats::var(res))              
                    upper <- fit$y + 2.0 * sigma * sqrt(fit$lev)   
                    lower <- fit$y - 2.0 * sigma * sqrt(fit$lev)  
		  }
	  } else if (smooth == "loess") {
        options(warn=-1)
          if(is.null(smooth.parm)) smooth.parm = 0.75		
		    fit <- stats::predict(stats::loess(y.pt ~ x.pt, span = smooth.parm), 
			                      se = TRUE, degree = 2)
		options(warn=0)
		fit.y <- fit$fit
		  if(prob == TRUE) {
            fit.y[fit.y < 0] <- 0
		    fit.y[fit.y > 1] <- 1
          }		  
	      if(conf == TRUE) {	
            lower <- fit$fit - stats::qt(0.975, fit$df) * fit$se 
            upper <- fit$fit + stats::qt(0.975, fit$df) * fit$se
          }		  
	  } else {
	    warning("Not a valid option for smoothing type, options are: spline or loess")
      }
		if(conf == FALSE) {	  
	      graphics::plot(x.pt, fit.y, type = "l", xlab=xlab, lwd=0.75, lty=2, ylab=ylab, 
		                 main=main)
		} else {
          graphics::plot(x.pt, y.pt, type = "n", xlab=xlab, ylab=ylab, main=main)
	        graphics::polygon(c(x.pt, rev(x.pt)), c(upper, rev(lower)),  
	                          col=grDevices::rgb(0.85, 0.85, 0.85, 0.5))
		    graphics::lines(x.pt, fit.y, lty=3, col="gray42")
		}				  
      } else {
	    fit.y <- y.pt
        graphics::plot(x.pt, y.pt, type = "l", xlab=xlab, ylab=ylab, main=main)
	  }  
    if( raw.line == TRUE) { graphics::lines(x.pt, y.pt, lty = 1, lwd=0.5, col="red") }
      if( pts == TRUE) { graphics::points(x.pt, y.pt, pch=19, cex=0.35, col="black") }
        if (rug == TRUE) {
          if (n.pt > 10) {
            graphics::rug(stats::quantile(xv, seq(0.1, 0.9, by=0.1)), side = 1)
              } else {
            graphics::rug(unique(xv, side = 1))
            }
	      }  
       }
    }	
  invisible(list(x=x.pt, y=fit.y))
}
