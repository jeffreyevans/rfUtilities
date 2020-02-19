#' @title Plot random forests cross-validation
#' @description Plot function for rf.cv object 
#'
#' @param  x        A rf.cv object
#' @param  stat     Which statistic to plot: classification: "users.accuracy", "producers.accuracy", 
#'                    "kappa", "oob", regression: "rmse", "mse", "var.exp", "mae", "mbe"  
#' @param  ...      Additional arguments passed to plot
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @method plot rf.cv 
#'
#' @export    	     
plot.rf.cv <- function(x, stat = "kappa", ...) {
  if(class(x)[2] == "classification") {
  plot.class <- function(x, ...) {
    dots <- as.list(match.call(expand.dots = TRUE)[-1])
    dots[["x"]] <- 1:nrow(x)
    dots[["y"]] <- sort(x[,1])
    dots[["type"]] <- "n"
    dots[["ylim"]] <- c(min(x),max(x)) 
  	if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <-  "index"
  	if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <-  ""
  	if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <- stat
  	  do.call("plot", dots)
  	    graphics::axis(side=1, at = as.numeric(as.factor(rownames(x$effect.size))), 
  	                   labels = rownames(x$effect.size))
          for(i in 1:ncol(x)) { graphics::lines(dots[["x"]], sort(x[,i]), col=i) }
  	        graphics::legend("bottomright", legend=colnames(x), col=1:nrow(x), 
  	                         lty = rep(1,nrow(x)), bg="white")   
  }   
  plot.kappa <- function(x, ...) {
    dots <- as.list(match.call(expand.dots = TRUE)[-1])
    dots[["x"]] <- stats::smooth.spline(1:nrow(x), sort(x[,"CV.kappa"]))$x
    dots[["y"]] <- stats::smooth.spline(1:nrow(x), sort(x[,"CV.kappa"]))$y
    dots[["type"]] <- "n"
  	if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <-  ""
  	if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <-  ""
  	if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <- "Kappa"
  	  do.call("plot", dots)
        graphics::lines(dots[["x"]], dots[["y"]]) 
  }
  plot.oob <- function(x, ...) {
    dots <- as.list(match.call(expand.dots = TRUE)[-1])
    dots[["x"]] <- 1:nrow(x)
    dots[["y"]] <- sort(x[,"CV.PCC"])
    dots[["type"]] <- "n"
  	if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <-  ""
  	if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <-  ""
  	if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <- "PCC"
  	  do.call("plot", dots)
        for(i in 1:(ncol(x)-1)) { 
		  graphics::lines(1:nrow(x), sort(x[,i]), col=i)
        }
      graphics::legend("bottomright", legend = colnames(x)[1:(ncol(x)-1)],  
  	                   col=1:nrow(x), lty = rep(1,ncol(x)-1), bg="white")   		
  }
    if(type == "cv" & stat == "users.accuracy") 
      { dat <- x$cv.users.accuracy 
    } else if(type == "cv" & stat == "producers.accuracy")
      {  dat <- x$cv.producers.accuracy 
    } else if(stat == "kappa" | stat == "oob")
      {  dat <- x$cv.oob } 	
  
    if( stat == "users.accuracy" | stat == "producers.accuracy" ) {
      plot.class( dat, ...)
    } else if( stat == "oob") {
      plot.oob(dat, ...)
    } else if( stat == "kappa") {
      plot.kappa(dat, ...)
    }
  } 
  if(class(x)[2] == "regression") {
   if( stat != "rmse" & stat != "mse" & stat != "var.exp" & stat != "mbe" & stat != "mae") 
     stat = "rmse"   
    if(stat == "rmse") 
        { dat <- x[["y.rmse"]]
		  slab = "Cross-validated Root Mean Squared Error"
      } else if(stat == "mse")
        {  dat <- x[["model.mse"]]
           slab = "Model Mean Square Error"
           fit <- x[["fit.mse"]] 		   
     }  else if(stat == "var.exp")
        {  dat <- dat <- x[["model.varExp"]]
           slab = "Model percent variance explained"
		   fit <- x[["fit.var.exp"]] 
     }  else if(stat == "mbe")
        {  dat <- x[["y.mbe"]]
		   slab = "Cross-validated Mean Bias Error"
	 }  else if(stat == "mae")
        {  dat <- x[["y.mae"]]
		   slab = "Cross-validated Mean Absolute Error"
		}  else if(stat == "ks")
        {  dat <- x[["D"]]
		   slab = "Model Kolmogorov-Smirnov statistic"
		}   
    plot.reg <- function(x, s = stat, ...) {
      dots <- as.list(match.call(expand.dots = TRUE)[-1])
        dots[["x"]] <- stats::smooth.spline(1:length(x), sort(x))$x
        dots[["y"]] <- stats::smooth.spline(1:length(x), sort(x))$y
        dots[["type"]] <- "n"
    	  if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <-  ""
    	  if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <-  ""
    	  if (is.null(dots[["main"]]) & "main" %in% names(dots) == FALSE) dots[["main"]] <- slab 
        do.call("plot", dots)
      graphics::lines(dots[["x"]], dots[["y"]])
        if(stat == "mse" | stat == "var.exp") graphics::abline(h=fit, col="blue")
    }
    plot.reg(dat, s = stat, ...)
  } 
}
