#' @title Bivariate partial-dependency plot
#' @description Bivariate partial dependence provides a graphical depiction of the marginal effect of two variables on the class probability (classification) or response (regression)
#'
#' @param x             random forest object
#' @param pred.data     data.frame of independent variables used in model
#' @param v1            Variable 1 used in partial dependency
#' @param v2            Variable 2 used in partial dependency
#' @param grid.size     Number of grid cells (NxN) to integrate partial dependency for 
#' @param which.class   Index of class probability (only if classification) 
#' @param plot          (TRUE/FALSE) Plot 3D surface
#' @param col.ramp      Colors used in building color ramp
#' @param ncols         Number of colors in color ramp
#' @param ...           Arguments passed to persp
#'
#' @return A list object with vectors of v1 (p1) and v2 (p2) and a matrix (estimate), estimate of the averaged estimates.
#'
#' @note In deriving the partial-dependence, at each plotted point, the background variables are held at their median values
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references Friedman, J.H. (2001) Greedy Function Approximation:  A Gradient Boosting Machine. Annals of Statistics 19(1)
#' @references Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer 
#' @references Baruch-Mordo, S., J.S. Evans, J. Severson, J. D. Naugle, J. Kiesecker, J. Maestas, & M.J. Falkowski (2013) Saving sage-grouse from the trees: A proactive solution to reducing a key threat to a candidate species Biological Conservation 167:233-241  
#'
#' @examples
#'  library(randomForest)
#'    data(iris)
#'    iris$Species <- ifelse( iris$Species == "versicolor", 1, 0 ) 
#'     
#'  # Add some noise
#'  idx1 <- which(iris$Species %in% 1)
#'  idx0 <- which( iris$Species %in% 0)
#'    iris$Species[sample(idx1, 2)] <- 0
#'    iris$Species[sample(idx0, 2)] <- 1
#'   
#'  # Specify model  
#'  y = iris[,"Species"] 
#'  x = iris[,1:4]
#'  
#'  set.seed(4364)  
#'  ( rf.mdl1 <- randomForest(x=x, y=factor(y)) )
#'  
#'  
#'  ( bvpd <- bivariate.partialDependence(rf.mdl1, iris, 
#'                    v1 = "Petal.Length", v2 = "Petal.Width", shade = 0.6,
#'                    grid.size = 20, ncols=100, border=NA, col.ramp=c("green","blue") ) ) 
#'		
#' @seealso \code{\link[graphics]{persp}} for persp ... plotting options
#'
#' @export
bivariate.partialDependence <- function(x, pred.data, v1, v2, grid.size = 20, which.class = 2, 
                                        plot = TRUE, col.ramp = c("#ffffff", "#2a2a2a"), 
										ncols = 20, ...) {
  if(!any(class(x) %in% c("randomForest","list"))) stop("x is not a randomForest object")
    	dots <- as.list(match.call(expand.dots = TRUE)[-1])
    s1 <- seq(from = min(pred.data[,v1]), to = max(pred.data[,v1]),
               by = (max(pred.data[,v1]) - min(pred.data[,v1]))/(grid.size-1))
    s2 <- seq(from = min(pred.data[,v2]), to = max(pred.data[,v2]),
                     by = (max(pred.data[,v2]) - min(pred.data[,v2]))/(grid.size-1))
      v <- expand.grid(s1, s2)
      v <- v[with(v, order(Var1, Var2)),]
      vrep <- pred.data[rep(1:nrow(pred.data), nrow(v)),]
        vrep[,v1] <- rep(v$Var1, each = nrow(pred.data))
        vrep[,v2] <- rep(v$Var2, each = nrow(pred.data))	
    if(x$type == "classification") {	
      vrep$pred <- stats::predict(x, vrep[,which(names(vrep) %in% rownames(x$importance))], 
	                              type="prob")[,which.class]
    } else {
      vrep$pred <- stats::predict(x, vrep[,which(names(vrep) %in% rownames(x$importance))])
    }  
    idx <- sort(rep(1:nrow(v), length.out=nrow(vrep)))   
      idx.med <- as.numeric(tapply(vrep$pred, idx, FUN=stats::median)) 
    z <- matrix(idx.med, nrow = length(s1), byrow = TRUE) 
	if(plot == TRUE) { 
      dots[["x"]] <- s1
	  dots[["y"]] <- s2
	  dots[["z"]] <- z
    if (is.null(dots[["xlab"]]) & "xlab" %in% names(dots) == FALSE) dots[["xlab"]] <- deparse(substitute(v1))
      if (is.null(dots[["ylab"]]) & "ylab" %in% names(dots) == FALSE) dots[["ylab"]] <- deparse(substitute(v2))
  	    if (is.null(dots[["zlab"]]) & "zlab" %in% names(dots) == FALSE) dots[["zlab"]] <- "\nPredicted Value"
	    if (is.null(dots[["theta"]]) & "theta" %in% names(dots) == FALSE) dots[["theta"]] <- -45
      if (is.null(dots[["cex.lab"]]) & "cex.lab" %in% names(dots) == FALSE) dots[["cex.lab"]] <- 1
	    if (is.null(dots[["col"]]) & "col" %in% names(dots) == FALSE) {
          jet.colors <- grDevices::colorRampPalette( col.ramp ) 
            color <- jet.colors(ncols)
            zfacet <- z[-1, -1] + z[-1, -1 * length(s1)] + 
                      z[-1 * length(s2), -1] + z[-1 * length(s1), -1 * length(s2)]
          facetcol <- cut(zfacet, ncols)
          dots[["col"]] <- color[facetcol]
	    }
	  options(warn=-1)
	    do.call("persp", dots)
	  options(warn=0)
    }	
  return( invisible(list( p1 = s1, p2 = s2, estimate = z )) )
}
