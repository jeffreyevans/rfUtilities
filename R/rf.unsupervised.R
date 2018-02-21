#' @title Unsupervised Random Forests
#' @description Performs an unsupervised Random Forests for returning clustering, based on dissimilarity, and optional neighbor distance. 
#'
#' @param x                A matrix/data/frame object to cluster
#' @param n                Number of clusters
#' @param proximity        (FALSE/TRUE) Return matrix of neighbor distances based on proximity 
#' @param silhouettes      (FALSE/TRUE) Return adjusted silhouette values 
#' @param clara            (FALSE/TRUE) Use clara partitioning, for large data
#' @param ...              Additional Random Forests arguments 
#'
#' @return A vector of clusters or list class object  of class "unsupervised", containing the following components:
#' @return  distances              Scaled proximity matrix representing dissimilarity neighbor distances  
#' @return  k                      Vector of cluster labels using adjusted silhouettes  
#' @return  silhouette.values      Adjusted silhouette cluster labels and silhouette values 
#'
#' @note Clusters (k) are derived using the random forests proximity matrix, treating it as dissimilarity neighbor distances. 
#' @note The clusters are identified using a Partitioning Around Medoids where negative silhouette values are assigned to the nearest neighbor. 
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references Rand, W.M. (1971) Objective Criteria for the Evaluation of Clustering Methods. Journal of the American Statistical Association, 66:846-850.
#' @references Shi, T., Seligson, D., Belldegrun, A.S., Palotie, A., and Horvath, Ss (2005) Tumor Classification by Tissue Microarray Profiling: Random Forest Clustering Applied to Renal Cell Carcinoma. Modern Pathology, 18:547-557.  
#' 
#' @examples 
#'  library(randomForest) 
#'  data(iris)
#'  n = 4
#'  clust.iris <- rf.unsupervised(iris[,1:4], n=n, proximity = TRUE, 
#'                                silhouettes = TRUE)
#'  clust.iris$k
#'
#'  mds <- stats:::cmdscale(clust.iris$distances, eig=TRUE, k=n)
#'    colnames(mds$points) <- paste("Dim", 1:n)
#'    mds.col <- ifelse(clust.iris$k == 1, rainbow(4)[1],
#'                 ifelse(clust.iris$k == 2, rainbow(4)[2],
#'  			     ifelse(clust.iris$k == 3, rainbow(4)[3],
#'  				   ifelse(clust.iris$k == 4, rainbow(4)[4], NA))))
#'  plot(mds$points[,1:2],col=mds.col, pch=20) 				   
#'  pairs(mds$points, col=mds.col, pch=20)
#'   
#' @seealso \code{\link[randomForest]{randomForest}} for ... options
#' @seealso \code{\link[cluster]{pam}} for details on Partitioning Around Medoids (PAM)  
#' @seealso \code{\link[cluster]{clara}} for details on Clustering Large Applications (clara) 
#' 
#' @exportClass unsupervised 
#' @export
rf.unsupervised <- function(x, n = 2, proximity = FALSE, silhouettes = FALSE, 
                            clara = FALSE, ...) {
  if( !class(x) == "matrix" & !class(x) == "data.frame") 
      stop("x must be data.frame or matrix object") 
  silhouette <- function(x, k, d = inherits(x, "dist"), large = FALSE, ...) {
      if (d) {
        if (!is.null(attr(x, "Labels"))) { original.row.names <- attr(x, "Labels")}
          names(x) <- as.character(c(1:attr(x, "Size")))
      } else {
        if(!is.null(dimnames(x)[[1]])) { original.row.names <- dimnames(x)[[1]]}
          row.names(x) <- as.character(c(1:dim(x)[[1]]))
      }
	  if(large == TRUE) {
	    message("Running Clustering Large Applications version of PAM")
	    pam1 <- cluster::pam(x, k, diss = d, ...)
	  } else {
        pam1 <- cluster::pam(x, k, diss = d, ...)
	  }
       label2 <- pam1$clustering
        silinfo1 <- pam1$silinfo$widths
         index1 <- as.numeric(as.character(row.names(silinfo1)))
         silinfo2 <- silinfo1[order(index1),]
       labelnew <- ifelse(silinfo2[,3]<0, silinfo2[,2], silinfo2[,1])
       names(labelnew) <- original.row.names
     return( list(k = labelnew, sil = silinfo2) )    
    }
  
  unlabled.rf <- randomForest::randomForest(x, ...)  
  krf <- silhouette( (1 - unlabled.rf$proximity), k = n, large = clara)
  k <- krf$k  
    if(proximity == TRUE) k <- list(distances = (1 - unlabled.rf$proximity), k = k)
      if(silhouettes == TRUE) {
	    if(class(k) == "numeric") {
		  k <- list(krf$sil, k = k) 
		    } else {
	      k[["silhouette.values"]] <- krf$sil
        }
      }
    class(k) <- "unsupervised"
  return(k)
}
