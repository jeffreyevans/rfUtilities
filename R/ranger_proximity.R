#' @title random forest ranger proximity matrix
#' @description Calculates OOB proximities from ranger object
#'
#' @param y       A range object with keep.inbag = TRUE 
#' @param x       Independent data used in model
#' @param sparse  (FALSE/TRUE) Output a sparse matrix 
#'
#' @return An n x n matrix where position i, j gives the proportion 
#'         of times observation i and j are in the same terminal 
#'         node across all trees. If sparse = TRUE, object is a 
#'         sparse matrix dsCMatrix class from the Matrix package 
#'
#' @author Jeffrey S. Evans  <jeffrey_evans<at>tnc.org
#'
#' @examples
#' library(ranger)
#' fit <- ranger(Species ~ ., iris, keep.inbag = TRUE)
#' p <- ranger_proximity(fit, iris[,-5]) 
#' 
#' @export ranger_proximity 
ranger_proximity <- function(y, x, sparse = FALSE) {
  if (class(y) != "ranger")
    stop("y is not a ranger object")
  if (is.null(y$inbag.counts))
    stop("call ranger with keep.inbag = TRUE")
  pred <- ranger:::predict.ranger(y, x, type = "terminalNodes")$predictions
    prox <- matrix(NA, nrow(pred), nrow(pred))
      ntree = ncol(pred)
    n = nrow(prox)   
  inbag = simplify2array(y$inbag.counts)
    prox <- outer(1:n, 1:n,
      Vectorize(function(i,j) {
        tree_idx <- inbag[i, ] == 0 & inbag[j, ] == 0
        prox[i, j] <- sum(pred[i, tree_idx] == pred[j, tree_idx]) / sum(tree_idx)
      }))
	if(sparse){ 
	  prox <- Matrix::Matrix(prox, sparse = TRUE)
	  message("Output is a sparse Matrix::dsCMatrix object")
    } 	
  return(prox)
}
