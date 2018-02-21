#' @title Combine Random Forests Ensembles
#' @description Combine two more more random forests models into a single ensemble.
#'
#' @param ...  two or more randomForest class objects
#'
#' @return An object of class randomForest 
#'
#' @note The confusion, err.rate, mse and rsq components (as well as the corresponding components in the test component, if exist) are averaged across ensembles
#' @note This is a modification of the randomForest \code{\link[randomForest]{combine}} function that returns averaged validation statistics
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @examples
#' library(randomForest)
#' data(iris)
#' 
#' c1 <- randomForest(Species ~ ., iris, ntree=50, norm.votes=FALSE)
#' c2 <- randomForest(Species ~ ., iris, ntree=50, norm.votes=FALSE)
#' c3 <- randomForest(Species ~ ., iris, ntree=50, norm.votes=FALSE)
#' 
#' ( class.combine <- rf.combine(c1,c2,c3) )
#' 
#' data(airquality)
#' set.seed(131)
#' r1 <- randomForest(Ozone ~ ., data=airquality, mtry=3,
#'                    importance=TRUE, na.action=na.omit)
#' r2 <- randomForest(Ozone ~ ., data=airquality, mtry=3,
#'                    importance=TRUE, na.action=na.omit)
#' r3 <- randomForest(Ozone ~ ., data=airquality, mtry=3,
#'                    importance=TRUE, na.action=na.omit)
#' 
#' ( regress.combine <- rf.combine(r1,r2,r3) )				   
#' 
#' @seealso \code{\link[randomForest]{randomForest}} for randomForest details
#' @seealso \code{\link[randomForest]{combine}} for original combine function details
#'
#' @exportClass rf.ensembles
#' @export
rf.combine <- function(...) {
   pad0 <- function(x, len) c(x, rep(0, len-length(x)))
   padm0 <- function(x, len) rbind(x, matrix(0, nrow=len-nrow(x),ncol=ncol(x)))
   rflist <- list(...)
     areForest <- sapply(rflist, function(x) inherits(x, "randomForest")) 
   if (any(!areForest)) stop("Argument must be a list of randomForest objects")
     rf <- rflist[[1]]
       classRF <- rf$type == "classification"
         trees <- sapply(rflist, function(x) x$ntree)
           ntree <- sum(trees)
         rf$ntree <- ntree
       nforest <- length(rflist)
     haveTest <- ! any(sapply(rflist, function(x) is.null(x$test)))
   vlist <- lapply(rflist, function(x) rownames(randomForest::importance(x)))
     numvars <- sapply(vlist, length)
   if (! all(numvars[1] == numvars[-1]))
       stop("Unequal number of predictor variables in the randomForest objects.")
    for (i in seq_along(vlist)) {
        if (! all(vlist[[i]] == vlist[[1]]))
            stop("Predictor variables are different in the randomForest objects.")
    }	
  haveForest <- sapply(rflist, function(x) !is.null(x$forest))
    if (all(haveForest)) {
      nrnodes <- max(sapply(rflist, function(x) x$forest$nrnodes))
      rf$forest$nrnodes <- nrnodes
      rf$forest$ndbigtree <- unlist(sapply(rflist, function(x) x$forest$ndbigtree))
      rf$forest$nodestatus <- do.call("cbind", lapply(rflist, function(x)
                                      padm0(x$forest$nodestatus, nrnodes)))
      rf$forest $bestvar <- do.call("cbind", lapply(rflist, function(x)
                                    padm0(x$forest$bestvar, nrnodes)))
      rf$forest$xbestsplit <- do.call("cbind", lapply(rflist, function(x)
                                      padm0(x$forest$xbestsplit, nrnodes)))
      rf$forest$nodepred <- do.call("cbind", lapply(rflist, function(x)
                                    padm0(x$forest$nodepred, nrnodes)))
      tree.dim <- dim(rf$forest$treemap)
      if (classRF) {
        rf$forest$treemap <- array(unlist(lapply(rflist, function(x) 
	                               apply(x$forest$treemap, 2:3, pad0, nrnodes))),
                                         c(nrnodes, 2, ntree))
      } else {
        rf$forest$leftDaughter <- do.call("cbind", lapply(rflist, function(x)
                                          padm0(x$forest$leftDaughter, nrnodes)))
        rf$forest$rightDaughter <- do.call("cbind", lapply(rflist, function(x)
                                           padm0(x$forest$rightDaughter, nrnodes)))
      }
    rf$forest$ntree <- ntree
      if (classRF) rf$forest$cutoff <- rflist[[1]]$forest$cutoff
    } else {
      rf$forest <- NULL
    }	
   if (classRF) {
       rf$votes <- 0
       rf$oob.times <- 0
       areVotes <- all(sapply(rflist, function(x) any(x$votes > 1, na.rf=TRUE)))
       if (areVotes) {
         for(i in 1:nforest) {
           rf$oob.times <- rf$oob.times + rflist[[i]]$oob.times
           rf$votes <- rf$votes +
               ifelse(is.na(rflist[[i]]$votes), 0, rflist[[i]]$votes)
         }
       } else {
         for(i in 1:nforest) {
           rf$oob.times <- rf$oob.times + rflist[[i]]$oob.times            
           rf$votes <- rf$votes +
           ifelse(is.na(rflist[[i]]$votes), 0, rflist[[i]]$votes) * rflist[[i]]$oob.times
           }
           rf$votes <- rf$votes / rf$oob.times
       }
       rf$predicted <- factor(colnames(rf$votes)[max.col(rf$votes)],
                              levels=levels(rf$predicted))
       if(haveTest) {
         rf$test$votes <- 0
         if (any(rf$test$votes > 1)) {
           for(i in 1:nforest) rf$test$votes <- rf$test$votes + rflist[[i]]$test$votes
         } else {
           for (i in 1:nforest)
             rf$test$votes <- rf$test$votes + rflist[[i]]$test$votes * rflist[[i]]$ntree
         }
         rf$test$predicted <- factor(colnames(rf$test$votes)[max.col(rf$test$votes)],
                                     levels=levels(rf$test$predicted))
       }
    } else {
      rf$predicted <- 0
      for (i in 1:nforest) rf$predicted <- rf$predicted + rflist[[i]]$predicted * rflist[[i]]$ntree
       rf$predicted <- rf$predicted / ntree
       if (haveTest) {
         rf$test$predicted <- 0
           for (i in 1:nforest) rf$test$predicted <- rf$test$predicted + rflist[[i]]$test$predicted * rflist[[i]]$ntree
         rf$test$predicted <- rf$test$predicted / ntree
      }
    }
   have.imp <- !any(sapply(rflist, function(x) is.null(x$importance)))
    if (have.imp) {
      rf$importance <- rf$importanceSD <- 0
      for(i in 1:nforest) {
        rf$importance <- rf$importance + rflist[[i]]$importance * rflist[[i]]$ntree
        rf$importanceSD <- rf$importanceSD + rflist[[i]]$importanceSD^2 * rflist[[i]]$ntree
      }
    rf$importance <- rf$importance / ntree
    rf$importanceSD <- sqrt(rf$importanceSD / ntree)
    haveCaseImp <- !any(sapply(rflist, function(x) is.null(x$localImportance)))
      if (haveCaseImp) {
        rf$localImportance <- 0
        for(i in 1:nforest) {
          rf$localImportance <- rf$localImportance + rflist[[i]]$localImportance * rflist[[i]]$ntree
        }
        rf$localImportance <- rf$localImportance / ntree
      }
    }
   have.prox <- !any(sapply(rflist, function(x) is.null(x$proximity)))
   if (have.prox) {
     rf$proximity <- 0
       for(i in 1:nforest)
         rf$proximity <- rf$proximity + rflist[[i]]$proximity * rflist[[i]]$ntree
       rf$proximity <- rf$proximity / ntree
    }
	hasInBag <- all(sapply(rflist, function(x) !is.null(x$inbag)))
   	  if (hasInBag) rf$inbag <- do.call(cbind, lapply(rflist, "[[", "inbag"))
	  
   	if (classRF) {
	  confusion <- lapply(rflist, function(x) x$confusion)
        rf$confusion <- round(Reduce('+', confusion) / length(confusion), 3)
        rf$err.rate <- stats::median(unlist(lapply(rflist, function(x) x$err.rate[,1])))
      if (haveTest) {
	    tconfusion <- lapply(rflist, function(x) x$test$confusion)	  
        rf$test$confusion <- round(Reduce('+', tconfusion) / length(tconfusion), 3)
        rf$test$err.rate <- stats::median(unlist(lapply(rflist, function(x) x$test$err.rate[,1])))  
      }
   	} else {
	  rsq <- unlist(lapply(rflist, function(x) x$rsq))
	    rf$rsq <- round(100*rsq[length(rsq)], digits=2) 
	  mse <- unlist(lapply(rflist, function(x) x$mse))
        rf$mse <- mse[length(mse)] 
	  if (haveTest) {
        rsq <- unlist(lapply(rflist, function(x) x$test$rsq))
		  rf$test$rsq <- round(100*rsq[length(rsq)], digits=2) 
        mse <- unlist(lapply(rflist, function(x) x$test$mse))
          rf$test$mse <- mse[length(mse)] 
	  } else {
        rf$test$rsq <- NULL
		rf$test$mse <- NULL
      }	  
   	}
	  rf$nrf <- length(rflist)
	class(rf) <- c("rf.ensembles", class(rf))
  return( rf )
}
