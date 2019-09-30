#' @title Accuracy
#' @description Classification accuracy measures for pcc, kappa, users accuracy, producers accuracy
#'    
#' @param x   vector of predicted data or table/matrix contingency table
#' @param y   vector of observed data, if x is not table/matrix contingency table 
#'
#' @return A list class object with the following components:
#' \itemize{ 
#' \item   PCC                    percent correctly classified (accuracy)
#' \item   auc                    Area Under the ROC Curve
#' \item   users.accuracy         The users accuracy  
#' \item   producers.accuracy     The producers accuracy
#' \item   kappa                  Cohen's Kappa (chance corrected accuracy)
#' \item   true.skill             Hanssen-Kuiper skill score (aka true score statistic)
#' \item   sensitivity            Sensitivity (aka, recall)
#' \item   specificity            Specificity 
#' \item   plr                    Positive Likelihood Ratio   
#' \item   nlr                    Negative Likelihood Ratio  
#' \item   typeI.error            Type I error (omission)
#' \item   typeII.error           Type II error (commission)
#' \item   gini                   Gini entropy index
#' \item   f.score                F-score
#' \item   gain                   Information gain (aka precision)
#' \item   mcc                    Matthew's correlation 
#' \item   confusion              A confusion matrix 
#'  }
#' 
#' @note
#'   \itemize{ 
#'   \item   sensitivity = true positives / ( true positives + false positives ) 
#'   \item   specificity = true negatives / ( true negatives + false positives ) 
#'   \item   Type I error = 1 - specificity
#'   \item   Type II error = 1 - sensitivity
#'   \item   Positive Likelihood Ratio  = sensitivity / (1 - specificity) 
#'   \item   Negative Likelihood Ratio  = (1 - sensitivity) / specificity
#'   \item   gain  = sensitivity / ( (true positives + true negatives) / n )
#'   \item   auc = (tpr - fpr + 1) / 2
#'   \item   F-Score = 2 * (precision * recall) / (precision + recall) 
#'   \item   Hanssen-Kuiper skill score (aka true score statistic) = [(tp * tn) - (fp * fn)] / [(tp + fn) + (fp + tn)], The true skill score has an expected -1 to +1, with 0 representing no discrimination.   
#'  } 
#' @note Using the table function matrix positions for a 2x2 confusion matrix are TP(1), FN(3), FP(2), TN(4)
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references
#' Cohen, J. (1960) A coefficient of agreement for nominal scales. Educational and Psychological Measurement 20 (1):37-46
#' Cohen, J. (1968) Weighted kappa: Nominal scale agreement with provision for scaled disagreement or partial credit. Psychological Bulletin 70 (4):213-220  
#' Powers, D.M.W., (2011). Evaluation: From Precision, Recall and F-Measure to ROC, Informedness, Markedness & Correlation. Journal of Machine Learning Technologies 2(1):37-63.
#'
#' @examples 
#'  # Two classes (vector)
#'  observed <- sample(c(rep("Pres",50),rep("Abs",50)), 100, replace=TRUE )
#'  accuracy(observed[sample(1:length(observed))], observed)
#'
#'  # Two classes (contingency table)
#' accuracy(cbind(c(15,11), c(2,123)))
#'
#'  # Multiple classes
#'  accuracy(iris[sample(1:150),]$Species, iris$Species)
#'
#' @export 
accuracy <- function (x, y) {
  if(inherits(x, c("table", "matrix"))) {
    t.xy <- x 
  } else {
    t.xy <- table(x, y)
  }
    if(inherits(t.xy, "matrix")) {
      if(is.null(rownames(t.xy))) {	
        rownames(t.xy) <- c("0","1")
	    colnames(t.xy) <- c("0","1")
	  }
    }
	tc <- match(colnames(t.xy), rownames(t.xy))
	mtc <- matrix(ncol = ncol(t.xy), nrow = length(tc[tc == "NA"]), 0)
        nrn <- colnames(t.xy)[is.na(tc) == TRUE]
          rownames(mtc) <- nrn
            t1 <- rbind(t.xy, mtc)
              tr <- match(rownames(t1), colnames(t1))
                mtr <- matrix(nrow = nrow(t1), ncol = length(tr[tr == "NA"]), 0)
                  ncn <- rownames(t1)[is.na(tr) == TRUE]
                colnames(mtr) <- ncn
              t2 <- cbind(t1, mtr)
            sr <- sort(rownames(t2))
          mr <- match(sr, rownames(t2))
        t3 <- t(t2[mr, ])
      sc <- sort(rownames(t3))
    mc <- match(sc, rownames(t3))
    t4 <- t(t3[mc, ])
      agree <- diag(t4)
        prod1 <- apply(t4, 1, sum)
          prod2 <- agree / prod1
             user1 <- apply(t4, 2, sum)
            user2 <- agree / user1
          N <- sum(t4)
        k1 <- sum(agree)
      k2 <- sum(prod1 * user1)
    khat <- abs(((N * k1) - k2) / (N^2 - k2))	
    if( dim(t.xy)[1] == 2 ) {
      # TP(1)  FN(3)  
      # FP(2)  TN(4)  	
      n = sum(t.xy)  # N	
      TP <- t.xy[1]  # True Positives, Power 
	  FP <- t.xy[2]  # False Positives, Type-I error 
	  FN <- t.xy[3]  # False Negatives, Type-II error 
	  TN <- t.xy[4]  # True Negatives
      # prevalence <- TP / n
      precision <- TP / (TP + FP)  	  
	  tpr <- TP / (TP + FN)  # true positive rate (aka, sensitivity, recall)
	  tnr <- TN / (TN + FP)  # true negative rate (aka, specificity, selectivity)	 
	  fpr <- FP / (FP + TN) 
	  fnr <- FN / (FN + TP)  # Beta       
	  type1.error <- 1 - tnr                     
      type2.error <- 1 - tpr                     
      plr <- tpr / (1 - tnr)                     
      nlr <- (1 - tpr) / tnr 
	  auc <- (tpr - fpr + 1) / 2
	  gini <- 2 * auc - 1  
	  f.score <- 2 * (precision * tpr) / (precision + tpr)
	  true.skill <- ( (t.xy[1] * t.xy[4]) - (t.xy[3] * t.xy[2]) ) / 
	                ( (t.xy[1] + t.xy[2]) * (t.xy[3] + t.xy[4]) )
	  gain <- precision / ( (t.xy[1] + t.xy[4]) / n )
	  mcc <- (TP * TN - FP * FN) / sqrt( (TP + FP) * (TP + FN) * 
			 (TN + FP) * (TN + FN) )			 
		confusion <- matrix(c(paste0("True positive(", TP, ")"),
		         paste0("False positive(", FP, ")"),
		         paste0("False negative(", FN, ")"),
		         paste0("True negative(", TN, ")")),
				 nrow=2, byrow=TRUE)
          rownames(confusion) <- rownames(t.xy)
		  colnames(confusion) <- colnames(t.xy)
	  acc <- list( PCC = (sum(diag(t4))/sum(t4)) * 100,
                   auc = auc, 	
	               users.accuracy = round(user2 * 100, 1),  
	               producers.accuracy = round(prod2 * 100, 1),
                   kappa = round(khat, 4),
				   true.skill = true.skill, 
				   sensitivity = tpr,
				   specificity = tnr,
				   plr = plr,
				   nlr = nlr,
				   typeI.error = type1.error,
				   typeII.error = type2.error,
				   gini = gini,
				   f.score = f.score,
				   gain = gain,
                   matthews = mcc,
				   confusion = confusion )	
    } else {
	  acc <- list( PCC = (sum(diag(t4))/sum(t4)) * 100, 
	             users.accuracy = round(user2 * 100, 1),  
	             producers.accuracy = round(prod2 * 100, 1),
                 kappa = round(khat, 4),
                 confusion = t.xy )
	}			 
	class(acc) <- c("accuracy", "list") 			   
          return( acc )
}
