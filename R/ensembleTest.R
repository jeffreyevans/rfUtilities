#' Tests the correlation across a Bootstrap ensemble to ensure
#' a quasi-independent sample and uncorrelated ensemble. 
#'
#' @parm x            A matrix or data.frame of the data to test
#' @parm perm         Number of Bootstraps (eg., same as random forests ntree)        
#' @parm n            Sample proportion
#' @parm p            Accept/Reject p-value
#' @parm type         Type of matrix similarity statistic
#' @parm replacement  (TRUE/FALSE) Sample with replacement 
#'
#' @return 
#'
#' @note 
#' type options are:
#'   \itemize{ 
#'   \item mcv - multivariate correlation (RV-coefficient) 
#'   \item cov - covariance equivalence
#'   \item psi - Procrustes Similarity Index 
#'   \item pairwise - Averaged column pairwise comparison  
#'   }
#'
#' @author Jeffrey S. Evans   <jeffrey_evans<at>tnc.org>
#'
#' @references 
#' Robert, P., Escoufier, Y. (1976). A Unifying Tool for Linear Multivariate 
#'   Statistical Methods: The RV-Coefficient. Applied Statistics 25(3):257-265.
#'
#' Smilde, AK., Kiers, HA., Bijlsma, S., Rubingh, CM., van Erk, MJ (2009)
#'   Matrix correlations for high-dimensional data: the modified RV-coefficient. 
#'   Bioinformatics 25(3): 401-5.
#' 
#' library(raster)
#' r <- stack("C:/evans/India/Chennai/data/2019_OLI8.tif")
#'   x <- sampleRandom(r, 5000)
#'  
#'  
#' @export ensembleTest  
ensembleTest <- function(x, perm = 99, p = 0.05, sampling = c("boot", "elimination"), 
                         type = c("mvc", "cov", "psi", "pairwise"), 
						 n = NULL, replacement = TRUE) {		 
	    if(!any(class(x)[1] == c("matrix", "data.frame")))
		  stop(deparse(substitute(x)), "must be a matrix or data.frame")
        if(ncol(x) < 3)
          stop("Must have more than two paramters")		
		if(any(unique(sapply(x, is.factor)) == TRUE))
		  stop("Cannot have factors in", deparse(substitute(x)))
		if(any(unique(sapply(x, is.character)) == TRUE))
		  stop("Cannot have non-numeric data in", deparse(substitute(x)))  
		if(is.null(n)) {
		  if(sampling == "boot") { n = nrow(x) } else { n = nrow(x) * 0.34}
		} else {
		  n = round(nrow(x) * n, 0)
		}	
    cov.equivalence <- function(m1, m2, pVal, verbose = FALSE) {
       k = 2
        p = 2
         n1 = nrow(m1)
          n2 = nrow(m2) 
           n = n1 + n2
            s1 <- crossprod(m1[1:nrow(m1)])
             s2 <- crossprod(m2[1:nrow(m2)])
              c1 = (1/(n1-1)) * s1
              c2 = (1/(n2-1)) * s2
             c3 = (s1+s2)/(n-k)
            d = det(c3)
            d1 = det(c1)
           d2 = det(c2) 
          m = ( (n - k) * log(d) ) - ( (n1 - 1) * log(d1) + (n2 - 1) * log(d2) )
         h = 1 - ((2 * p * p + 3 * p - 1) / (6 * (p + 1) * (k - 1)) * 
		         (1 / (n1 - 1) + 1 / (n2 - 1) + 1 / (n - k)))
        chi = round(abs(m * h),digits=6)
        dfree = p * (p + 1) * (k - 1) / 2
		if(verbose){
		  if(missing(pVal)) pVal = 0.05
          cat("Equivalence p =", chi, "with", dfree, "degrees of freedom", "\n")
            if ( (chi <= pVal ) == TRUE & ( dfree > 2) | (dfree > 20)  == TRUE ) { 
              cat("You can reject the null hypothesis", "\n")
            } else {
              cat("The null hypothesis cannot be rejected", "\n")
            }
		}	
	  return(chi)
    }
  cat("testing ensemble using", sampling, "and", type, "test statistic", "\n") 	
  sample.cor <- vector()
    for(i in 1:perm) {
      if(i == 1) {
	    if(sampling == "boot") {
          m2 <- x[sample(1:nrow(x), n, replace=replacement),]
		} else if(sampling == "elimination") {
		  m2 <- x[-sample(1:nrow(x), n),] 
	    }
      next } 
        if(exists("m2")) m1 <- m2  
	      if(sampling == "boot") {
            m2 <- x[sample(1:nrow(x), n, replace=replacement),]
		  } else if(sampling == "elimination") {
		    m2 <- x[-sample(1:nrow(x), n),] 
		  }
        if(type == "mvc") {
          sample.cor <- append(sample.cor, MatrixCorrelation::RV2(m1, m2))
		} else if(type == "cov") {
		  sample.cor <- append(sample.cor, cov.equivalence(cov(m1), cov(m2)))
		} else if(type == "psi") {		  
          sample.cor <- append(sample.cor, MatrixCorrelation::PSI(m1, m2))
		} else if(type == "pairwise") {
         sample.cor <- append(sample.cor, mean(cancor(m1, m2)$cor))
		}
	  }
  return( sample.cor )   
}
