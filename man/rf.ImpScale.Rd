% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf.ImpScale.R
\name{rf.ImpScale}
\alias{rf.ImpScale}
\title{Scaling of Random Forests importance values}
\usage{
rf.ImpScale(x, scaling = c("mir", "se", "p"), n = 99, sort = FALSE)
}
\arguments{
\item{x}{randomForest or ranger object with a classification or 
regression instance}

\item{scaling}{Type of importance scaling, options are ("mir","se", "p")}

\item{n}{If scaling = "p" how many permutations}

\item{sort}{Sort output by importance}
}
\value{
A data.frame with:
\itemize{  
  \item   {"parameter"} {Name of the parameter}
  \item   {"importance"} {scaled importance values}
  \item   {"pvalue"} {Optional p-value for ranger objects}
}
}
\description{
Various scaling approaches for scaling permuted parameter 
             importance values
}
\details{
The "mir" scale option performs a row standardization and the "se" option 
performs normalization using the "standard errors" of the permutation-based 
importance measure. Both options result in a 0-1 range but, "se" sums to 1.

The scaled importance measures are calculated as: 
  mir = i/max(i) and se = (i / se) / ( sum(i) / se).

The "p" scale option calculates a p-value using the Janitza et al, (2018)
  method where, assuming that noise variables vary randomly around zero,
  uses negative values to build a null distribution.
}
\examples{
library(randomForest)
library(ranger)

#### Classification
data(iris)
iris$Species <- as.factor(iris$Species)

# Random Forests using x, y arguments with index
( rf.class <- randomForest(x = iris[,1:4], y = iris[,"Species"], 
                           importance = TRUE) )

  rf.ImpScale(rf.class, scaling="mir")
  rf.ImpScale(rf.class, scaling="se")

# ranger using formula interface 
( rgr.class <- ranger(Species ~., data=iris, 
                      importance = "permutation") )
				   
  rf.ImpScale(rgr.class, scaling="mir")
  rf.ImpScale(rgr.class, scaling="p", n=10)

#### Regression
data(airquality)
  airquality <- na.omit(airquality)
  
( rf.reg <- randomForest(x=airquality[,2:6], y=airquality[,1],
                         importance = TRUE) )
  rf.ImpScale(rf.reg, scaling="mir")
  rf.ImpScale(rf.reg, scaling="se")

( rgr.reg <- ranger(x=airquality[,2:6], y=airquality[,1],
                    importance = "permutation") )
  rf.ImpScale(rgr.reg, scaling="mir")
  rf.ImpScale(rgr.reg, scaling="p", n=10)

}
\references{
Altmann, A., Tolosi, L., Sander, O. & Lengauer, T. (2010). Permutation importance: 
  a corrected feature importance measure, Bioinformatics 26:1340-1347.

Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas 
    connectivity in Yellowstone National Park with landscape genetics. 
    Ecology 91:252-261

Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species 
    distribution and change using Random Forests CH.8 in Predictive Modeling 
    in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer

Janitza, S., E. Celik, and A.-L. Boulesteix (2018). A computationally fast 
    variable  importance test for random forests for high-dimensional data. 
    Advances in Data Analysis and Classification 12(4):885-915
}
\seealso{
\code{\link[ranger]{importance_pvalues}} details on Altmann's p-value method
}
\author{
Jeffrey S. Evans  <jeffrey_evans@tnc.org>
}
