#' Calculate relative risk based on ratio of median unbiased estimators.
#' 
#' This function calculates an estimate of relative risk based on the ratio of two median unbiased estimates of proportions based on the work by \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2918902/}{Carter et al (2010)}.  The relative risk will be in the order of Pr(Group 1) / Pr(Group 2).
#' 
#' @param n1  Sample size for group 1. 
#' @param y1  Number of events in group 1.
#' @param n2  Sample size for group 2.
#' @param y2  Number of events in group 2.
#' @param alpha  The significance level for the confidence interval. Default value is 0.05.
#' @return  A dataframe with the various components generated during the estimation along with the MUE-based estimate of the relative risk.
#' @export
#' @examples
#' mue_rr(9,1,11,0)
#' mue_rr(9,1,11,0,0.05)
#' mue_rr(3,0,4,0,0.05)
#' mue_rr(12,1,15,1,0.15)




mue_rr <- function (n1, y1, n2, y2, alpha=0.05){
  p1_mue<-dobeta(n1,y1)
  p2_mue<-dobeta(n2,y2)
  mue_rr_estimate <- p1_mue/p2_mue
  ci<-mue_confidence_interval(n1, p1_mue ,n2, p2_mue, alpha)
  prop1 <- c(paste(y1,n1, sep="/"))
  prop2 <- c(paste(y2,n2,sep="/"))
  mue_rr <- data.frame(prop1,p1_mue,prop2,p2_mue, mue_rr_estimate,ci)
  return(mue_rr)
}
