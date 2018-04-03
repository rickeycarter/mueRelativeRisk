#' Generate the exact bootstrap interval for the MUE estimate.
#' 
#' This function enumerates the exact bootstrap confidence interval for the MUE-based estimate of the relative risk.
#' 
#' @param n1  Sample size for group 1. 
#' @param p1  MUE estimate calculated by the dobeta function for group 1.
#' @param n2  Sample size for group 2.
#' @param p2  MUE estimate calculated by the dobeta function for group 2.
#' @param alpha  The significance level for the confidence interval.
#' @return  Returns a matrix of the upper and lower limits for the confidence interval.
#' @export

mue_confidence_interval <- function(n1,p1, n2, p2, alpha) {
  # Parameter information:
  #    n1 and n2: Sample sizes for the two groups
  #    p1 and p2: these are the MUE based estimates calculated by the dobeta function
  
  
  # expand the grid to prepare for the exact bootstrap confidence interval
  # this grid will include all possible success numbers for fixed n1 and n2
  # the binomial proportions are assumed constant across the samples
  
  ci_long <- as.data.frame(expand.grid(y1 = seq(0,n1), n1=n1, p1=p1, y2=seq(0,n2) ,n2 = n2, p2=p2 ))
  ci_long <- as.data.frame(apply(ci_long,2,as.numeric))
  
  # call the dobeta function to return the MUE estimate for p1 and p2
  mue1<-dobeta(ci_long$n1,ci_long$y1)
  mue2<-dobeta(ci_long$n2,ci_long$y2)
  # define the MUE-based estimate of the RR as the ratio of MUE estimates for p1 and p2
  rr_mue <- mue1/mue2
  
  #compute the product binomial probability for the individual sample configuration
  prodbin <- stats::dbinom(x=ci_long$y1,size=ci_long$n1,prob=ci_long$p1) * stats::dbinom(x=ci_long$y2,size=ci_long$n2,prob=ci_long$p2)
  
  # piece the vectors back together and return the estimates 
  ci_long <-cbind(ci_long,prodbin,mue1,mue2,rr_mue)
  
  
  
  
  ## Need rr_mue may not be unique, so we need to aggregate over values to sum the probabilities
  ## After aggregating to produce the emprical probability mass function, order and sort to yield
  ## the cumulative distribution function. This function will be used to determine the alpha/2 tails
  
  ci_expanded_sum <-stats::aggregate(prodbin ~ rr_mue,data=ci_long,FUN=sum)
  ci_expanded_sum_sorted <-ci_expanded_sum[order(ci_expanded_sum$rr_mue),]
  ci_expanded_sum_sorted$prob_sum <- cumsum(ci_expanded_sum_sorted$prodbin)
  ci_expanded_sum_sorted$prob_sum_comp <- 1 - ci_expanded_sum_sorted$prob_sum + ci_expanded_sum_sorted$prodbin
  
  alpha_2 <- rep(alpha/2,nrow(ci_expanded_sum_sorted))
  
  ci_expanded_sum_sorted<-cbind(ci_expanded_sum_sorted, alpha_2)
  
  ci_expanded_sum_sorted$ilower <-ifelse(ci_expanded_sum_sorted$prob_sum < ci_expanded_sum_sorted$alpha_2, 1, 0)
  ci_expanded_sum_sorted$iupper <-ifelse(ci_expanded_sum_sorted$prob_sum_comp > ci_expanded_sum_sorted$alpha_2, 1, 0)
  
  
  #### user linear interpolation for the endpoints of the CI
  ### Establish lower Limit
  #pulls of the row where indicator changes from 0 to 1
  index_low<-which(!duplicated(ci_expanded_sum_sorted$ilower))[2]
  
  #If index_low = 1, then the first possible RR estimate has probability > alpha/2. In this case, set the value
  # for the lower limit to be 0
  
  if (index_low[1] == 1) {
    mue_ci_lower <- 0
  } else {
    endpoints_mue<-ci_expanded_sum_sorted$rr_mue[(index_low-1):index_low]
    endpoints_prob <- ci_expanded_sum_sorted$prob_sum[(index_low-1):index_low]
    mue_ci_lower <- (
      endpoints_mue[1]*(endpoints_prob[2]- ci_expanded_sum_sorted$alpha_2[1]) + 
        endpoints_mue[2]*(ci_expanded_sum_sorted$alpha_2[1] - endpoints_prob[1])
    ) / 
      (endpoints_prob[2] - endpoints_prob[1])
  }
  
  
  
  ### Establish Upper limit
  index_high<-which(!duplicated(ci_expanded_sum_sorted$iupper,fromLast=T))[1]
  
  #If index_high = 1, then the largest possible RR estimate has probability > alpha/2. In this case, set the value
  # for the upper limit to be large - 10**100
  
  if (index_high[1] == NROW(ci_expanded_sum_sorted$rr_mue)) {
    mue_ci_lower <- 10**100
  } else {
    endpoints_mue<-ci_expanded_sum_sorted$rr_mue[(index_high+1):index_high]
    endpoints_prob <- ci_expanded_sum_sorted$prob_sum_comp[(index_high+1):index_high]
    
    mue_ci_upper <- (
      endpoints_mue[1]*(endpoints_prob[2]- ci_expanded_sum_sorted$alpha_2[1]) +
        endpoints_mue[2]*(ci_expanded_sum_sorted$alpha_2[1] - endpoints_prob[1])
    ) / 
      (endpoints_prob[2] - endpoints_prob[1])
  }
  
  mue_ci <-  cbind(mue_ci_lower, mue_ci_upper)
  return(mue_ci)
}
