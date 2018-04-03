#' Function to determine the MUE estimate based on the beta distribution.
#' 
#' This function calculates the median unbiased estimate for a single proportion.
#' 
#' @param n  Number of independent trials. 
#' @param y  Number of successes.
#' @return  A number with the MUE estimate
#' @export


dobeta <- function(n,y){
  #   browser()
  lower <- ifelse(y==0, 0, ifelse(y==n,.5 ** (1/n) , stats::qbeta(.5,y,n-y+1) ))
  upper <- ifelse(y==0, 1-.5**(1/n), ifelse(y==n, 1 , stats::qbeta(.5,y+1,n-y) ))
  mue <- (lower + upper) / 2
  return(mue)
}