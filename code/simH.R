#' Simulate H given its expectation and covariance
#' 
#' Simulate H given its expectation and covariance
#' 
#' 
#' @usage simH(EH, s1 = 1)
#' @param EH expected value of H
#' @param s1 attribute variance
#' @return a simulated value of H
#' @author Selena Wang
#' @export simH
simH <-
function(EH,s1=1)
{ 
  EH+rnorm(nrow(EH)*ncol(EH), sd=sqrt(s1)) 
} 






