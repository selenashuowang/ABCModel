#' Simulate Z given its expectation and covariance
#' 
#' Simulate Z given its expectation and covariance
#' 
#' 
#' @usage simZ(EZ, rho, s2 = 1)
#' @param EZ expected value of Z
#' @param s2 dyadic variance
#' @return a simulated value of Z
#' @author Selena Wang
#' @export simZ
simZ <-
function(EZ,s2=1)
{ 
  V<-nrow(EZ)
  tmp<-rnorm(V*(V-1)/2, sd=sqrt(s2))
  EC<-matrix(0,nrow(EZ),nrow(EZ))
  EC[upper.tri(EC,diag = FALSE)]<-tmp
  EZ+  t(EC) + EC    
} 

