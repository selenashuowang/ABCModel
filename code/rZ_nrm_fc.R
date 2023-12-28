#' Simulate missing values in a normal AME model
#' 
#' Simulates missing values of a sociomatrix under a normal AME model
#' 
#' 
#' @usage rZ_nrm(Z, EZ, rho,s2, X)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @param X square relational matrix
#' @return a square matrix, equal to X at non-missing values
#' @author Selena Wang
#' @export rZ_nrm
rZ_nrm<-function(Z,EZ,rho,s2,X)
{
  ZS<-simX_nrm(EZ,rho,s2)
  diag(ZS)<-rnorm(nrow(X),diag(EZ),sqrt(s2*(1+rho)))
  Z[is.na(X)]<-ZS[is.na(X)]  # this isn't quite right if there is asymmetric missingness. 
  Z
}

