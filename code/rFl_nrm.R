#' Simulate missing values in a normal connectivity model
#' 
#' Simulates missing values of a sociomatrix under a normal connectivity model
#' 
#' 
#' @usage rFl_nrm(Z, EZ, s2, X)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param s2 dyadic variance
#' @param X square relational matrix
#' @return a square matrix, equal to  at non-missing values
#' @author Selena Wang
#' @export rFl_nrm
rFl_nrm<-function(Z,EZ,s2,X)
{
  ZS<-simX_nrm(EZ,s2)
  if(is.na(diag(EZ)[1])){
    diag(ZS)<-rnorm(nrow(X),rep(0,nrow(X)),sqrt(s2))
  }else{
    diag(ZS)<-rnorm(nrow(X),diag(EZ),sqrt(s2))
    
  }
  Z[is.na(X)]<-ZS[is.na(X)]  # this isn't quite right if there is asymmetric missingness. 
  Z
}

