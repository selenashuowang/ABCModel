#' Simulate missing values in a normal attribute model
#' 
#' Simulates missing values under an attribute model
#' 
#' 
#' @usage rTl_nrm(H, EH,s1, Y)
#' @param H a V x P matrix, the current value of H
#' @param EH expected value of H
#' @param s1 attribute variance
#' @param Y attribute matrix
#' @return an attribute matrix, equal to  at non-missing values
#' @author Selena Wang
#' @export rTl_nrm
rTl_nrm<-function(H,EH,s1,Y)
{
  HS<-simY_nrm(EH,s1)
  H[is.na(Y)]<-HS[is.na(Y)]  # this isn't quite right if there is asymmetric missingness. 
  H
}

