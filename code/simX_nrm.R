#' Simulate a normal connectivity matrix
#' 
#' Simulates a normal connectivity matrix
#' 
#' 
#' @usage simX_nrm(EX, s2)
#' @param EX square matrix giving the expected value of the connectivity matrix
#' @param s2 dyadic variance
#' @return a square matrix
#' @author Selena Wang
#' @export simX_nrm
simX_nrm <-
function(EX,s2) 
{
  XS<-simZ(EX,s2) 
  diag(XS)<-NA
  XS
}
