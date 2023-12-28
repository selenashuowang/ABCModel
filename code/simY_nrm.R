#' Simulate a normal attribute matrix
#' 
#' Simulates a normal attribute matrix
#' 
#' 
#' @usage simY_nrm(EY, s1)
#' @param EY matrix giving the expected value of the attribute matrix
#' @param s1 variance
#' @return a V by P matrix
#' @author Selena Wang
#' @export simY_nrm
simY_nrm <-
function(EY,s1) 
{
  YS<-simH(EY,s1) 
  YS
}
