#' Simulate a binary network
#' 
#' Simulates a binary network
#' 
#' 
#' @usage simX_bin(EZ, rho)
#' @param EX square matrix giving the expected value of the latent Z matrix
#' @param rho dyadic correlation
#' @return a square binary matrix
#' @author Selena Wang
#' @export simX_bin
simX_bin <-
function(EZ,rho)
{
  ZS<-simZ(EZ,rho) 
  Xs<-1*(ZS>0) ; diag(Xs)<-NA
  Xs
}
