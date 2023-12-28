#' Simulate a binary item response matrix
#' 
#' Simulates a binary item response matrix
#' 
#' 
#' @usage simY_bin(EH)
#' @param EH N x M matrix giving the expected value of the latent H matrix
#' @return an item response matrix
#' @author Selena Wang
#' @export simY_bin
simY_bin <-
function(EH)
{
  ZS<-simH(EH,s1 =1 ) 
  YS<-1*(ZS>0) ; diag(YS)<-NA
YS
}
