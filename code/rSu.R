#' Gibbs update for multiplicative effects covariance
#' 
#' @usage rSu(U,Su, Su0=NULL,etau=NULL) 
#' @param U n X K current value of U
#' @param Su0 prior (inverse) scale matrix for the prior distribution
#' @param etau prior degrees of freedom for the prior distribution
#' @author Selena Wang
#' @export rSu
#'
rSu<-function(U, Su0=NULL,etau=NULL) 
{
  if(is.null(Su0)){ Su0<-diag(ncol(U))  } 
  if(is.null(etau)){ etau<-ncol(U)+2 }
  Q=U
  S=solve(rwish(solve(etau*Su0+t(Q) %*% Q), etau+nrow(U)))
  
    list("Su" = S)


}


