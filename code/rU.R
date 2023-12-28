#' Gibbs sampling of U
#' 
#' A Gibbs sampler for updating U.
#' 
#' @usage rU(Fl,U,Theta, Stheta, Sutheta, Su, s2=1, offset=offset)
#' @param Fl a list of V X V normal relational matrix
#' @param U V X K matrix containing current value of U
#' @param Theta D X V current value of Theta
#' @param Stheta D X D covariance of Theta
#' @param Sutheta D X K covariance between U and Theta
#' @param Su K X K matrix containing covariance of U
#' @param s2 dyadic variance
#' @param offset a list of the same dimension as Fl. It is assumed that 
#' Fl-offset follows a SRRM, so the offset should contain any multiplicative 
#' effects (such as \code{U\%*\% t(U) } )
#' @return \item{U}{a new value of U}
#' @author Selena Wang
#' 
#' 
#' @export rU
rU<-function(Fl,U,Theta, Stheta, Sutheta, Su, s2, offset=offset)
{
  K<-ncol(U) ; V<-nrow(U) ; N <-length(Fl)
  
  invUTheta <- (-1)*solve(Stheta - Sutheta %*% solve(Su) %*% t(Sutheta)) %*% Sutheta %*% solve(Su)
  
  invThetaU <- (-1)* solve( Su - t(Sutheta) %*% solve(Stheta) %*% Sutheta) %*% t(Sutheta) %*% solve(Stheta)
  
  
  ivU <- solve(Su - t(Sutheta) %*% solve(Stheta) %*% Sutheta)
  
  ## decorrelate
  to<-as.numeric(sqrt(solve(s2)))
  
  Es=sapply(1:length(Fl), function(x) (Fl[[x]]-offset[[x]])*to, simplify = FALSE)
  Es=lapply(Es, function(x) { if(is.na(diag(x)[1])){diag(x) <- 0}; x})
  
  if(K==1){U=matrix(sapply(1:V, function(x) rU_each(i=x, Es,U,Theta, invUTheta, invThetaU, ivU, to, K,N), simplify = TRUE))}
  if(K>1){U=t(sapply(1:V, function(x) rU_each(i=x, Es,U,Theta, invUTheta, invThetaU, ivU, to, K,N), simplify = TRUE))}
  
  
  
  
  
  
  # for(i in rep(sample(1:V),2))
  # {
  # 
  #   
  #   tmp.U=sapply(Es, function(x) apply(U*x[i,],2,sum,na.rm=TRUE), simplify = FALSE)
  #   
  #   tmp.U1=Reduce('+', tmp.U)
  #   
  #   l<- tmp.U1*(to) - .5*t(invUTheta) %*% Theta[i,]- .5*invThetaU %*% Theta[i,]
  #   iQ<- solve( ivU +    N*(to)^2*( crossprod(U) - U[i,]%*%t(U[i,]) ) ) 
  #   U[i,]<- iQ%*%l + t(chol(iQ))%*%rnorm(K) 
  # }
  
  
  ## update each U[i,]
  
  
  U
}







