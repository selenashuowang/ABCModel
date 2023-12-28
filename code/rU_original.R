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
rU_original<-function(Fl,U,Theta, Stheta, Sutheta, Su, s2, offset=offset)
{
  K<-ncol(U) ; V<-nrow(U) ; N <-length(Fl)
  
  invUTheta <- (-1)*solve(Stheta - Sutheta %*% solve(Su) %*% t(Sutheta)) %*% Sutheta %*% solve(Su)
  
  invThetaU <- (-1)* solve( Su - t(Sutheta) %*% solve(Stheta) %*% Sutheta) %*% t(Sutheta) %*% solve(Stheta)
  
  
  ivU <- solve(Su - t(Sutheta) %*% solve(Stheta) %*% Sutheta)
  
  ## decorrelate
  to<-as.numeric(sqrt(solve(s2)))
  
  
  for(i in rep(sample(1:V),4))
  {
    tmp=0
    for(n in 1: length(Fl)){
      E<-Fl[[n]]-offset[[n]]
      Es<-to*E
      #rnorm(length(diag(Es)),mean=0,sd=sqrt(s2))
      
      if(is.na(diag(Es)[1])){diag(Es)=0}
      
      tmp= tmp + apply(U*Es[i,],2,sum,na.rm=TRUE) -  U[i,]*Es[i,i] 
    }
    l<- tmp*(to) - .5*t(invUTheta) %*% Theta[i,]- .5*invThetaU %*% Theta[i,]
    iQ<- solve( ivU +    N*(to)^2*( crossprod(U) - U[i,]%*%t(U[i,]) ) ) 
    U[i,]<- iQ%*%l + t(chol(iQ))%*%rnorm(K) 
  }
  
  
  ## update each U[i,]
  
  
  U
}







