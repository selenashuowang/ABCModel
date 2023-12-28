#' Gibbs sampling of Theta
#' 
#' A Gibbs sampler for updating the Person latent effect Theta.
#' 
#' @usage rTheta(Tl, beta, Alpha, Theta ,U, Stheta, Sutheta, Su, s1)
#' @param Tl V X P normal item response matrix
#' @param Theta V X P current value of Theta
#' @param U V X K matrix containing current value of U
#' @param s1 item response variance
#' @param Stheta P X P covariance of Theta
#' @param Sutheta P X K covariance between U and Theta
#' @param Su K X K matrix containing covariance of U
#' @param offset_Y a list of the same dimension as Tl. 
#' @return \item{Theta}{a new value of Theta}
#' @author Selena Wang
#' @export rTheta
rTheta<-function(Tl, Theta ,U, Stheta, Sutheta, Su, s1, offset_Y)
{
  P<-ncol(Theta) ; V<-nrow(Theta); N <- length(Tl)
  

  
  ivTheta <- solve(Stheta - Sutheta %*% solve(Su) %*% t(Sutheta))
  
  invUTheta <- (-1)*solve(Stheta - Sutheta %*% solve(Su) %*% t(Sutheta)) %*% Sutheta %*% solve(Su)
  
  invThetaU <- (-1)* solve( Su - t(Sutheta) %*% solve(Stheta) %*% Sutheta) %*% t(Sutheta) %*% solve(Stheta)
  
  if(P==1){Theta=matrix(sapply(1:V, function(x) rTheta_each(i=x, Tl,offset_Y, U, invUTheta, invThetaU, ivTheta, s1, P,N), simplify = TRUE))}
  if(P>1){Theta=t(sapply(1:V, function(x) rTheta_each(i=x, Tl,offset_Y, U, invUTheta, invThetaU, ivTheta, s1, P,N), simplify = TRUE))}
  
  # 
  # for(i in rep(sample(1:V),2))
  # {
  #   tmp=0
  #   
  #   
  #   for(n in 1: length(Tl)){
  #     Tl[[n]][is.na(Tl[[n]])]=0
  #     t<-Tl[[n]] - offset_Y[[n]]
  # 
  #     tmp= tmp + t[i,]
  #   }
  #   l<-(tmp /s1 - .5*invUTheta %*% matrix(U[i,])
  #       - .5*t(invThetaU) %*% matrix(U[i,]))
  #   
  #   iQ<- solve(  ivTheta + ( N*diag(1,nrow = P) )/s1  )
  #   
  #   Theta[i,]<- iQ%*%l + t(chol(iQ))%*%rnorm(P) 
  #   
  # }

  Theta 
}







