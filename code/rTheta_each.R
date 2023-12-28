
rTheta_each<-function(i, Tl,offset_Y, U, invUTheta, invThetaU, ivTheta, s1, P,N){
  
  tmp.t=sapply(1:length(Tl), function(x) (Tl[[x]] - offset_Y[[x]])[i,], simplify = FALSE)
  
  a=lapply(tmp.t, function(x) ifelse(is.na(x),0,x))
  

  
  tmp.t1=Reduce('+', a)
  
  l<-(tmp.t1 /s1 - .5*invUTheta %*% matrix(U[i,])
      - .5*t(invThetaU) %*% matrix(U[i,]))
  
  iQ<- solve(  ivTheta + ( N*diag(1,nrow = P) )/s1  )
  
  return(iQ%*%l + t(chol(iQ))%*%rnorm(P) )
  ###########

  
}
