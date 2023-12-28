
rU_con_each<-function(i, Es,U,ivU, to, K,N){
  
  
  
  tmp.U=sapply(Es, function(x) apply(U*x[i,],2,sum,na.rm=TRUE)-  U[i,]*x[i,i] , simplify = FALSE)
  
  tmp.U1=Reduce('+', tmp.U)
  
  l<- tmp.U1*(to) 
  iQ<- solve( ivU +    N*(to)^2*( crossprod(U) - U[i,]%*%t(U[i,]) ) ) 
  return(iQ%*%l + t(chol(iQ))%*%rnorm(K) )
}


