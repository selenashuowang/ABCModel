
getdata <- function (D, K, N, M,cor,error_s, error_i){

 
  sim.intercept=0
  

  Sim.utheta=diag(.5,2*K+D)
  Sim.utheta[upper.tri(Sim.utheta)]=c(0,cor,0)
  Sim.utheta=t(Sim.utheta)+Sim.utheta
  
  sim.U.Theta=mvrnorm(n=N,mu=rep(0,K*2+D),Sigma=Sim.utheta)
  sim.item=mvrnorm(n=M,mu=c(rep(1,D), 0),Sigma=diag(x=.5,nrow=D+1))
  sim.Alpha=sim.item[,1:D]
  sim.Beta=sim.item[,ncol(sim.item)]
  sim.U=sim.U.Theta[,1:K]
  sim.V=sim.U.Theta[,(K+1):(K+K)]
  sim.Theta=sim.U.Theta[,(2*K+1):(2*K+D)]
  sim.rho=0
  

  
  sim.P.i=sim.intercept +sim.U %*% t(sim.V)+ matrix(rnorm(n=N*N, mean=0, sd=sqrt(error_s)), N, N)

  ones=matrix(1,N,1)
  
  sim.P.ia=ones %*% sim.Beta + sim.Theta %*% t(sim.Alpha) + matrix(rnorm(n=N*M, mean=0, sd=sqrt(error_i)), N, M)
  

  
  X=ifelse(sim.P.i>0,1,0)
  Y=ifelse(sim.P.ia>0,1,0)

  
  return(list("X"=X,"Y"=Y,"Z"=sim.P.i, "H"=sim.P.ia, "U" =sim.U, "V"=sim.V,"Theta"=sim.Theta))
}
