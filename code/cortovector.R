cortovector<- function(S, K,P,N){

  Su=S[1:K,1:K]
  Stheta=S[(K+1):(K+P),(K+1):(K+P) ]
  Sutheta=S[(K+1):(K+P),1:K ]



  return(list("vc"=   c( Su[upper.tri(Su, diag = T)],
                         Stheta[upper.tri(Stheta, diag = T)]) , 
              "cov" =as.vector(Sutheta)))
}
