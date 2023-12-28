
get_lambda_distribution<-function(K,D,N, N_sim=100000, sig = 0.05){
  
  beta_values<- matrix(NA, nrow = N_sim, ncol = D)
  
  for (d in 1:D){
    beta_values[,d]=rbeta(N_sim, shape1 = (N-K-D+d)/2, shape2 = K/2)
    
  }
  
  lambda_dis <- apply(beta_values, 1, prod)

  return(lambda_dis)
  
}

