covtocor<- function(cov_est, vc_est, K,P,N){
  Su=matrix(0, nrow = K,ncol=K )
  Stheta=matrix(0, nrow = P,ncol=P )
  Sutheta=matrix(cov_est, nrow = P,ncol=K )
  
  Su[upper.tri(Su, diag = T)] = vc_est[1:length(Su[upper.tri(Su, diag = T)])]
  Su=t(Su)+Su
  diag(Su)=diag(Su)/2
  Stheta[upper.tri(Stheta, diag = T)] = vc_est[(length(Su[upper.tri(Su, diag = T)])+1):(length(Su[upper.tri(Su, diag = T)])+sum(upper.tri(Stheta, diag = T)))]
  Stheta=t(Stheta)+Stheta
  diag(Stheta)=diag(Stheta)/2
  
  S=matrix(NA, nrow = K+P,ncol=K+P )
  colnames(S)=c(paste("U",seq(1,K), sep = ""), paste("Theta", seq(1,P), sep = ""))
  rownames(S)=c(paste("U",seq(1,K), sep = ""), paste("Theta", seq(1,P), sep = ""))
  
  S[1:K,1:K]=Su
  S[(K+1):(K+P),(K+1):(K+P) ]=Stheta
  S[(K+1):(K+P),1:K ]=Sutheta
  S[(1):(K),(K+1):(K+P) ]=t(Sutheta)
  
  cor_est=cov2cor(S)
  
  return(cor_est)
}
