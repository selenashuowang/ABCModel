leaveoneout_bipartite <-function(all.super){
  N=nrow(all.super)
  M=ncol(all.super)
  dat=list()
  
  for(i in 1:nrow(all.super)){
    tmp.training=all.super
    tmp.test=matrix(NA, nrow = N, ncol=M)
    tmp.training[i,]=NA
    tmp.test[i,]=all.super[i,]
    
    dat[[i]]=list("test"=tmp.test, "training"=tmp.training)
  }
  
  
  
  return(dat)
}
