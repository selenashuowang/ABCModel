leaveoneout_network <-function(all.super){
  N=nrow(all.super)
  dat=list()

  for(i in 1:N){
    tmp.training=all.super
    tmp.test=matrix(NA, nrow = N, ncol=N)
    tmp.training[i,]=NA
    tmp.test[i,]=all.super[i,]
  
    dat[[i]]=list("test"=tmp.test, "training"=tmp.training)
  }



  return(dat)
}
