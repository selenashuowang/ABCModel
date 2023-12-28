
get_training_bipartite <-function(all.super,K){

  N=nrow(all.super)

  M = ncol(all.super)

  ID=seq(1,N*M)

  folds=sample(cut(ID,breaks=K,labels=FALSE))
  dat=list()

  for(i in 1:K){
    testID=which(folds==i)
    trainingID <- ID [!ID %in% testID]




    rowindex=seq(1,N)
    colindex=seq(1,M)

    df_index=expand.grid(rowindex,colindex)

    all.train=all.super
    all.test=all.super
    for(each in testID){
      all.train[df_index[each,1],df_index[each,2]]=NA
    }

    for(each in trainingID){
      all.test[df_index[each,1],df_index[each,2]]=NA
    }

    dat[[i]]=list("test"=all.test, "training"=all.train)

  }



  return(dat)
}



