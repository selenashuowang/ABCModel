


cv_Auc_separate<-function(K=NULL, D=NULL,X,Y,nscan){
  

  

  if((!is.null(X)) & (is.null(Y))){
    
    auc = matrix(NA, nrow = 2, ncol = nrow(X))
    colnames(auc)=paste("iter",seq(1,nrow(X)),sep="")
    rownames(auc)=c("network","responses")
  
    for (n in 1:nrow(X)){

      X.t = leaveoneout_network(X)[[n]]
      
      if((sd(c(X.t$test), na.rm = TRUE) !=0)){
        tmp<- ame(Y=X.t$training, family = "bin", R=K, rvar = FALSE ,
                  cvar = FALSE,  dcor = TRUE,
                  intercept=TRUE,
                  symmetric=FALSE, seed = 1, nscan = nscan, burn = 1, odens = 1, plot=FALSE, print = FALSE, gof=TRUE,
                  prior=list())

   
        
        a=roc(c(X.t$test),c(tmp$EZ),auc.polygon=FALSE, grid=FALSE,
              ylab="true positive rate",xlab="false positive rate",xlim=c(1,0),plot=TRUE,auc=TRUE,main=NULL, font.axis = 2,bty="n", cex.axis =1.2,ann=FALSE,  col="black",legacy.axes = TRUE, lwd=5)
        
        
        auc[1,n] = a$auc
      }
    }
      
      
  }
  
  
  if((is.null(X)) & (!is.null(Y))){
    
    auc = matrix(NA, nrow = 2, ncol = nrow(Y))
    colnames(auc)=paste("iter",seq(1,nrow(Y)),sep="")
    rownames(auc)=c("network","responses")
    
    for (n in 1:nrow(Y)){
      Y.t = leaveoneout_bipartite(Y)[[n]]

      
      if((sd(c(Y.t$test), na.rm = TRUE) !=0)){
        tmp<- lvirt(Y.t$training,family_responses="nrm", D=D, 
                    seed = 1, nscan = nscan, burn = 1, odens = 1,
                    print = FALSE, gof=TRUE, indices_irt = NULL,
                    prior=list())
        
        
    
        auc[2,n] = mean(abs(Y.t$test-tmp$EH)/Y, na.rm = TRUE)
      }
    }
    
    
  }
  
  
  
  
  return(auc)
  
}


