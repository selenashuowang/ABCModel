#' @title Assess the fit of the LSM
#'
#' @description assess the fit of the model using ROC curves and auc values
#'
#' @param Y.i N by N matrix containing the binary item response matrix
#' @param model object of class LSM from lvm4net
#'
#' @return scalar containing:
#'  \itemize{
#'  \item \code{Yi.auc} scaler of the area under the curve for the social network
#'  }
#'
#' @export


Goflsm<-function(model,Y.i){
  
  est.alpha.1 = model$xiT
  Z.i = model$lsmEZ
  
  
  
  D=nrow(model$lsmEZ)
  N=nrow(Z.i)
  
  
  
  
  
  est.P.i=Predictlsm(model)
  
  
  a=roc(c(Y.i),c(est.P.i),auc.polygon=FALSE, grid=FALSE,
        ylab="true positive rate",xlab="false positive rate",xlim=c(1,0),plot=TRUE,auc=TRUE,main=NULL, font.axis = 2,bty="n", cex.axis =1.2,ann=FALSE,  col="black",legacy.axes = TRUE, lwd=5)
  text(x=0.45, y=.1,cex=1.6, labels=paste("AUC = ", as.character(round(as.numeric(a$auc),digits = 4)),sep = ""), font= 2, col="black")
  
  

  return("Yi.auc" = a$auc)
  
}