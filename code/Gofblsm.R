#' @title Assess the fit of the BLSM
#'
#' @description assess the fit of the model using ROC curves and auc values
#'
#' @param Y.ia N by M matrix containing the binary item response matrix
#' @param model object of class BLSM
#'
#' @return scalar containing:
#'  \itemize{
#'  \item \code{Ya.auc} scaler of the area under the curve for the multivariate covariates
#'  }
#'
#' @export
#'
#' @examples
#' attach(french)
#' a=blsm(Niter=5,Y.ia,D=2)
#' Gofblsm(a,Y.ia)

Gofblsm<-function(model,Y.ia){

  est.alpha.1 = model$lsmbAlpha.1
  Z.i = model$lsmbEZ.i
  Z.a = model$lsmbEZ.a

  M=ncol(Y.ia)
  N=nrow(Y.ia)




  est.P.ia=Predictblsm(model)


  c=roc(c(Y.ia),c(est.P.ia),auc.polygon=FALSE, grid=FALSE,
        ylab="true positive rate",xlab="false positive rate",xlim=c(1,0),plot=TRUE,auc=TRUE,main=NULL, font.axis = 2,bty="n", cex.axis =1.2,ann=FALSE,  col="black",legacy.axes = TRUE, lwd=5)
  text(x=0.45, y=.1,cex=1.6, labels=paste("AUC = ", as.character(round(as.numeric(c$auc),digits = 4)),sep = ""), font= 2, col="black")

  return("Yia.auc" = c$auc)

}
