#' @title Assess the fit of the APLSM
#'
#' @description assess the fit of the model using ROC curves and auc values
#'
#' @param Y.i N by N matrix containing the binary social network
#' @param Y.ia N by M matrix containing the binary multivariate attributes
#' @param type character indicating the types of model. It could be "DD", distance by distance model, "DV", distance by vector model,
#'  "VV", vector by vector model
#' @param model object of class the APLSM
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{Yi.auc} scaler of the area under the curve for the social network
#'  \item \code{Ya.auc} scaler of the area under the curve for the multivariate covariates
#'  }
#'
#' @export
#'
#' @examples
#' attach(french)
#' b=aplsm(Niter=3,Y.i, Y.ia,D=2, type="DD")
#' GOFaplsm(b, "DD",Y.i, Y.ia)


Gofaplsm<-function(model, type,Y.i, Y.ia){

  est.alpha.0 = model$lsmhAlpha.0
  est.alpha.1 = model$lsmhAlpha.1
  Z.i = model$lsmhEZ.i
  Z.a = model$lsmhEZ.a
  D=nrow(model$lsmhVZ.0)
  M=nrow(model$lsmhEZ.a)
  N=nrow(model$lsmhEZ.i)


  Ps=Predictaplsm(model,type)


  #a=roc(c(Y.i),c(Ps[[1]]),auc.polygon=FALSE, grid=FALSE, plot=FALSE,auc=TRUE)

  a=roc(c(Y.i),c(Ps[[1]]),auc.polygon=FALSE, grid=FALSE,
        ylab="true positive rate",xlab="false positive rate",xlim=c(1,0),plot=TRUE,auc=TRUE,main=NULL, font.axis = 2,bty="n", cex.axis =1.2,ann=FALSE,  col="black",legacy.axes = TRUE, lwd=5)
  text(x=0.45, y=.1,cex=1.6, labels=paste("AUC = ", as.character(round(as.numeric(a$auc),digits = 4)),sep = ""), font= 2, col="black")



  #c=roc(c(Y.ia),c(Ps[[2]]),auc.polygon=FALSE, grid=FALSE,plot=FALSE,auc=TRUE)

  c=roc(c(Y.ia),c(Ps[[2]]),auc.polygon=FALSE, grid=FALSE,
      ylab="true positive rate",xlab="false positive rate",xlim=c(1,0),plot=TRUE,auc=TRUE,main=NULL, font.axis = 2,bty="n", cex.axis =1.2,ann=FALSE,  col="black",legacy.axes = TRUE, lwd=5)
  text(x=0.45, y=.1,cex=1.6, labels=paste("AUC = ", as.character(round(as.numeric(c$auc),digits = 4)),sep = ""), font= 2, col="black")


  return(c("Yi.auc" = a$auc,"Yia.auc" = c$auc))

}

