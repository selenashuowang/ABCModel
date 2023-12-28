#' @title Assess the fit of the latent variable model when the data is binary
#'
#' @description assess the fit of the model using ROC curves and auc values
#'
#' @param model a model returned by the jnirt, lvnm or lvirt function
#' @param X a social network, if applicable
#' @param Y an item response matrix, if applicable
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{aucs} a list of the area under the curve
#'  }
#'
#' @export
#'
#' @examples



Auc<-function(model,X=NULL,Y=NULL){
  
  if (!is.null(model$family_network) & !is.null(model$family_responses)){

    a=roc(c(X),c(model$EZ),auc.polygon=FALSE, grid=FALSE,
          ylab="true positive rate",xlab="false positive rate",xlim=c(1,0),plot=TRUE,auc=TRUE,main=NULL, font.axis = 2,bty="n", cex.axis =1.2,ann=FALSE,  col="black",legacy.axes = TRUE, lwd=5)
    text(x=0.45, y=.1,cex=1.6, labels=paste("AUC = ", as.character(round(as.numeric(a$auc),digits = 4)),sep = ""), font= 2, col="black")
    
    
    c=roc(c(Y),c(model$EH),auc.polygon=FALSE, grid=FALSE,
          ylab="true positive rate",xlab="false positive rate",xlim=c(1,0),plot=TRUE,auc=TRUE,main=NULL, font.axis = 2,bty="n", cex.axis =1.2,ann=FALSE,  col="black",legacy.axes = TRUE, lwd=5)
    text(x=0.45, y=.1,cex=1.6, labels=paste("AUC = ", as.character(round(as.numeric(c$auc),digits = 4)),sep = ""), font= 2, col="black")

    aucs<-c("network.auc" = a$auc,"responese.auc" = c$auc)
  }

  

  
  if (!is.null(model$family_network) & is.null(model$family_responses)){


    
    a=roc(c(X),c(model$EZ),auc.polygon=FALSE, grid=FALSE,
          ylab="true positive rate",xlab="false positive rate",xlim=c(1,0),plot=TRUE,auc=TRUE,main=NULL, font.axis = 2,bty="n", cex.axis =1.2,ann=FALSE,  col="black",legacy.axes = TRUE, lwd=5)
    text(x=0.45, y=.1,cex=1.6, labels=paste("AUC = ", as.character(round(as.numeric(a$auc),digits = 4)),sep = ""), font= 2, col="black")
    
  
    aucs<-c("network.auc" = a$auc)
    
  }
  
  if (is.null(model$family_network) & !is.null(model$family_responses)){



    c=roc(c(Y),c(model$EH),auc.polygon=FALSE, grid=FALSE,
          ylab="true positive rate",xlab="false positive rate",xlim=c(1,0),plot=TRUE,auc=TRUE,main=NULL, font.axis = 2,bty="n", cex.axis =1.2,ann=FALSE,  col="black",legacy.axes = TRUE, lwd=5)
    text(x=0.45, y=.1,cex=1.6, labels=paste("AUC = ", as.character(round(as.numeric(c$auc),digits = 4)),sep = ""), font= 2, col="black")
    
    aucs<-c("responese.auc" = c$auc)
    
  }
  aucs

}

