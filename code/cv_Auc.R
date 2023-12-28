#' @title Assess the fit of the latent variable model when the data is binary using cross validation
#'
#' @description assess the out-of-sample fit of the model using ROC curves and auc values
#'
#' @param model a model returned by the jnirt, lvnm or lvirt function
#' @param fold number of folds for the K-fold cross validation
#' @param Iter number of iterations for the cross validation
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{aucs} a list of the area under the curve
#'  }
#'
#' @export
#'
#' @examples



cv_Auc<-function(model, fold = 10, K=NULL, D=NULL){
  
  
  if (!is.null(model$family_network) & !is.null(model$family_responses)){
    auc = matrix(NA, nrow = 2, ncol = fold)
    colnames(auc)=paste("iter",seq(1,fold),sep="")
    rownames(auc)=c("network","responses")
    
    X=model$X
    Y=model$Y
    
    for (n in 1:fold){
      Y.t = get_training_bipartite(Y,fold)[[n]]
      X.t = get_training_network(X,fold)[[n]]

      tmp<-jnirt_cross(X=X.t$training, Y=Y.t$training,family_network=model$family_network,
            family_responses=model$family_responses, K=K,D=D, model = model$input$model,
            rvar = model$input$rvar , cvar = model$input$cvar, dcor = model$input$dcor, 
            seed = 1, nscan = model$input$nscan, burn = model$input$burn, odens = model$input$odens,
            print = FALSE, gof=FALSE, 
            prior=model$input$prior)
      auc[,n] =Auc(tmp,X=X.t$test,Y=Y.t$test)
        }
      }


  
  if (!is.null(model$family_network) & is.null(model$family_responses)){
    X=model$X
    auc = matrix(NA, nrow = 1, ncol = fold)
    colnames(auc)=paste("iter",seq(1,fold),sep="")
    rownames(auc)=c("network")
    
    for (n in 1:fold){
      X.t = get_training_network(X,fold)[[n]]
      
      tmp<-lvnm_cross(X=X.t$training, family_network=model$family_network,indices =model$input$indices,
                K=K, model = model$input$model,independent = model$input$independent,
                 rvar = model$input$rvar , cvar = model$input$cvar, dcor = model$input$dcor, 
                 seed = 1, nscan = model$input$nscan, burn = model$input$burn, odens = model$input$odens,
                 print = FALSE, gof=FALSE, 
                 prior=model$input$prior)
      auc[,n] =Auc(tmp,X=X.t$test,Y=NULL)
    }
  }
    
  
  if (is.null(model$family_network) & !is.null(model$family_responses)){

    Y=model$Y
    auc = matrix(NA, nrow = 1, ncol = fold)
    colnames(auc)=paste("iter",seq(1,fold),sep="")
    rownames(auc)=c("responese")

    for (n in 1:fold){
      Y.t = get_training_bipartite(Y,fold)[[n]]

      tmp<-lvirt_cross(Y=Y.t$training,
                 family_responses=model$family_responses,D=D, 
                 seed = 1, nscan = model$input$nscan, burn = model$input$burn, odens = model$input$odens,
                 print = FALSE, gof=FALSE, 
                 prior=model$input$prior)
      auc[,n] =Auc(tmp,X=NULL,Y=Y.t$test)
    }
    
  }
  
  apply(auc,1,mean)

}

