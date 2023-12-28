#' @title Select the number of dimensions using cross validation
#'
#' @description select the numebr of dimensions in the data based on out-of-sample auc values
#'
#' @param model a model returned by the jnirt, lvnm or lvirt function
#' @param fold number of folds for the K-fold cross validation
#' @param D maximum number of dimensions for the social network if applicable
#' @param K maximum number of dimensions for the item responses if applicable
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{aucs} an array of the area under the curve
#'  }
#'
#' @export
#'
#' @examples

rselect_dimension<-function(model,D=NULL,K=NULL,fold=5){
  
  if (!is.null(model$family_network) & !is.null(model$family_responses)){
    auc_a=array(NA, dim = c(K,D,2), dimnames = list(paste("K = ",seq(1,K),sep=""),
                                                    paste("D = ",seq(1,D),sep=""),
                                                    c("network","responses")))
    
    for(d in 1:D){
      for(k in 1:K){
        auc_a[k,d,]=cv_Auc(model, fold = fold, K=k, D=d)
      }
    }
    
  }
  
  if (!is.null(model$family_network) & is.null(model$family_responses)){
    auc_a=array(NA, dim = c(K,1), dimnames = list(paste("K = ",seq(1,K),sep=""),
                                                    c("network")))
    

      for(k in 1:K){
        auc_a[k,]=cv_Auc(model, fold = fold, K=k)
      }
    }
    
  
  
  if (is.null(model$family_network) & !is.null(model$family_responses)){
    auc_a=array(NA, dim = c(D,1), dimnames = list(paste("D = ",seq(1,D),sep=""),
                                                  c("responses")))
    
    
    for(d in 1:D){
      auc_a[d,]=cv_Auc(model, fold = fold, D=d)
    }
  }

  auc_a
}




