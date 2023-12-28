#' get goodness of fit for a latent variable model
#' 
#' Goodness of fit statistics 
#' 
#' 
#' @usage get_gof(model, N)
#' @param model a latent variable model
#' @param N number of simulated data
#' @return goodness of fit statistics 
#' @author Selena Wang
#' @examples
#' 
#' 
#' 
#' @export get_gof_ame
get_gof_ame <- function(X, N, family, model){
  
    gofXY<-c(gofstats_network(X))
    GOF<-matrix(gofXY,1,length(gofXY))  
    rownames(GOF)<-"obs"
    colnames(GOF)<-names(gofXY)
    
    for (i in 1:N){
      if(family=="bin") { Xs<-simX_bin(model$EZ,mean(model$VC[,"rho"])) }
      if(family=="nrm") { Xs<-simX_nrm(model$EZ,
                                               rho=mean(model$VC[,"rho"]),
                                               s2=mean(model$VC[,"ve"])) }
      
      Xs[is.na(X)]<-NA ; GOF<-rbind(GOF,c(gofstats_network(Xs))) 
      
    }

  GOF
  
  
  
}

