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
#' @export get_gof
get_gof <- function(model, N){
  
  if (!is.null(model$family_network) & !is.null(model$family_responses)){
    X=model$X
    Y=model$Y

    gofXY<-c(gofstats_network(X), gofstats_responses(Y))
    
    GOF<-matrix(gofXY,1,length(gofXY))  
    rownames(GOF)<-"obs"
    colnames(GOF)<-names(gofXY)
    
    for (i in 1:N){
      df=sim_XY(model)
      Xs=df$X
      Ys=df$Y
      Xs[is.na(X)]<-NA ; Ys[is.na(Y)]<-NA; 
      GOF<-rbind(GOF,c(gofstats_network(Xs),gofstats_responses(Ys))) 
     
    }
  }
    
  if (!is.null(model$family_network) & is.null(model$family_responses)){
    X=model$X
    family_network=model$family_network
    gofXY<-c(gofstats_network(X))
    GOF<-matrix(gofXY,1,length(gofXY))  
    rownames(GOF)<-"obs"
    colnames(GOF)<-names(gofXY)
    
    for (i in 1:N){
      if(family_network=="bin") { Xs<-simX_bin(model$EZ,mean(model$VC[,"rho"])) }
      if(family_network=="nrm") { Xs<-simX_nrm(model$EZ,
                                               rho=mean(model$VC[,"rho"]),
                                               s2=mean(model$VC[,"ve_network"])) }
      
      Xs[is.na(X)]<-NA ; GOF<-rbind(GOF,c(gofstats_network(Xs))) 
      
    }
    
  }
  
  if (is.null(model$family_network) & !is.null(model$family_responses)){
    Y=model$Y
    gofXY<-c( gofstats_responses(Y))
    family_responses=model$family_responses

    GOF<-matrix(gofXY,1,length(gofXY))  
    rownames(GOF)<-"obs"
    colnames(GOF)<-names(gofXY)
    
    for (i in 1:N){
      
      if(family_responses=="bin") { Ys<-simY_bin(model$EH) }
      if(family_responses=="nrm") { Ys<-simY_nrm(model$EH,
                                                 s1=mean(model$VC[,"ve_responses"])) }

      Ys[is.na(Y)]<-NA ; GOF<-rbind(GOF,c(gofstats_responses(Ys))) 
      
    }
    
  }
  
  GOF
  

  
}

