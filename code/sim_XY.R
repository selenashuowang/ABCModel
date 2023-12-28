#' simulate from the JNIRT model
#' 
#' 
#' 
#' @usage sim_XY(model,  family_network,family_responses)
#' @param model a latent variable model

#' @return a list of simulated data
#' @author Selena Wang
#' 
#' 
#' 
#' @export simXY
sim_XY <- function(model){
  EZ<-model$EZ
  EH<-model$EH
  family_network<-model$family_network
  family_responses<-model$family_responses
  # 

    if(family_network=="bin") { Xs<-simX_bin(EZ,rho=mean(model$VC[,"rho"])) }
    if(family_network=="nrm") { Xs<-simX_nrm(EZ,rho=mean(model$VC[,"rho"]),s2=mean(model$VC[,"ve_network"])) }
    
    if(family_responses=="bin") { Ys<-simY_bin(EH) }
    if(family_responses=="nrm") { Ys<-simY_nrm(EH,s1=mean(model$VC[,"ve_responses"])) }

  list(X=Xs, Y=Ys)
  

}

