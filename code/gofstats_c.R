#' Goodness of fit statistics
#' 
#' Goodness of fit statistics evaluating second and third-order dependence
#' patterns
#' 
#' 
#' @usage gofstats_c(Y)
#' @param Y a relational data matrix
#' @return a vector of gof statistics
#' @author Selena Wang
#' @examples
#' 
#' 
#' gofstats_c(YX_nrm$Y) 
#' 
#' 
#' @export gofstats_c
gofstats_c<-function(Y)
{

 
    E<-Y-mean(Y,na.rm=TRUE)
    D<-1*(!is.na(E)) ; E[is.na(E)]<-0
    triad.dep<-c( 
      sum(diag(E%*%E%*%E))/( sum(diag(D%*%D%*%D))*sd(c(Y),na.rm=TRUE)^3), 
      sum(diag(E%*%t(E)%*%E))/( sum(diag(D%*%t(D)%*%D))*sd(c(Y),na.rm=TRUE)^3) ) 
    
    gof<-c(mean(Y, na.rm=TRUE),sd(Y, na.rm=TRUE), triad.dep ) 
    
    gof[is.na(gof)]<-0 
    
    names(gof)<-c("mean","sd","cycle.dep","trans.dep")
    
    gof

  }

 




