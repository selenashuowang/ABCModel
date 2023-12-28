#' Goodness of fit statistics
#' 
#' Goodness of fit statistics evaluating second and third-order dependence
#' patterns
#' 
#' 
#' @usage gofstats_c(Y)
#' @param Y an attribute data matrix
#' @return a vector of gof statistics
#' @author Selena Wang
#' @examples
#' 
#' 
#' gofstats_c(YX_nrm$Y) 
#' 
#' 
#' @export gofstats_c
gofstats_a<-function(Y)
{

 
    
    gof<-c(mean(Y, na.rm=TRUE),sd(Y, na.rm=TRUE) ) 
    
    gof[is.na(gof)]<-0 
    
    names(gof)<-c("mean","sd")
    
    gof

  }

 




