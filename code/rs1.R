#' Gibbs update for attribute variance
#' 
#' Gibbs update for attribute variance
#' 
#' 
#' @usage rs1(Tl, offset=0,nu1=NULL,s10=NULL)
#' @param Tl a list of V X P normal attribute matrix
#' @param nu1 prior degrees of freedom 
#' @param s10 prior estimate of s1
#' @return a new value of s1
#' @author Selena Wang
#' @export rs1
rs1 <-function(Tl,offset = offset, nu1=NULL,s10=NULL)
  { 
  
  tmp=sapply(1:length(Tl), function(x) Tl[[x]]-offset[[x]], simplify = FALSE)
  tmp.1=sum(sapply(tmp, function(x) sum(as.numeric(x)^2,na.rm = TRUE), simplify = TRUE))
  # 
  # tmp=0
  # for (i in 1:length(Tl)){
  #   R=Tl[[i]]-offset[[i]]
  #   #Es=rbind(Es,matrix(as.numeric(R)))
  #   tmp=tmp+sum(as.numeric(R)^2,na.rm = TRUE)
  #   
  #   
  # }  
  
    N <- length(Tl)
    P <- ncol(Tl[[1]])
    V <- nrow(Tl[[1]])
    if(is.null(nu1)){ nu1<-1 } 
    if(is.null(s10)){ s10<-1 } 

    1/rgamma(1, (N*V*P+s10)/2 , (tmp.1+nu1*s10)/2 )
  }
