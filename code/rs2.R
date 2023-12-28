#' Gibbs update for connectivity variance
#' 
#' Gibbs update for connectivity variance
#' 
#' 
#' @usage rs2(Fl, offset=0,nu2=NULL,s20=NULL)
#' @param Fl a list of V X V normal connectivity matrix
#' @param nu2 prior degrees of freedom 
#' @param s20 prior estimate of s2
#' @return a new value of s2
#' @author Selena Wang
#' @export rs2
rs2 <-function(Fl,offset = offset, nu2=NULL,s20=NULL)
  { 

  tmp=sapply(1:length(Fl), function(x) Fl[[x]]-offset[[x]], simplify = FALSE)
  tmp.1=sum(sapply(tmp, function(x) sum(x[upper.tri(x, diag = FALSE)]^2, na.rm = TRUE), simplify = TRUE))
  #sapply(tmp.1, sum, na.rm=TRUE )
  
  # for (i in 1:length(Fl)){
  #   R=Fl[[i]]-offset[[i]]
  #   
  #   tmp=tmp+ sum(R[upper.tri(R, diag = FALSE)]^2, na.rm = TRUE)
  #   #Es=rbind(Es,matrix(R[upper.tri(R, diag = FALSE)]))
  #   
  #   
  # }  
  
    N <- length(Fl)
    V <- nrow(Fl[[1]])
    if(is.null(nu2)){ nu2<-1 } 
    if(is.null(s20)){ s20<-1 } 

    1/rgamma(1, (N*V*(V-1)/2+s20)/2 , (tmp.1+nu2*s20)/2 )
  }
