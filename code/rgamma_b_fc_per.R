#' Conditional simulation of intercept and regression coefficients
#' 
#' Simulates from the joint full conditional distribution of (gamma,intercept)
#' in a brain connectivity model
#' 
#' @param Tl a list of V X P normal attribute matrix
#' @param H N x Q1 covariate matrix for the attributes, the first column should be 1 indicating the global intercept
#' @param s1 attribute variance
#' @param offset a list of the same dimension as Tl. 
#' @param S1 prior precision matrix for regression parameters
#' @param gamma0 prior mean vector for regression parameters 
#' @param ivB prior inverse variance for the intercept parameters 
#' 
#' @return \item{gamma}{regression coefficients} \item{b}{subject-specific intercept for attributes}
#' @author Selena Wang
#' @export rgamma_b_fc
rgamma_b_fc_per<-
  function(Tl,H=NULL,s1=1,offset=offset,ivB=NULL,gamma0=NULL,S1=NULL)
  {
    V<-nrow(Tl[[1]]) ; N <-length(Tl) ; P<-ncol(Tl[[1]])
    
    
    if(is.null(ivB))
    { 
      ivB <- matrix(1, nrow = 1, ncol = 1)
      
    } 
    


    one_vector <- matrix(1,nrow = V*P, ncol = 1)
    to<-as.numeric(sqrt(solve(s1)))
    b=matrix(NA, nrow = N, ncol = 1)
    
    for (i in 1:length(Tl)){
      R=Tl[[i]]-offset
      r=matrix(as.numeric(R))
      r[is.na(r)]=0
      
      rs<- to*r
      iQB = solve( to^2*t(one_vector) %*% one_vector + ivB) 
      
      m<-rs 
      lB <- to*t(one_vector)%*%m
      b[i,] <- iQB%*%lB + t(chol(iQB))%*%rnorm(1) 
      
      
      
    }

    
    
    

      
      return(list(gamma=NULL,b=b ) )
    
    
    
    
  }
