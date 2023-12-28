#' Conditional simulation of intercept and regression coefficients
#' 
#' Simulates from the joint full conditional distribution of (beta,intercept)
#' in a brain connectivity model
#' 
#' @param Fl a list of V X V normal connectivity matrix
#' @param W N x Q covariate matrix
#' @param s2 dyadic variance
#' @param offset a V by V matrix, UU 
#' @param S0 prior precision matrix for regression parameters
#' @param beta0 prior mean vector for regression parameters 
#' @param ivA prior inverse variance for the intercept parameters 
#' 
#' @return \item{beta}{regression coefficients} \item{a}{subject-specific intercept}
#' @author Selena Wang
#' @export rbeta_a_fc
rbeta_a_fc_per<-
  function(Fl,W=NULL,s2=1,offset=offset,ivA=NULL,beta0=NULL,S0=NULL)
  {
    V<-nrow(Fl[[1]]) ; N <-length(Fl) 
    
    if(is.null(ivA))
    { 
      ivA <- matrix(1, nrow = 1, ncol = 1)
      
    } 
    

    
    one_vector <- matrix(1,nrow = V*(V-1)/2, ncol = 1)
    #G <- t(kronecker(diag(1, nrow = N), t(one_vector)))
    
    to<-as.numeric(sqrt(solve(s2)))
    a=matrix(NA, nrow = N, ncol = 1)
    for (i in 1:length(Fl)){
      R=Fl[[i]]-offset
      r=matrix(R[upper.tri(R, diag = FALSE)])
      r[is.na(r)]=0
      rs<- to*r
      
      
      m<-rs 
      iQA = solve( to^2*t(one_vector) %*% one_vector + ivA) 
      
      lA <- to*t(one_vector)%*%m
      a[i,] <- iQA%*%lA + t(chol(iQA))%*%rnorm(1) 

      
    }
    

      return(list(beta=NULL,a=a ) )
    
    
    
    
  }
