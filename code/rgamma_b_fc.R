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
rgamma_b_fc<-
  function(Tl,H=NULL,s1=1,offset=offset,ivB=NULL,gamma0=NULL,S1=NULL)
  {
    V<-nrow(Tl[[1]]) ; N <-length(Tl) ; P<-ncol(Tl[[1]])
    
    
    if(is.null(ivB))
    { 
      ivB <- diag(1,N)
      
    } 
    
    
    ### set priors 
    if(!is.null(H) &is.null(S1))
    { 
      # g-prior plus small ridge in case XX is singular
      S1 <- diag(1,ncol(H))
    } 
    
    if( ! is.null(H) &is.null(gamma0))
    { 
      gamma0 <- matrix(rep(0,ncol(H)), nrow = ncol(H),ncol=1)
      
    } 
    
    one_vector <- matrix(1,nrow = V*P, ncol = 1)
    G <- t(kronecker(diag(1, nrow = N), t(one_vector)))
    
    to<-as.numeric(sqrt(solve(s1)))
    r=matrix(NA, nrow=0, ncol=1)
    for (i in 1:length(Tl)){
      R=Tl[[i]]-offset
      r=rbind(r,matrix(as.numeric(R)))
      
      
    }
    r[is.na(r)]=0
    
    rs<- to*r
    
    iQB = solve( to^2*t(G) %*% G + ivB) 
    
    
    
    ### update gamma
    
    
    if(!is.null(H)){
      S <- G %*% iQB %*% t(G)
      
      
      Hg <- kronecker(H,one_vector)
      iQG <- solve(to^2*t(Hg) %*% S %*% Hg + S1)
      lG <- (to^2*t(Hg) %*% S %*% r + S1 %*% gamma0)
      gamma <- iQG%*%lG + t(chol(iQG))%*%rnorm(ncol(H)) 
      
      m<-rs - to*kronecker(H %*% gamma, one_vector)
      lB <- to*t(G)%*%m
      b <- iQB%*%lB + t(chol(iQB))%*%rnorm(N) 
      
      
      return(list(gamma=gamma,b=b ) )
    }else{
      m<-rs 
      lB <- to*t(G)%*%m
      b <- iQB%*%lB + t(chol(iQB))%*%rnorm(N) 
      
      
      
      return(list(gamma=NULL,b=b ) )
    }
    
    
    
  }
